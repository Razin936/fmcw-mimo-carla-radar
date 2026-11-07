"""
carla_radar_sender.py
=====================================================================
 CARLA  ->  MATLAB  |  True FMCW MIMO Radar Simulation
=====================================================================

 MIMO layout  (TI AWR1843-style):
   numTX = 2  (spacing = numRX * lambda/2)
   numRX = 4  (half-wavelength spacing)
   numVirtual = 8  uniform lambda/2 virtual ULA

 UDP fragmented protocol  (MAGIC-identified, robust):
 --------------------------------------------------------
  META packet  (12 bytes):
    uint32  MAGIC_META  = 0xCAFEBEEF
    uint32  num_chunks
    uint32  total_bytes

  CHUNK packets  (16 + payload bytes each):
    uint32  MAGIC_CHUNK = 0xDEADF00D
    uint32  frame_id
    uint32  chunk_idx   (0-based)
    uint32  chunk_len
    uint8   payload[chunk_len]

  Every datagram carries its own magic number so MATLAB can
  always identify packet type regardless of arrival order.

 Full frame payload layout  (assembled from chunks):
 --------------------------------------------------------
  Bytes  0-3    uint32  numVirtual
  Bytes  4-7    uint32  numADC
  Bytes  8-11   uint32  numChirps
  Bytes 12-19   double  timestamp
  Bytes 20-23   uint32  numDet
  Bytes 24-24+numDet*16   float32[numDet][4]  [vel,az,alt,dep]
  Remaining               complex64 IQ cube   [numVirt x numADC x numChirps]

HOW TO RUN
----------
  1. Start CARLA simulator
  2. python carla_radar_sender.py
  3. Run carla_radar_full_live.m in MATLAB
"""

import carla
import random
import time
import socket
import struct
import csv
import numpy as np

# ---- UDP destination -----------------------------------------------
MATLAB_IP   = "127.0.0.1"
MATLAB_PORT = 5005

# ---- CARLA radar sensor parameters --------------------------------
RADAR_RANGE             = 60.0
RADAR_HORIZONTAL_FOV    = 30.0
RADAR_VERTICAL_FOV      = 10.0
RADAR_POINTS_PER_SECOND = 1500
RADAR_SENSOR_TICK       = 0.1

# ---- FMCW MIMO waveform parameters  (77 GHz TI AWR1843-style) ----
C           = 3e8
FC          = 77e9
BW          = 150e6
T_CHIRP     = 10e-6
SLOPE       = BW / T_CHIRP
LAMBDA      = C / FC

NUM_TX      = 2
NUM_RX      = 4
NUM_VIRTUAL = NUM_TX * NUM_RX   # 8
NUM_ADC     = 256
NUM_CHIRPS  = 128

D_RX        = LAMBDA / 2
D_TX        = NUM_RX * D_RX    # TX spacing = 4 * lambda/2 -> uniform virtual ULA

NOISE_SIGMA = 0.05

# ---- Detection filter applied before cube synthesis ---------------
# Same gates as the MATLAB display filter, so the cube and the plots
# are built from the same set of detections.
MIN_DEPTH   = 4.0     # m   removes near-field self-returns
MAX_DEPTH   = 60.0    # m
MIN_ABS_VEL = 0.5     # m/s removes static returns (|v| <= MIN_ABS_VEL)
MAX_ABS_VEL = 40.0    # m/s

# ---- Heuristic reflectivity (CARLA provides no RCS) ---------------
SLOW_SPEED_LIMIT = 2.0   # m/s  targets slower than this ...
SLOW_RCS_WEIGHT  = 0.1   #      ... get this amplitude weight (pedestrian-like)

# ---- Protocol magic numbers ---------------------------------------
MAGIC_META  = 0xCAFEBEEF
MAGIC_CHUNK = 0xDEADF00D
CHUNK_SIZE  = 60000

# ---- Scene settings -----------------------------------------------
FOLLOW_CAMERA        = True
NUM_NPC_VEHICLES     = 12
SPAWN_LEAD_VEHICLE   = True
LEAD_DISTANCE_METERS = 25.0
LEAD_MIN_CLEARANCE   = 10.0
SPAWN_WALKERS        = True
NUM_WALKERS          = 20
WALKER_RUNNING_PROB  = 0.10
WALKER_CROSSING_PROB = 0.20
WALKER_SPAWN_TRIES   = 120

# ---- Pre-compute virtual element positions ------------------------
FS     = NUM_ADC / T_CHIRP
DT     = 1.0 / FS
T_AXIS = np.arange(NUM_ADC) * DT                   # [NUM_ADC]

RX_POS   = np.arange(NUM_RX) * D_RX                # [4]
TX_POS   = np.array([0.0, D_TX])                   # [2]
VIRT_POS = np.array(
    [tx + rx for tx in TX_POS for rx in RX_POS]
)                                                   # [8] uniform lambda/2 ULA

frame_id_counter = 0


# -------------------------------------------------------------------
def filter_detections(dets):
    """Boolean mask of detections used for the cube: finite,
    forward-facing, MIN_DEPTH <= R <= MAX_DEPTH and
    MIN_ABS_VEL < |v| <= MAX_ABS_VEL."""
    v, az, R = dets[:, 0], dets[:, 1], dets[:, 3]
    finite = np.isfinite(dets).all(axis=1)
    with np.errstate(invalid="ignore"):
        forward = R * np.cos(az) >= 0.0
        return (finite & forward
                & (R >= MIN_DEPTH) & (R <= MAX_DEPTH)
                & (np.abs(v) > MIN_ABS_VEL) & (np.abs(v) <= MAX_ABS_VEL))


def amplitude_weights(v, R):
    """Voltage amplitude per detection: heuristic RCS weight times the
    square root of an R^-4 power factor normalised to the nearest target."""
    rcs = np.where(np.abs(v) < SLOW_SPEED_LIMIT, SLOW_RCS_WEIGHT, 1.0)
    power = R ** -4.0
    power = power / power.max()
    return rcs * np.sqrt(power)


def build_rdc(detections_raw):
    """
    Synthesise a true FMCW MIMO Radar Data Cube.

    For each CARLA detection (R, v, az_rad):
      fr          = 2*R*slope/c              range beat frequency
      fd          = 2*v*fc/c                 Doppler frequency
      fbeat       = fr + fd                  total beat frequency
      phi_spatial = 2pi * d_v(n) * sin(az) / lambda   per virtual element
      per-chirp Doppler phase = 2pi * fd * chirp_idx * T_CHIRP

    Only detections that pass filter_detections() are used, each scaled
    by its amplitude_weights() value.

    Returns
    -------
    RDC : np.ndarray  complex64  [NUM_VIRTUAL, NUM_ADC, NUM_CHIRPS]
    """
    RDC = np.zeros((NUM_VIRTUAL, NUM_ADC, NUM_CHIRPS), dtype=np.complex64)

    dets = np.array([[d.velocity, d.azimuth, d.altitude, d.depth]
                     for d in detections_raw], dtype=np.float64).reshape(-1, 4)
    dets = dets[filter_detections(dets)]
    amps = amplitude_weights(dets[:, 0], dets[:, 3]) if len(dets) else []

    for (v, az, _alt, R), amp in zip(dets, amps):
        fr    = 2.0 * R * SLOPE / C
        fd    = 2.0 * v * FC    / C
        fbeat = fr + fd

        # Spatial phase per virtual element  [NUM_VIRTUAL]
        phi_spatial = (2.0 * np.pi / LAMBDA) * VIRT_POS * np.sin(az)

        # Fast-time envelope  [NUM_ADC]
        beat_time = np.exp(1j * 2.0 * np.pi * fbeat * T_AXIS)

        # Chirp loop
        chirp_idx_arr = np.arange(NUM_CHIRPS)
        phi_dop_arr   = 2.0 * np.pi * fd * chirp_idx_arr * T_CHIRP  # [NUM_CHIRPS]

        for v_idx in range(NUM_VIRTUAL):
            # [NUM_ADC, 1] x [1, NUM_CHIRPS] -> [NUM_ADC, NUM_CHIRPS]
            sig = amp * beat_time[:, np.newaxis] * np.exp(
                1j * (phi_dop_arr[np.newaxis, :] + phi_spatial[v_idx])
            )
            RDC[v_idx] += sig.astype(np.complex64)

    # Add complex AWGN
    noise = (NOISE_SIGMA / np.sqrt(2.0)) * (
        np.random.randn(*RDC.shape) + 1j * np.random.randn(*RDC.shape)
    )
    RDC += noise.astype(np.complex64)

    return RDC


def send_frame(sock, detections_raw, timestamp):
    """
    Build full frame payload and send as magic-identified fragments.
    """
    global frame_id_counter

    num_det = len(detections_raw)

    # Raw detection bytes  [vel, az, alt, dep]  float32 x4
    det_bytes = bytearray()
    for det in detections_raw:
        det_bytes += struct.pack("<ffff",
            float(det.velocity), float(det.azimuth),
            float(det.altitude), float(det.depth))

    # Build IQ cube
    RDC = build_rdc(detections_raw)   # [NUM_VIRTUAL, NUM_ADC, NUM_CHIRPS]
    rdc_bytes = RDC.astype(np.complex64).tobytes()

    # Fixed header  (28 bytes)
    header = struct.pack("<IIIdI",
        NUM_VIRTUAL, NUM_ADC, NUM_CHIRPS, timestamp, num_det)

    payload = bytes(header) + bytes(det_bytes) + rdc_bytes
    total   = len(payload)

    num_chunks = (total + CHUNK_SIZE - 1) // CHUNK_SIZE

    # Send META packet
    meta_pkt = struct.pack("<III", MAGIC_META, num_chunks, total)
    sock.sendto(meta_pkt, (MATLAB_IP, MATLAB_PORT))
    time.sleep(0.001)

    # Send CHUNK packets
    fid = frame_id_counter & 0xFFFFFFFF
    for i in range(num_chunks):
        chunk_data = payload[i * CHUNK_SIZE : (i+1) * CHUNK_SIZE]
        c_len      = len(chunk_data)
        chunk_hdr  = struct.pack("<IIII", MAGIC_CHUNK, fid, i, c_len)
        sock.sendto(chunk_hdr + chunk_data, (MATLAB_IP, MATLAB_PORT))
        time.sleep(0.0002)

    frame_id_counter += 1


def radar_callback(radar_data, sock, csv_writer=None):
    detections = list(radar_data)
    send_frame(sock, detections, radar_data.timestamp)

    if csv_writer is not None:
        for det in detections:
            csv_writer.writerow([
                radar_data.timestamp, len(detections),
                float(det.velocity), float(det.azimuth),
                float(det.altitude), float(det.depth)])

    if detections:
        nearest  = min(d.depth    for d in detections)
        min_vel  = min(d.velocity for d in detections)
        max_vel  = max(d.velocity for d in detections)
        mean_vel = sum(d.velocity for d in detections) / len(detections)
        print(f"[Radar] t={radar_data.timestamp:.2f}s | "
              f"pts={len(detections)} | nearest={nearest:.2f}m | "
              f"vel=[{min_vel:.2f}, {max_vel:.2f}] | mean={mean_vel:.2f} m/s")
    else:
        print(f"[Radar] t={radar_data.timestamp:.2f}s | pts=0")


# ---- Scene helpers -------------------------------------------------

def is_spawn_far_enough(candidate_spawn, existing_actors, min_distance=18.0):
    loc = candidate_spawn.location
    for actor in existing_actors:
        try:
            if actor.is_alive and loc.distance(actor.get_location()) < min_distance:
                return False
        except Exception:
            pass
    return True


def choose_vehicle_blueprint(vehicle_bps):
    safe = [bp for bp in vehicle_bps
            if bp.has_attribute("number_of_wheels")
            and int(bp.get_attribute("number_of_wheels").as_int()) == 4]
    bp = random.choice(safe) if safe else random.choice(vehicle_bps)
    for attr in ("color", "driver_id"):
        if bp.has_attribute(attr):
            bp.set_attribute(attr,
                random.choice(bp.get_attribute(attr).recommended_values))
    if bp.has_attribute("role_name"):
        bp.set_attribute("role_name", "autopilot")
    return bp


def setup_vehicle_behavior(vehicle, tm, is_ego=False):
    vehicle.set_autopilot(True, tm.get_port())
    tm.auto_lane_change(vehicle, False)
    tm.distance_to_leading_vehicle(vehicle, 6.0)
    tm.vehicle_percentage_speed_difference(vehicle, 35.0 if is_ego else 30.0)
    tm.ignore_lights_percentage(vehicle, 0.0)
    tm.ignore_signs_percentage(vehicle, 0.0)
    tm.ignore_vehicles_percentage(vehicle, 0.0)
    tm.ignore_walkers_percentage(vehicle, 0.0)
    try:
        tm.random_left_lanechange_percentage(vehicle, 0.0)
        tm.random_right_lanechange_percentage(vehicle, 0.0)
    except Exception:
        pass


def spawn_lead_vehicle(world, blueprint_library, vehicle_bps, tm,
                       ego_vehicle, distance_m=25.0,
                       min_clearance_m=10.0, max_tries=10):
    carla_map = world.get_map()
    ego_loc   = ego_vehicle.get_transform().location
    wp = carla_map.get_waypoint(
        ego_loc, project_to_road=True, lane_type=carla.LaneType.Driving)
    if wp is None:
        return None
    step = max(min_clearance_m, 5.0)
    for i in range(max_tries):
        dist = distance_m + i * step
        nxt  = wp.next(dist)
        if not nxt:
            continue
        tf = nxt[0].transform
        tf.location.z += 0.3
        lead = world.try_spawn_actor(choose_vehicle_blueprint(vehicle_bps), tf)
        if lead is None:
            continue
        setup_vehicle_behavior(lead, tm, is_ego=False)
        print(f"Lead vehicle spawned {dist:.1f}m ahead: {lead.type_id}")
        return lead
    return None


def _rand_nav(world, tries):
    for _ in range(tries):
        loc = world.get_random_location_from_navigation()
        if loc is not None:
            return loc
    return None


def spawn_walkers(world, blueprint_library, count):
    walker_bps    = blueprint_library.filter("walker.pedestrian.*")
    controller_bp = blueprint_library.find("controller.ai.walker")
    walkers, controllers = [], []
    if not walker_bps:
        print("No walker blueprints found.")
        return walkers, controllers
    for _ in range(count):
        nav_loc = _rand_nav(world, WALKER_SPAWN_TRIES)
        if nav_loc is None:
            continue
        bp     = random.choice(walker_bps)
        walker = world.try_spawn_actor(bp, carla.Transform(nav_loc))
        if walker is None:
            continue
        ctrl = world.try_spawn_actor(
            controller_bp, carla.Transform(), attach_to=walker)
        if ctrl is None:
            walker.destroy()
            continue
        walkers.append(walker)
        controllers.append(ctrl)
    for walker, ctrl in zip(walkers, controllers):
        ctrl.start()
        ctrl.set_max_speed(3.5 if random.random() < WALKER_RUNNING_PROB else 1.4)
        dest = _rand_nav(world, WALKER_SPAWN_TRIES) or walker.get_location()
        ctrl.go_to_location(dest)
    print(f"Spawned {len(walkers)} walkers.")
    return walkers, controllers


# ---- Main ----------------------------------------------------------

def main():
    actors_to_destroy = []

    sock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
    sock.setsockopt(socket.SOL_SOCKET, socket.SO_SNDBUF, 4 * 1024 * 1024)

    log_file   = open("radar_log.csv", "w", newline="")
    csv_writer = csv.writer(log_file)
    csv_writer.writerow([
        "timestamp", "num_points", "velocity", "azimuth", "altitude", "depth"])

    client = carla.Client("localhost", 2000)
    client.set_timeout(20.0)

    world             = client.get_world()
    blueprint_library = world.get_blueprint_library()
    tm                = client.get_trafficmanager(8000)

    tm.set_synchronous_mode(False)
    tm.set_global_distance_to_leading_vehicle(6.0)
    tm.global_percentage_speed_difference(30.0)
    try:
        tm.set_hybrid_physics_mode(True)
        tm.set_hybrid_physics_radius(70.0)
    except Exception:
        pass

    try:
        vehicle_bps  = blueprint_library.filter("vehicle.*")
        spawn_points = world.get_map().get_spawn_points()
        if not spawn_points:
            print("No spawn points found.")
            return
        random.shuffle(spawn_points)

        # Ego vehicle
        ego = world.try_spawn_actor(
            choose_vehicle_blueprint(vehicle_bps), spawn_points.pop(0))
        if ego is None:
            print("Could not spawn ego vehicle.")
            return
        actors_to_destroy.append(ego)
        print(f"Ego spawned: {ego.type_id}")
        setup_vehicle_behavior(ego, tm, is_ego=True)

        # Lead vehicle
        if SPAWN_LEAD_VEHICLE:
            lead = spawn_lead_vehicle(
                world, blueprint_library, vehicle_bps, tm, ego,
                distance_m=LEAD_DISTANCE_METERS,
                min_clearance_m=LEAD_MIN_CLEARANCE)
            if lead:
                actors_to_destroy.append(lead)

        # NPC vehicles
        npc_count = 0
        for sp in spawn_points:
            if npc_count >= NUM_NPC_VEHICLES:
                break
            if not is_spawn_far_enough(sp, actors_to_destroy, 18.0):
                continue
            npc = world.try_spawn_actor(
                choose_vehicle_blueprint(vehicle_bps), sp)
            if npc is None:
                continue
            actors_to_destroy.append(npc)
            setup_vehicle_behavior(npc, tm, is_ego=False)
            npc_count += 1
        print(f"Spawned {npc_count} NPC vehicles.")

        # Walkers
        if SPAWN_WALKERS:
            walkers, controllers = spawn_walkers(
                world, blueprint_library, NUM_WALKERS)
            actors_to_destroy.extend(controllers)
            actors_to_destroy.extend(walkers)
            try:
                tm.set_pedestrians_cross_factor(WALKER_CROSSING_PROB)
            except Exception:
                pass

        # Radar sensor
        radar_bp = blueprint_library.find("sensor.other.radar")
        radar_bp.set_attribute("horizontal_fov",    str(RADAR_HORIZONTAL_FOV))
        radar_bp.set_attribute("vertical_fov",      str(RADAR_VERTICAL_FOV))
        radar_bp.set_attribute("range",             str(RADAR_RANGE))
        radar_bp.set_attribute("points_per_second", str(RADAR_POINTS_PER_SECOND))
        radar_bp.set_attribute("sensor_tick",       str(RADAR_SENSOR_TICK))

        radar_tf = carla.Transform(
            carla.Location(x=2.2, z=1.0),
            carla.Rotation(pitch=0.0, yaw=0.0, roll=0.0))
        radar = world.spawn_actor(radar_bp, radar_tf, attach_to=ego)
        actors_to_destroy.append(radar)

        print("FMCW MIMO radar attached to ego vehicle.")
        print(f"  TX={NUM_TX}  RX={NUM_RX}  Virtual={NUM_VIRTUAL}")
        print(f"  ADC={NUM_ADC}  Chirps={NUM_CHIRPS}")
        print(f"  fc={FC/1e9:.0f}GHz  BW={BW/1e6:.0f}MHz  T={T_CHIRP*1e6:.0f}us")
        print("Streaming IQ data cube to MATLAB on UDP port 5005 ...")
        print("Press Ctrl+C to stop.")

        radar.listen(lambda data: radar_callback(data, sock, csv_writer))

        spectator = world.get_spectator()
        while True:
            if FOLLOW_CAMERA and ego.is_alive:
                tf  = ego.get_transform()
                fwd = tf.get_forward_vector()
                cam_loc = tf.location - carla.Location(
                    x=8.0*fwd.x, y=8.0*fwd.y, z=-4.0)
                cam_loc.z += 4.0
                spectator.set_transform(carla.Transform(
                    cam_loc,
                    carla.Rotation(pitch=-15.0, yaw=tf.rotation.yaw, roll=0.0)))
            time.sleep(0.05)

    except KeyboardInterrupt:
        print("\nStopping...")

    finally:
        print("Destroying actors...")
        for actor in actors_to_destroy:
            try:
                if actor.is_alive:
                    actor.destroy()
            except Exception:
                pass
        try:
            log_file.close()
        except Exception:
            pass
        sock.close()
        print("Done.")


if __name__ == "__main__":
    main()