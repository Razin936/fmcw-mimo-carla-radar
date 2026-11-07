"""
carla_radar_sender.py
=====================================================================
 CARLA  ->  MATLAB  |  Synthetic 77 GHz FMCW MIMO radar data cube
=====================================================================

 What this script does
 ---------------------
  1. Builds a CARLA traffic scene (ego vehicle with autopilot, a lead
     vehicle, NPC vehicles and AI pedestrians) and attaches CARLA's
     ray-cast radar sensor to the ego vehicle.
  2. On every radar tick, filters the CARLA detections and synthesises
     a complex IQ radar data cube [numVirtual x numADC x numChirps]
     from them (FMCW beat signal, Doppler phase and per-element
     spatial phase).
  3. Streams the raw detections plus the cube to MATLAB over UDP,
     fragmented into magic-tagged packets.

 Threads
 -------
  The CARLA sensor callback only copies the detections and returns.
  A sender thread builds and transmits the cube (keeping only the
  newest frame if it falls behind), and a separate logger thread
  writes radar_log.csv, so neither blocks the simulator.

 Array model (2 TX x 4 RX -> 8 virtual elements)
 -----------------------------------------------
  RX spacing lambda/2, TX spacing 4*lambda/2, giving a uniform
  8-element virtual ULA at [0, d, ..., 7d], d = lambda/2.
  NOTE: every virtual channel is synthesised on every chirp (ideal
  MIMO). TDM transmit scheduling and its time offset are NOT modelled,
  and the cube carries no elevation phase. See README "Known
  limitations".

 UDP protocol (all integers little-endian)
 -----------------------------------------
  META  packet (12 bytes):
    uint32 MAGIC_META = 0xCAFEBEEF | uint32 num_chunks | uint32 total_bytes
  CHUNK packet (16 + chunk_len bytes):
    uint32 MAGIC_CHUNK = 0xDEADF00D | uint32 frame_id |
    uint32 chunk_idx (0-based)      | uint32 chunk_len | payload

  Assembled frame payload:
    bytes  0-3   uint32  numVirtual
    bytes  4-7   uint32  numADC
    bytes  8-11  uint32  numChirps
    bytes 12-19  double  timestamp (CARLA simulation time, s)
    bytes 20-23  uint32  numDet
    bytes 24-    float32 [numDet x 4]  raw detections [vel, az, alt, depth]
    then         complex64 cube [numVirtual x numADC x numChirps] (C order)

 Simulation timing
 -----------------
  With USE_SYNC_MODE the script switches CARLA to synchronous mode with
  a fixed 50 ms step and paces the steps to real time. The scene then
  moves evenly and every radar frame sees the same time step. The
  original world settings are restored when the script stops.

 How to run
 ----------
  1. Start the CARLA simulator (server on localhost:2000)
  2. python carla_radar_sender.py
  3. Run matlab/carla_radar_full_live.m in MATLAB
"""

import carla
import csv
import queue
import random
import socket
import struct
import threading
import time

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

# ---- Simulation timing ----------------------------------------------
USE_SYNC_MODE    = True    # False = CARLA's default asynchronous mode
SYNC_FIXED_DELTA = 0.05    # s per simulation step (radar fires every 2nd step)

# ---- FMCW MIMO waveform parameters  (77 GHz) ----------------------
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
D_TX        = NUM_RX * D_RX    # -> uniform lambda/2 virtual ULA

NOISE_SIGMA = 0.05             # receiver noise added in Python (MATLAB adds AWGN at SNR_dB on top)

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

# ---- Protocol ------------------------------------------------------
MAGIC_META  = 0xCAFEBEEF
MAGIC_CHUNK = 0xDEADF00D
CHUNK_SIZE  = 60000
HEADER_FMT  = "<IIIdI"                     # 24 bytes, no padding
HEADER_BYTES = struct.calcsize(HEADER_FMT)

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

# ---- Pre-computed axes and array geometry -------------------------
FS        = NUM_ADC / T_CHIRP
T_AXIS    = np.arange(NUM_ADC) / FS                 # fast time  [NUM_ADC]
CHIRP_IDX = np.arange(NUM_CHIRPS)                   # slow time  [NUM_CHIRPS]

RX_POS   = np.arange(NUM_RX) * D_RX                 # [4]
TX_POS   = np.array([0.0, D_TX])                    # [2]
VIRT_POS = np.array([tx + rx for tx in TX_POS for rx in RX_POS])   # [8]


# -------------------------------------------------------------------
# Signal model
# -------------------------------------------------------------------
def detections_to_array(radar_data):
    """Copy a CARLA RadarMeasurement into a float32 array [N x 4]
    with columns [velocity, azimuth_rad, altitude_rad, depth_m]."""
    return np.array(
        [[d.velocity, d.azimuth, d.altitude, d.depth] for d in radar_data],
        dtype=np.float32,
    ).reshape(-1, 4)


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


def build_rdc(dets, rng=np.random):
    """
    Synthesise the FMCW MIMO radar data cube from CARLA detections.

    For each kept detection i (range R, radial velocity v, azimuth az):
      f_r   = 2 R slope / c          range beat frequency
      f_d   = 2 v f_c / c            Doppler frequency
      f_b   = f_r + f_d              beat frequency within a chirp
      phi_d = 2 pi f_d m T_chirp     slow-time Doppler phase, chirp m
      phi_s = 2 pi/lambda p_q sin(az)  spatial phase, virtual element q

      cube[q, n, m] = sum_i a_i exp(j(2 pi f_b t_n + phi_d + phi_s)) + noise

    The sum is evaluated as one matrix product instead of a Python
    loop over detections and elements.

    Returns complex64 array [NUM_VIRTUAL, NUM_ADC, NUM_CHIRPS].
    """
    dets = np.asarray(dets, dtype=np.float64).reshape(-1, 4)
    dets = dets[filter_detections(dets)]

    if len(dets):
        v, az, R = dets[:, 0], dets[:, 1], dets[:, 3]
        amp = amplitude_weights(v, R)

        f_d = 2.0 * v * FC / C
        f_b = 2.0 * R * SLOPE / C + f_d

        spatial = amp[:, None] * np.exp(
            1j * (2.0 * np.pi / LAMBDA) * VIRT_POS[None, :] * np.sin(az)[:, None])  # [K, V]
        fast = np.exp(1j * 2.0 * np.pi * f_b[:, None] * T_AXIS[None, :])          # [K, A]
        slow = np.exp(1j * 2.0 * np.pi * f_d[:, None] * CHIRP_IDX[None, :] * T_CHIRP)  # [K, C]

        sv_fast = (spatial[:, :, None] * fast[:, None, :]).reshape(len(dets), -1)  # [K, V*A]
        rdc = (sv_fast.T @ slow).reshape(NUM_VIRTUAL, NUM_ADC, NUM_CHIRPS)
    else:
        rdc = np.zeros((NUM_VIRTUAL, NUM_ADC, NUM_CHIRPS), dtype=np.complex128)

    noise = (NOISE_SIGMA / np.sqrt(2.0)) * (
        rng.standard_normal(rdc.shape) + 1j * rng.standard_normal(rdc.shape))
    return (rdc + noise).astype(np.complex64)


# -------------------------------------------------------------------
# UDP transport
# -------------------------------------------------------------------
def build_payload(dets, timestamp):
    """Header + raw detections + IQ cube, as one bytes object."""
    dets = np.asarray(dets, dtype=np.float32).reshape(-1, 4)
    header = struct.pack(HEADER_FMT, NUM_VIRTUAL, NUM_ADC, NUM_CHIRPS,
                         float(timestamp), len(dets))
    det_bytes = dets.astype("<f4").tobytes()
    cube_bytes = build_rdc(dets).astype("<c8").tobytes()
    return header + det_bytes + cube_bytes


def send_frame(sock, dets, timestamp, frame_id):
    """Send one frame as a META packet followed by CHUNK packets."""
    payload = build_payload(dets, timestamp)
    total = len(payload)
    num_chunks = (total + CHUNK_SIZE - 1) // CHUNK_SIZE

    sock.sendto(struct.pack("<III", MAGIC_META, num_chunks, total),
                (MATLAB_IP, MATLAB_PORT))
    time.sleep(0.001)

    fid = frame_id & 0xFFFFFFFF
    for i in range(num_chunks):
        chunk = payload[i * CHUNK_SIZE:(i + 1) * CHUNK_SIZE]
        hdr = struct.pack("<IIII", MAGIC_CHUNK, fid, i, len(chunk))
        sock.sendto(hdr + chunk, (MATLAB_IP, MATLAB_PORT))
        time.sleep(0.0002)   # pacing so the loopback receive buffer is not overrun


class RadarStreamer:
    """Decouples the CARLA sensor callback from cube synthesis, UDP
    dispatch and CSV logging.

    - on_radar() runs in CARLA's sensor thread: it copies the detections
      and returns immediately.
    - The sender thread keeps only the newest frame (older queued frames
      are dropped and counted) and sends it.
    - The logger thread writes every frame to the CSV file.
    """

    def __init__(self, sock, csv_writer):
        self.sock = sock
        self.csv_writer = csv_writer
        self.send_q = queue.Queue(maxsize=1)
        self.log_q = queue.Queue()
        self.stop_event = threading.Event()
        self.frame_id = 0
        self.dropped = 0
        self.sent = 0
        self._threads = [
            threading.Thread(target=self._send_loop, name="radar-sender", daemon=True),
            threading.Thread(target=self._log_loop, name="radar-logger", daemon=True),
        ]

    def start(self):
        for t in self._threads:
            t.start()

    def on_radar(self, radar_data):
        item = (float(radar_data.timestamp), detections_to_array(radar_data))
        try:
            self.send_q.put_nowait(item)
        except queue.Full:
            try:
                self.send_q.get_nowait()       # discard the older, unsent frame
                self.dropped += 1
            except queue.Empty:
                pass
            try:
                self.send_q.put_nowait(item)
            except queue.Full:
                self.dropped += 1
        self.log_q.put(item)

    def _send_loop(self):
        while not self.stop_event.is_set():
            try:
                ts, dets = self.send_q.get(timeout=0.2)
            except queue.Empty:
                continue
            try:
                send_frame(self.sock, dets, ts, self.frame_id)
                self.frame_id += 1
                self.sent += 1
            except OSError as exc:
                print(f"[Sender] UDP error: {exc}")

    def _log_loop(self):
        while not (self.stop_event.is_set() and self.log_q.empty()):
            try:
                ts, dets = self.log_q.get(timeout=0.2)
            except queue.Empty:
                continue
            if self.csv_writer is not None:
                for v, az, alt, dep in dets:
                    self.csv_writer.writerow([ts, len(dets), v, az, alt, dep])
            if len(dets):
                v, dep = dets[:, 0], dets[:, 3]
                print(f"[Radar] t={ts:.2f}s | pts={len(dets)} | nearest={dep.min():.2f}m | "
                      f"vel=[{v.min():.2f}, {v.max():.2f}] | mean={v.mean():.2f} m/s | "
                      f"sent={self.sent} dropped={self.dropped}")
            else:
                print(f"[Radar] t={ts:.2f}s | pts=0 | sent={self.sent} dropped={self.dropped}")

    def stop(self):
        self.stop_event.set()
        for t in self._threads:
            t.join(timeout=2.0)


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


def spawn_lead_vehicle(world, vehicle_bps, tm, ego_vehicle,
                       distance_m=25.0, min_clearance_m=10.0, max_tries=10):
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
    if USE_SYNC_MODE:
        world.tick()           # walker controllers need one tick before start()
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
    radar = None
    streamer = None
    original_settings = None

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

    tm.set_global_distance_to_leading_vehicle(6.0)
    tm.global_percentage_speed_difference(30.0)
    try:
        tm.set_hybrid_physics_mode(True)
        tm.set_hybrid_physics_radius(70.0)
    except Exception:
        pass

    try:
        if USE_SYNC_MODE:
            original_settings = world.get_settings()
            settings = world.get_settings()
            settings.synchronous_mode = True
            settings.fixed_delta_seconds = SYNC_FIXED_DELTA
            world.apply_settings(settings)
            tm.set_synchronous_mode(True)
            print(f"[SYNC] synchronous mode, fixed step {SYNC_FIXED_DELTA} s")
        else:
            tm.set_synchronous_mode(False)

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
        if USE_SYNC_MODE:
            world.tick()       # register the ego vehicle before spawning around it

        # Lead vehicle
        if SPAWN_LEAD_VEHICLE:
            lead = spawn_lead_vehicle(
                world, vehicle_bps, tm, ego,
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
            try:
                # carla.World method; must be set before the walkers are spawned
                world.set_pedestrians_cross_factor(WALKER_CROSSING_PROB)
            except Exception:
                pass
            walkers, controllers = spawn_walkers(
                world, blueprint_library, NUM_WALKERS)
            actors_to_destroy.extend(controllers)
            actors_to_destroy.extend(walkers)

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

        streamer = RadarStreamer(sock, csv_writer)
        streamer.start()
        radar.listen(streamer.on_radar)

        spectator = world.get_spectator()
        next_t = time.perf_counter()
        while True:
            if USE_SYNC_MODE:
                world.tick()
                next_t += SYNC_FIXED_DELTA            # pace steps to real time
                delay = next_t - time.perf_counter()
                if delay > 0:
                    time.sleep(delay)
                else:
                    next_t = time.perf_counter()      # slower than real time: no catch-up burst
            else:
                time.sleep(0.05)

            if FOLLOW_CAMERA and ego.is_alive:
                tf  = ego.get_transform()
                fwd = tf.get_forward_vector()
                cam_loc = tf.location - carla.Location(
                    x=8.0 * fwd.x, y=8.0 * fwd.y, z=-4.0)   # 8 m behind, 4 m above
                spectator.set_transform(carla.Transform(
                    cam_loc,
                    carla.Rotation(pitch=-15.0, yaw=tf.rotation.yaw, roll=0.0)))

    except KeyboardInterrupt:
        print("\nStopping...")

    finally:
        if radar is not None:
            try:
                radar.stop()
            except Exception:
                pass
        if streamer is not None:
            streamer.stop()
        if original_settings is not None:
            try:
                tm.set_synchronous_mode(False)
                world.apply_settings(original_settings)   # never leave CARLA frozen in sync mode
                print("[SYNC] world settings restored")
            except Exception:
                pass
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
