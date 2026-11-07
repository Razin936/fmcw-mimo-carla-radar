"""
cut_in_demo.py
=====================================================================
 Scripted cut-in scene for the CARLA -> MATLAB radar pipeline
=====================================================================

 A short, repeatable scene for demonstrations and screen recordings.
 It uses the same radar sensor, cube synthesis and UDP stream as
 carla_radar_sender.py, so the MATLAB receiver runs unchanged.

 Timeline (all values can be changed below)
 ------------------------------------------
  1. Free road     The ego vehicle drives alone in its lane. A second
                   car follows in the left lane, behind the ego and
                   outside the radar field of view.
  2. Overtake      The second car speeds up and passes the ego.
  3. Cut-in        Close in front of the ego it changes into the ego lane.
  4. Braking       It slows down. The ego brakes with a simple distance
                   controller, so the radar sees a closing target
                   (negative relative velocity).
  5. Pull away     The car ahead speeds up again and the gap opens
                   (positive relative velocity).
  6. End           Both cars stop and the script exits.

 Detections with |v| <= MIN_ABS_VEL (0.5 m/s) are removed by the
 sender's filter. A car that drives at exactly the ego speed therefore
 disappears from the plots, which is why the scene ends with the car
 pulling away rather than with a steady gap.

 The scene runs on a straight two-lane stretch (Town04 by default),
 without other traffic or pedestrians. Both cars are driven by this
 script, not by the Traffic Manager, so every run looks the same.
 The phases are printed with the simulation time, which matches the
 t= value in the MATLAB frame log.

 How to run
 ----------
  1. Start the CARLA simulator (server on localhost:2000)
  2. Run matlab/carla_radar_full_live.m in MATLAB
  3. python cut_in_demo.py
  When you are done, return to the map of the normal scene with
     python cut_in_demo.py --load-map Town10HD_Opt
  (or restart CARLA) before running carla_radar_sender.py again.
"""

import argparse
import csv
import math
import socket
import time

import carla
import numpy as np

import carla_radar_sender as radar_io

# ---- Scene -----------------------------------------------------------
DEMO_TOWN        = "Town04"   # long straight multi-lane roads; None = keep the current map
EGO_SPEED        = 10.0       # m/s  ego cruise speed (36 km/h)
START_BEHIND_M   = 10.0       # m    the second car starts this far behind the ego
FREE_DRIVE_S     = 6.0        # s    free road before the overtake starts
OVERTAKE_SPEED   = 13.0       # m/s  speed of the second car while overtaking
CUT_IN_AT_M      = 8.0        # m    start the lane change this far ahead (centre to centre)
LANE_CHANGE_S    = 2.5        # s    duration of the lane change
CUT_IN_SPEED     = 5.0        # m/s  speed of the second car after the cut-in
CUT_IN_DECEL     = 3.0        # m/s^2 how quickly it slows down
BRAKE_S          = 4.0        # s    braking phase
PULL_AWAY_SPEED  = 13.0       # m/s  speed of the car ahead when it pulls away
PULL_AWAY_ACCEL  = 2.0        # m/s^2
PULL_AWAY_S      = 7.0        # s    pull-away phase
FOLLOW_GAP_M     = 15.0       # m    distance the ego keeps behind it (centre to centre)
CUT_IN_SIDE      = "left"     # lane the second car comes from: "left" or "right"
SHOW_LABELS      = True       # phase text above the second car in the CARLA window

EGO_BLUEPRINT    = "vehicle.tesla.model3"
CUT_IN_BLUEPRINT = "vehicle.audi.a2"

# ---- Camera (chase view behind the ego) ----------------------------
CAM_BEHIND_M = 10.0
CAM_ABOVE_M  = 5.0
CAM_PITCH    = -15.0

# ---- Route search ----------------------------------------------------
PATH_SPACING_M = 2.0          # waypoint spacing of the precomputed route
PATH_LENGTH_M  = 400.0        # route length needed for the whole scene
HEADING_LIMITS_DEG = (30.0, 90.0, 180.0)   # prefer straight roads, accept curves if needed


def wrap_deg(a):
    return (a + 180.0) % 360.0 - 180.0


def smoothstep(x):
    x = min(max(x, 0.0), 1.0)
    return x * x * (3.0 - 2.0 * x)


# ---- Route ----------------------------------------------------------

def side_lane(wp, side):
    lane = wp.get_left_lane() if side == "left" else wp.get_right_lane()
    if lane is None or lane.lane_type != carla.LaneType.Driving:
        return None
    if lane.lane_id * wp.lane_id <= 0:          # opposite direction
        return None
    return lane


def build_route(start_wp, side, max_heading_deg):
    """Ego-lane and side-lane points [N x 2] with matching indices, the
    road height [N], or None."""
    n = int(PATH_LENGTH_M / PATH_SPACING_M)
    yaw0 = start_wp.transform.rotation.yaw
    ego_pts, side_pts, z = [], [], []
    wp = start_wp
    for _ in range(n):
        if wp.is_junction:
            return None
        sw = side_lane(wp, side)
        if sw is None:
            return None
        if abs(wrap_deg(wp.transform.rotation.yaw - yaw0)) > max_heading_deg:
            return None
        ego_pts.append((wp.transform.location.x, wp.transform.location.y))
        side_pts.append((sw.transform.location.x, sw.transform.location.y))
        z.append(max(wp.transform.location.z, sw.transform.location.z))
        nxt = wp.next(PATH_SPACING_M)
        if not nxt:
            return None
        wp = nxt[0]
    return np.array(ego_pts), np.array(side_pts), np.array(z)


def find_route(carla_map, side, rng):
    candidates = [wp for wp in carla_map.generate_waypoints(10.0)
                  if wp.lane_type == carla.LaneType.Driving]
    rng.shuffle(candidates)
    for limit in HEADING_LIMITS_DEG:
        for wp in candidates[:3000]:
            route = build_route(wp, side, limit)
            if route is not None:
                return route
    return None


# ---- Vehicle control ---------------------------------------------

class RouteDriver:
    """Pure-pursuit steering and PI speed control along a route."""

    def __init__(self, vehicle, ego_pts, side_pts, dt):
        self.v = vehicle
        self.ego_pts, self.side_pts = ego_pts, side_pts
        self.dt = dt
        self.idx = 0
        self.lane_mix = 0.0               # 0 = side lane, 1 = ego lane
        self.i_term = 8.0                 # start near a cruising throttle
        try:
            pc = vehicle.get_physics_control()
            self.max_steer_deg = max(w.max_steer_angle for w in pc.wheels) or 70.0
        except Exception:
            self.max_steer_deg = 70.0
        self.wheelbase = 2.8

    def speed(self):
        vel = self.v.get_velocity()
        return math.hypot(vel.x, vel.y)

    def update_index(self):
        loc = self.v.get_transform().location
        lo, hi = max(self.idx - 5, 0), min(self.idx + 40, len(self.ego_pts))
        d = np.hypot(self.ego_pts[lo:hi, 0] - loc.x, self.ego_pts[lo:hi, 1] - loc.y)
        self.idx = lo + int(np.argmin(d))
        return self.idx

    def remaining(self):
        return len(self.ego_pts) - 1 - self.idx

    def target_point(self, lookahead_m):
        i = min(self.idx + int(math.ceil(lookahead_m / PATH_SPACING_M)),
                len(self.ego_pts) - 1)
        m = self.lane_mix
        return (1.0 - m) * self.side_pts[i] + m * self.ego_pts[i]

    def step(self, v_ref):
        tf  = self.v.get_transform()
        v   = self.speed()
        ld  = min(max(0.8 * v + 4.0, 6.0), 15.0)
        tx, ty = self.target_point(ld)
        dx, dy = tx - tf.location.x, ty - tf.location.y
        alpha = math.radians(wrap_deg(math.degrees(math.atan2(dy, dx)) - tf.rotation.yaw))
        delta = math.atan2(2.0 * self.wheelbase * math.sin(alpha), max(math.hypot(dx, dy), 1.0))
        steer = float(np.clip(math.degrees(delta) / self.max_steer_deg, -1.0, 1.0))

        err = v_ref - v
        self.i_term = float(np.clip(self.i_term + err * self.dt, 0.0, 16.0))
        u = 0.3 * err + 0.05 * self.i_term
        if v_ref < 0.1 and v < 0.5:
            throttle, brake = 0.0, 1.0
        elif u >= 0.0:
            throttle, brake = min(u, 0.9), 0.0
        elif err < -0.5:
            throttle, brake = 0.0, float(min(-u, 0.8))
        else:
            throttle, brake = 0.0, 0.0
        self.v.apply_control(carla.VehicleControl(
            throttle=float(throttle), steer=steer, brake=float(brake)))


def spawn_on_route(world, bp_lib, name, role, xy, z, yaw):
    bps = list(bp_lib.filter(name))
    if not bps:                                    # blueprint missing: any car
        bps = [b for b in bp_lib.filter("vehicle.*")
               if b.has_attribute("number_of_wheels")
               and b.get_attribute("number_of_wheels").as_int() == 4]
    bp = bps[0]
    if bp.has_attribute("role_name"):
        bp.set_attribute("role_name", role)
    tf = carla.Transform(carla.Location(x=float(xy[0]), y=float(xy[1]), z=float(z) + 0.5),
                         carla.Rotation(yaw=float(yaw)))
    return world.try_spawn_actor(bp, tf)


def route_yaw(pts, i):
    j = min(i + 1, len(pts) - 1)
    i = j - 1
    return math.degrees(math.atan2(pts[j, 1] - pts[i, 1], pts[j, 0] - pts[i, 0]))


def banner(t, text):
    print("=" * 66)
    print(f"[DEMO] t={t:7.2f} s  {text}")
    print("=" * 66)


# ---- Main ----------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description="Scripted cut-in scene for the radar demo")
    ap.add_argument("--town", default=DEMO_TOWN,
                    help="map to use (default %(default)s; 'current' keeps the loaded map)")
    ap.add_argument("--side", default=CUT_IN_SIDE, choices=["left", "right"])
    ap.add_argument("--seed", type=int, default=7, help="route choice (change for another road)")
    ap.add_argument("--load-map", metavar="MAP",
                    help="only load this map (e.g. Town10HD_Opt) and exit")
    args = ap.parse_args()

    client = carla.Client("localhost", 2000)
    client.set_timeout(60.0)

    if args.load_map:
        print(f"Loading {args.load_map} ...")
        client.load_world(args.load_map)
        print("Done.")
        return

    world = client.get_world()
    current = world.get_map().name.split("/")[-1]
    if args.town and args.town != "current" and not current.startswith(args.town):
        print(f"Loading {args.town} (currently {current}) ...")
        world = client.load_world(args.town)
        current = world.get_map().name.split("/")[-1]
    print(f"Map: {current}")

    dt = radar_io.SYNC_FIXED_DELTA
    original_settings = world.get_settings()
    actors = []
    radar = streamer = None

    sock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
    sock.setsockopt(socket.SOL_SOCKET, socket.SO_SNDBUF, 4 * 1024 * 1024)
    log_file = open("radar_log_cut_in.csv", "w", newline="")
    csv_writer = csv.writer(log_file)
    csv_writer.writerow(["timestamp", "num_points", "velocity", "azimuth", "altitude", "depth"])

    try:
        settings = world.get_settings()
        settings.synchronous_mode = True
        settings.fixed_delta_seconds = dt
        world.apply_settings(settings)
        print(f"[SYNC] synchronous mode, fixed step {dt} s")

        rng = np.random.default_rng(args.seed)
        route = find_route(world.get_map(), args.side, rng)
        if route is None:
            print("No long two-lane stretch found on this map. "
                  "Try --town Town04 or --side right.")
            return
        ego_pts, side_pts, road_z = route
        print(f"Route: {PATH_LENGTH_M:.0f} m, starting at "
              f"x={ego_pts[0, 0]:.0f} y={ego_pts[0, 1]:.0f}")

        bp_lib = world.get_blueprint_library()
        i_ego = int(round((START_BEHIND_M + 10.0) / PATH_SPACING_M))
        i_cut = i_ego - int(round(START_BEHIND_M / PATH_SPACING_M))
        ego = spawn_on_route(world, bp_lib, EGO_BLUEPRINT, "hero",
                             ego_pts[i_ego], road_z[i_ego], route_yaw(ego_pts, i_ego))
        cut = spawn_on_route(world, bp_lib, CUT_IN_BLUEPRINT, "cut_in",
                             side_pts[i_cut], road_z[i_cut], route_yaw(side_pts, i_cut))
        if ego is None or cut is None:
            print("Could not spawn the two vehicles on the route.")
            return
        actors += [ego, cut]
        print(f"Ego: {ego.type_id}   cut-in car: {cut.type_id} ({args.side} lane)")

        for _ in range(10):                  # let both cars settle on the road
            world.tick()

        ego_drv = RouteDriver(ego, ego_pts, side_pts, dt)
        cut_drv = RouteDriver(cut, ego_pts, side_pts, dt)
        ego_drv.lane_mix, cut_drv.lane_mix = 1.0, 0.0
        ego_drv.idx, cut_drv.idx = i_ego, i_cut
        for veh, pts, i in ((ego, ego_pts, i_ego), (cut, side_pts, i_cut)):
            yaw = math.radians(route_yaw(pts, i))
            veh.set_target_velocity(carla.Vector3D(
                EGO_SPEED * math.cos(yaw), EGO_SPEED * math.sin(yaw), 0.0))

        # Radar: same sensor settings and stream as carla_radar_sender.py
        radar_bp = bp_lib.find("sensor.other.radar")
        radar_bp.set_attribute("horizontal_fov",    str(radar_io.RADAR_HORIZONTAL_FOV))
        radar_bp.set_attribute("vertical_fov",      str(radar_io.RADAR_VERTICAL_FOV))
        radar_bp.set_attribute("range",             str(radar_io.RADAR_RANGE))
        radar_bp.set_attribute("points_per_second", str(radar_io.RADAR_POINTS_PER_SECOND))
        radar_bp.set_attribute("sensor_tick",       str(radar_io.RADAR_SENSOR_TICK))
        radar = world.spawn_actor(
            radar_bp,
            carla.Transform(carla.Location(x=2.2, z=1.0), carla.Rotation()),
            attach_to=ego)
        actors.append(radar)
        streamer = radar_io.RadarStreamer(sock, csv_writer)
        streamer.start()
        radar.listen(streamer.on_radar)
        print("Streaming IQ data cube to MATLAB on UDP port "
              f"{radar_io.MATLAB_PORT}. Press Ctrl+C to stop.")

        spectator = world.get_spectator()
        phase, t_phase = "free", 0.0
        cut_ref = OVERTAKE_SPEED
        t0 = world.get_snapshot().timestamp.elapsed_seconds
        banner(t0, f"FREE ROAD: no vehicle in front of the ego ({FREE_DRIVE_S:.0f} s)")
        next_t = time.perf_counter()

        while True:
            world.tick()
            next_t += dt
            delay = next_t - time.perf_counter()
            if delay > 0:
                time.sleep(delay)
            else:
                next_t = time.perf_counter()

            t = world.get_snapshot().timestamp.elapsed_seconds
            i_e, i_c = ego_drv.update_index(), cut_drv.update_index()
            gap = (i_c - i_e) * PATH_SPACING_M           # centre to centre, along the road
            v_ego, v_cut = ego_drv.speed(), cut_drv.speed()

            # ---- phase logic for the second car ----
            if phase == "free":
                v_cut_ref = v_ego + 0.3 * (-START_BEHIND_M - gap)   # hold position behind
                if t - t0 >= FREE_DRIVE_S:
                    phase, t_phase = "overtake", t
                    banner(t, "OVERTAKE: second car passes in the "
                              f"{args.side} lane")
            elif phase == "overtake":
                v_cut_ref = OVERTAKE_SPEED
                if gap >= CUT_IN_AT_M:
                    phase, t_phase = "cut_in", t
                    banner(t, f"CUT-IN: lane change {gap:.0f} m ahead of the ego")
            elif phase == "cut_in":
                v_cut_ref = OVERTAKE_SPEED
                cut_drv.lane_mix = smoothstep((t - t_phase) / LANE_CHANGE_S)
                if t - t_phase >= LANE_CHANGE_S:
                    cut_drv.lane_mix = 1.0
                    phase, t_phase = "brake", t
                    banner(t, f"BRAKING: car ahead slows to {CUT_IN_SPEED * 3.6:.0f} km/h, "
                              "ego brakes (closing target)")
            elif phase == "brake":
                cut_ref = max(CUT_IN_SPEED, cut_ref - CUT_IN_DECEL * dt)
                v_cut_ref = cut_ref
                if t - t_phase >= BRAKE_S:
                    phase, t_phase = "pull_away", t
                    banner(t, f"PULL AWAY: car ahead speeds up to {PULL_AWAY_SPEED * 3.6:.0f} km/h "
                              "(opening target)")
            elif phase == "pull_away":
                cut_ref = min(PULL_AWAY_SPEED, cut_ref + PULL_AWAY_ACCEL * dt)
                v_cut_ref = cut_ref
                if t - t_phase >= PULL_AWAY_S:
                    phase, t_phase = "stop", t
                    banner(t, "END: both cars stop")
            else:   # stop
                v_cut_ref = 0.0

            # ---- ego: cruise, distance keeping once the car is in its lane ----
            if phase == "stop":
                v_ego_ref = 0.0
            else:
                v_ego_ref = EGO_SPEED
                if cut_drv.lane_mix >= 0.5 and gap > 0.0:
                    v_ego_ref = float(np.clip(v_cut + 0.5 * (gap - FOLLOW_GAP_M), 0.0, EGO_SPEED))

            ego_drv.step(v_ego_ref)
            cut_drv.step(v_cut_ref)

            if SHOW_LABELS and phase != "free":
                label = {"overtake": "overtaking", "cut_in": "CUT-IN",
                         "brake": f"braking  gap {gap - 4.5:4.1f} m",
                         "pull_away": f"pulling away  gap {gap - 4.5:4.1f} m",
                         "stop": "stop"}[phase]
                world.debug.draw_string(
                    cut.get_transform().location + carla.Location(z=2.2), label,
                    draw_shadow=True, color=carla.Color(255, 200, 0), life_time=dt * 1.2)

            tf = ego.get_transform()
            yaw = math.radians(tf.rotation.yaw)
            spectator.set_transform(carla.Transform(
                tf.location + carla.Location(x=-CAM_BEHIND_M * math.cos(yaw),
                                             y=-CAM_BEHIND_M * math.sin(yaw),
                                             z=CAM_ABOVE_M),
                carla.Rotation(pitch=CAM_PITCH, yaw=tf.rotation.yaw)))

            if (phase == "stop" and (t - t_phase > 4.0 or v_ego < 0.2)) \
                    or min(ego_drv.remaining(), cut_drv.remaining()) < 10:
                if phase != "stop":
                    banner(t, "END: end of the route")
                break

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
        try:
            world.apply_settings(original_settings)
            print("[SYNC] world settings restored")
        except Exception:
            pass
        for actor in actors:
            try:
                if actor.is_alive:
                    actor.destroy()
            except Exception:
                pass
        log_file.close()
        sock.close()
        print(f"Done. CARLA is still on {current}. Before running carla_radar_sender.py, "
              "restart CARLA or run:  python cut_in_demo.py --load-map Town10HD_Opt")


if __name__ == "__main__":
    main()
