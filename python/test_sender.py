"""
Offline checks for carla_radar_sender.py (no CARLA or MATLAB needed).

    python test_sender.py

1. The vectorised cube equals a straightforward per-detection loop.
2. Filtering and amplitude weighting behave as documented.
3. A frame sent over UDP loopback reassembles to the same header,
   detections and cube (parsed the same way the MATLAB receiver does).
4. Timing of one cube synthesis.
"""
import socket
import struct
import sys
import time
import types

import numpy as np

sys.modules.setdefault("carla", types.ModuleType("carla"))   # stub: not needed offline
import carla_radar_sender as s  # noqa: E402


def reference_cube(dets):
    """Plain loop over detections and virtual elements (noise-free)."""
    dets = np.asarray(dets, dtype=np.float64)
    dets = dets[s.filter_detections(dets)]
    out = np.zeros((s.NUM_VIRTUAL, s.NUM_ADC, s.NUM_CHIRPS), dtype=np.complex128)
    if not len(dets):
        return out
    amp = s.amplitude_weights(dets[:, 0], dets[:, 3])
    for (v, az, _alt, R), a in zip(dets, amp):
        fd = 2 * v * s.FC / s.C
        fb = 2 * R * s.SLOPE / s.C + fd
        beat = np.exp(1j * 2 * np.pi * fb * s.T_AXIS)
        dop = 2 * np.pi * fd * s.CHIRP_IDX * s.T_CHIRP
        for q in range(s.NUM_VIRTUAL):
            ph = 2 * np.pi / s.LAMBDA * s.VIRT_POS[q] * np.sin(az)
            out[q] += a * beat[:, None] * np.exp(1j * (dop[None, :] + ph))
    return out


class ZeroNoise:
    def standard_normal(self, shape):
        return np.zeros(shape)


def random_dets(n, seed=0):
    r = np.random.default_rng(seed)
    return np.column_stack([
        r.uniform(-45, 45, n),                 # velocity (some outside the gate)
        r.uniform(-0.3, 0.3, n),               # azimuth (rad)
        r.uniform(-0.09, 0.09, n),             # altitude (rad)
        r.uniform(0.5, 70, n),                 # depth (some outside the gate)
    ]).astype(np.float32)


def test_cube_matches_reference():
    d = random_dets(40)
    got = s.build_rdc(d, rng=ZeroNoise())
    ref = reference_cube(d)
    err = np.abs(got - ref).max() / np.abs(ref).max()
    assert err < 1e-5, err
    print(f"ok  vectorised cube == loop reference   (max rel. error {err:.1e})")


def test_filter_and_weights():
    d = np.array([[-13, 0.1, 0, 20],     # kept
                  [0.2, 0.0, 0, 20],     # static -> dropped
                  [-13, 0.0, 0, 3.0],    # too close -> dropped
                  [-13, 0.0, 0, 65],     # too far -> dropped
                  [-45, 0.0, 0, 20],     # too fast -> dropped
                  [1.2, 0.0, 0, 10],     # slow walker -> kept, weight 0.1
                  [np.nan, 0, 0, 20]],   # non-finite -> dropped
                 dtype=np.float64)
    m = s.filter_detections(d)
    assert m.tolist() == [True, False, False, False, False, True, False], m
    a = s.amplitude_weights(d[m, 0], d[m, 3])
    # nearest (10 m, slow) -> 0.1 * 1 ; 20 m fast -> 1 * (10/20)^2
    assert np.allclose(a, [0.25, 0.1]), a
    assert not np.any(s.build_rdc(d[~m & np.isfinite(d).all(1)], rng=ZeroNoise())), \
        "filtered-out detections must not enter the cube"
    print("ok  filter gates and amplitude weights")


def test_udp_round_trip():
    rx = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
    rx.setsockopt(socket.SOL_SOCKET, socket.SO_RCVBUF, 8 * 1024 * 1024)
    rx.bind(("127.0.0.1", 0))
    rx.settimeout(2.0)
    s.MATLAB_PORT = rx.getsockname()[1]
    tx = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)

    d = random_dets(25, seed=3)
    np.random.seed(7)
    s.send_frame(tx, d, 545.31, frame_id=42)

    meta = rx.recv(65535)
    magic, n_chunks, total = struct.unpack("<III", meta)
    assert magic == s.MAGIC_META
    buf = bytearray(total)
    for _ in range(n_chunks):
        pkt = rx.recv(65535)
        mg, fid, idx, ln = struct.unpack("<IIII", pkt[:16])
        assert mg == s.MAGIC_CHUNK and fid == 42
        buf[idx * s.CHUNK_SIZE: idx * s.CHUNK_SIZE + ln] = pkt[16:16 + ln]

    nv, na, nc, ts, nd = struct.unpack(s.HEADER_FMT, buf[:24])
    assert (nv, na, nc, nd) == (8, 256, 128, 25) and abs(ts - 545.31) < 1e-12
    p = 24
    det_back = np.frombuffer(bytes(buf[p:p + nd * 16]), "<f4").reshape(-1, 4)
    assert np.array_equal(det_back, d)
    p += nd * 16
    flat = np.frombuffer(bytes(buf[p:]), "<f4")
    # Same steps as the MATLAB receiver: interleaved re/im, reshape
    # column-major to [chirps, adc, virt], then permute to [adc, chirps, virt].
    cplx = flat[0::2] + 1j * flat[1::2]
    cube_m = cplx.reshape((nc, na, nv), order="F").transpose(1, 0, 2)
    np.random.seed(7)
    expected = s.build_rdc(d)                     # [virt, adc, chirps]
    assert np.allclose(cube_m, expected.transpose(1, 2, 0))
    print(f"ok  UDP round trip ({n_chunks} chunks, {total} bytes) and MATLAB-order reshape")


def test_timing():
    d = random_dets(150, seed=5)
    d[:, 0] = np.random.default_rng(1).uniform(-15, -1, 150)   # all inside the gates
    d[:, 3] = np.random.default_rng(2).uniform(5, 59, 150)
    t = time.perf_counter()
    for _ in range(5):
        s.build_rdc(d)
    dt = (time.perf_counter() - t) / 5
    print(f"ok  build_rdc with 150 detections: {dt * 1000:.0f} ms per frame")


class FakeDet:
    def __init__(self, row):
        self.velocity, self.azimuth, self.altitude, self.depth = map(float, row)


class FakeMeasurement(list):
    def __init__(self, rows, ts):
        super().__init__(FakeDet(r) for r in rows)
        self.timestamp = ts


def test_streamer_threads():
    import csv
    import io
    rx = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
    rx.setsockopt(socket.SOL_SOCKET, socket.SO_RCVBUF, 32 * 1024 * 1024)
    rx.bind(("127.0.0.1", 0))
    s.MATLAB_PORT = rx.getsockname()[1]
    tx = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
    out = io.StringIO()
    st = s.RadarStreamer(tx, csv.writer(out))
    st.start()
    t = time.perf_counter()
    for k in range(10):                      # 10 frames in a burst, faster than real time
        st.on_radar(FakeMeasurement(random_dets(12, seed=k), 100.0 + 0.1 * k))
    callback_ms = (time.perf_counter() - t) * 1000 / 10
    time.sleep(1.5)
    st.stop()
    rows = out.getvalue().strip().splitlines()
    assert len(rows) == 10 * 12, len(rows)                   # every frame logged
    assert st.sent + st.dropped == 10 and st.sent >= 1, (st.sent, st.dropped)
    print(f"ok  streamer: callback {callback_ms:.2f} ms, sent={st.sent}, dropped(newest-only)={st.dropped}, csv rows={len(rows)}")


if __name__ == "__main__":
    test_cube_matches_reference()
    test_filter_and_weights()
    test_udp_round_trip()
    test_timing()
    test_streamer_threads()
    print("all sender checks passed")
