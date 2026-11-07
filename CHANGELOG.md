# Changelog

## v1.0.0: first public release

Changes since the initial version **`v0.1.0`**.

### MATLAB receiver (`matlab/carla_radar_full_live.m`)

- **Angle estimation: virtual-channel order.** The steering vectors used the element positions `[0 4 1 5 2 6 3 7]·d`, but the cube's channels are ordered `[0 1 2 … 7]·d`. Every angle estimate was therefore wrong. For example, a target at +10° appeared at +30° in the range-azimuth map and at −30° in CBF and MUSIC. The positions now follow the channel order.
- **Angle estimation: phase sign.** The steering-vector sign now matches the signal model. Previously CBF and MUSIC mirrored the angle even when the channel order was right.
- **Range-azimuth map.** The map now uses a real Doppler FFT and takes the slice at the dominant Doppler bin. Before, the Doppler-bin number was used as a chirp number, which selected one chirp.
- **CBF beampattern.** Each target now uses its own range-Doppler cell.
- **CARLA reference marker.** The dashed reference line now comes from the detection closest in both range and velocity.
- **CA-CFAR range limit.** Detection is limited to the 60 m sensor range. The noise-only region from 60 to 255 m produced most of the false alarms, and those took up the six MUSIC slots.
- **CA-CFAR ranking.** Detections are ranked by power before cluster suppression. Each cluster now keeps its peak cell, and the six strongest targets go to CBF and MUSIC.
- **Flush counter.** It is now measured; before, it was a constant 0. When two or more frames are queued, the receiver drops the queue and resynchronises on the next frame.
- **Packet handling.**
  - Chunks are locked to one `frame_id`; stale, foreign and duplicate chunks are discarded.
  - A new META packet abandons an unfinished frame instead of stalling for 3 s.
  - Queued packets are read in one batch, and a backlog is dropped with a single `flush`, so the receiver can keep up in real time.
  - Each frame log line shows the time spent receiving, processing and drawing.
- **Header size.** Corrected to 24 bytes (it was documented as 28).
- **Velocity axis.** It now matches the FFT bins exactly, so zero Doppler shows at 0 m/s instead of +0.77 m/s.
- **Waveform plot (Fig 1).** It now shows the baseband chirp. The 77 GHz carrier cannot be represented at 25.6 MHz sampling.
- **Window layout.** The ten figure windows are tiled around the CARLA window to fit the screen every time the script starts.
- **Display rate.** `DISPLAY_FPS` (default 2) caps how often MATLAB processes and redraws, so it no longer keeps the GPU busy all the time.
- **Fig 9 label.** It is now marked as a placeholder, because the cube carries no elevation phase.
- **Code clean-up.**
  - Removed the amplitude-weighting code that was computed but never used; it now lives in the Python sender.
  - Replaced `physconst`, so no toolbox is needed.

### Python sender (`python/carla_radar_sender.py`)

- **Detection filter.** The filter (4–60 m, 0.5 < |v| ≤ 40 m/s, forward-facing) is now applied *before* the cube is built. Previously every detection from 1–60 m went into the cube, including static returns.
- **Amplitude weighting.** The heuristic reflectivity weight (0.1 below 2 m/s) and the R⁻⁴ power factor are now applied. Previously every target had amplitude 1.
- **Non-blocking callback.** The CARLA sensor callback only copies the detections. A sender thread transmits the newest frame, dropping older ones if it falls behind. A separate logger thread writes `radar_log.csv`.
- **Synchronous mode.** CARLA runs in fixed 50 ms steps paced to real time, so the scene moves evenly and the radar sees a constant time step. The original settings are restored on exit.
- **Faster cube synthesis.** The cube is now built with one matrix product: about 15 ms instead of about 140 ms for 150 detections on the test machine.
- **Scene fixes.**
  - The spectator camera now sits 8 m behind and 4 m above the ego vehicle; it was 8 m above.
  - The pedestrian crossing factor is now set on `carla.World` before the walkers spawn. Before, the call went to the Traffic Manager and failed silently.
- **Tests.** Added `test_sender.py` for offline checks.
