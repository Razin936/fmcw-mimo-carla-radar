# Python sender

| File | What it is |
|---|---|
| `carla_radar_sender.py` | **v1**: base version (velocities show the ego speed twice) |
| `carla_radar_sender_v2.py` | **v2**: correct velocities |
| `carla_radar_sender_v3.py` | **v3**: v2 + TDM and elevation (3 TX, 12 channels) |
| `cut_in_demo.py` | Scripted cut-in scene for demos |
| `test_sender.py` | Offline tests, no CARLA needed |

Run a sender with the MATLAB receiver of the same version, for example
`carla_radar_sender_v3.py` with `matlab/carla_radar_full_live_v3.m`.
