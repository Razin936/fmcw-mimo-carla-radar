# MATLAB receiver

| File | What it is |
|---|---|
| `carla_radar_full_live.m` | **v1**: base version, processes about 2 frames per s |
| `carla_radar_full_live_v2.m` | **v2**: processes every frame (10 per s), redraws every 1 s |
| `carla_radar_full_live_v3.m` | **v3**: v2 + TDM compensation and elevation MUSIC (Fig 9) |

Run a receiver with the Python sender of the same version, for example
`carla_radar_full_live_v3.m` with `python/carla_radar_sender_v3.py`.
Start MATLAB first and wait until it says it is waiting for packets.
