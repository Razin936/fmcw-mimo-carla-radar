# Demo videos

Screen recordings of the CARLA scene (left) and the ten MATLAB windows (right).

| Folder | Version | What it shows |
|---|---|---|
| [`v3/`](v3) | v3: TDM and elevation | 3 TX with TDM, 12 virtual channels. Fig 9 shows the elevation estimate next to CARLA's altitude (within about 0.3°). |

| Video | What it shows |
|---|---|
| `v3/demo_v3_city_drive_1.mp4` | City drive: range-Doppler, CA-CFAR, MUSIC azimuth and elevation |
| `v3/demo_v3_city_drive_2.mp4` | Second city drive: elevation within about 0.3° of CARLA's altitude |
| `v3/demo_v3_tuned.mp4` | Final receiver settings: 97–98 % of frames processed at 10 Hz |

Run the matching sender and receiver to reproduce a video:
`python/carla_radar_sender_v3.py` with `matlab/carla_radar_full_live_v3.m`.