# Demo videos

Screen recordings of the CARLA scene (left) and the ten MATLAB windows (right).

| Folder | Version | What it shows |
|---|---|---|
| [`v1/`](v1) | v1: base pipeline | Range-Doppler, CA-CFAR, MUSIC and CBF azimuth. Velocities still show the ego speed twice; about 2 frames per s are processed. |
| [`v2/`](v2) | v2: correct velocities | Roadside objects at the real ego speed. Every 10 Hz frame is processed (100 %); the windows are redrawn every 1 s. |
| [`v3/`](v3) | v3: TDM and elevation | 3 TX with TDM, 12 virtual channels. Fig 9 shows the elevation estimate next to CARLA's altitude (within about 0.3°). |

| Video | What it shows |
|---|---|
| `v1/demo_v1_city_drive.mp4` | City drive with the original pipeline |
| `v2/demo_v2_city_drive.mp4` | City drive with correct velocities, every frame processed |
| `v3/demo_v3_city_drive_1.mp4` | City drive: range-Doppler, CA-CFAR, MUSIC azimuth and elevation |
| `v3/demo_v3_city_drive_2.mp4` | Second city drive: elevation within about 0.3° of CARLA's altitude |
| `v3/demo_v3_tuned.mp4` | Final receiver settings: 97–98 % of frames processed at 10 Hz |

Run the sender and receiver of the same version to reproduce a video, for example
`python/carla_radar_sender_v3.py` with `matlab/carla_radar_full_live_v3.m`.