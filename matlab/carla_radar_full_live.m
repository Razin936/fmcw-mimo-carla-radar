%% carla_radar_full_live.m
% =========================================================
%  Live CARLA Radar  +  TRUE FMCW MIMO Signal Processing
% =========================================================
%
%  TRUE MIMO layout  (TI AWR1843-style)
%    numTX = 2,  numRX = 4,  numVirtual = 8  uniform lambda/2 ULA
%
%  UDP protocol uses MAGIC numbers on every packet so arrival
%  order / stale packets never cause index errors.
%
%  EXTRAS MERGED FROM V5:
%   - Fig 9  : MUSIC Elevation Spectrum (per CFAR target)
%   - Fig 10 : FFT CBF Beampattern (normalised, GT azimuth marker)
%   - TDM-MIMO chirp scheduling (tx_now = mod(chirp-1, numTX))
%   - RCS weighting + path-loss amplitude per detection
%   - AWGN injected from SNR_dB (instead of fixed noise_sigma)
%   - conv2 fast CA-CFAR  (replaces slow double for-loop)
%   - Dominant Doppler bin for Range-Azimuth map
%   - 3D point cloud with altitude (Z axis)
%   - MDL rank estimator for MUSIC
%   - Elevation MUSIC (per CFAR target, RX sub-array)
%   - Frame budget timer: MUSIC skips if over budget
%   - FPS tracking + flush_count display
%   - Aggressive UDP flush (keep newest packet only)
%
%  10 FIGURE WINDOWS
%   Fig 1  - TX / RX / Beat waveform
%   Fig 2  - Range-Doppler Map
%   Fig 3  - CA-CFAR detections  (conv2 fast)
%   Fig 4  - FFT Range-Azimuth Map  (virtual aperture)
%   Fig 5  - MUSIC Azimuth Spectrum  (MDL rank)
%   Fig 6  - 3D Point Cloud  (X/Y/Z + velocity color)
%   Fig 7  - Raw CARLA detections  R vs V
%   Fig 8  - Info panel
%   Fig 9  - MUSIC Elevation Spectrum  (per CFAR target)
%   Fig 10 - FFT CBF Beampattern  (normalised, GT marker)
%
%  HOW TO RUN
%   1. Start CARLA simulator
%   2. python carla_radar_sender.py
%   3. Run THIS script in MATLAB

clear; clc; close all;

%% ---- FMCW MIMO radar parameters  (must match Python) ---------------
c          = physconst('LightSpeed');
fc         = 77e9;
BW         = 150e6;
T          = 10e-6;
slope      = BW / T;
lambda     = c / fc;

numTX      = 2;
numRX      = 4;
numVirt    = numTX * numRX;     % 8 virtual elements
numADC     = 256;
numChirps  = 128;

d_rx       = lambda / 2;
d_tx       = numRX * d_rx;
rx_pos     = (0:numRX-1).' * d_rx;
tx_pos     = [0; d_tx];
% Channel order must match the cube: (tx0,rx0..rx3), (tx1,rx0..rx3)
virt_pos   = reshape(rx_pos + tx_pos.', [], 1);   % [0 d 2d ... 7d]

dR   = c  / (2*BW);
Vmax = lambda / (4*T);
dV   = lambda / (2*numChirps*T);

Fs      = numADC / T;
dt      = 1 / Fs;
t_chirp = (0:numADC-1).' * dt;

R_ax = (0:numADC-1)   * dR;
V_ax = linspace(-Vmax, Vmax, numChirps);

% Azimuth steering matrix  [numVirt x 181]
ang_ax = (-90:1:90).';
% Sign follows the sender's model x_q = s*exp(+j*2*pi/lambda*p_q*sin(theta))
steer_az = exp(+1j * 2*pi/lambda * virt_pos * sind(ang_ax.'));

% Elevation axis + steering matrix  [numRX x 61]  -- FROM V5
el_ax  = (-30:1:30).';
d_el   = lambda / 2;
el_idx = (0:numRX-1).';
steer_el = exp(+1j * 2*pi/lambda * d_el * el_idx * sind(el_ax.'));

fprintf('FMCW MIMO Radar Parameters\n');
fprintf('  fc=%g GHz  BW=%g MHz  T=%g us  slope=%.2f MHz/us\n', ...
        fc/1e9, BW/1e6, T*1e6, slope/1e12);
fprintf('  TX=%d  RX=%d  Virtual=%d-element ULA\n', numTX, numRX, numVirt);
fprintf('  ADC=%d  Chirps=%d\n', numADC, numChirps);
fprintf('  dR=%.3f m  Vmax=%.2f m/s  dV=%.3f m/s\n\n', dR, Vmax, dV);

%% ---- Noise / SNR  -- FROM V5 ---------------------------------------
SNR_dB = 15;

%% ---- Filtering thresholds ------------------------------------------
MAX_ABS_VEL = 40;
MIN_DEPTH   = 4.0;
MAX_DEPTH   = 60.0;

%% ---- CA-CFAR parameters  (conv2 fast version)  -- FROM V5 ----------
numGuard   = 2;
numTrain   = 4;
P_fa       = 1e-4;
fullW      = 2*(numGuard+numTrain)+1;
guardW     = 2*numGuard+1;
Ncfar      = fullW^2 - guardW^2;
alpha_cfar = Ncfar * (P_fa^(-1/Ncfar) - 1);
MAX_TAR    = 6;

% conv2 filter masks  -- FROM V5
W_full  = ones(fullW, fullW);
W_guard = zeros(fullW, fullW);
W_guard(numTrain+1:numTrain+guardW, numTrain+1:numTrain+guardW) = 1;
W_train = W_full - W_guard;

%% ---- Frame budget timer  -- FROM V5 --------------------------------
FRAME_BUDGET_MS = 48;   % total ms per frame  (just under 50 ms = 20 Hz)
MAX_MUSIC_MS    = 25;   % MUSIC gets at most this slice

%% ---- UDP setup -----------------------------------------------------
UDP_PORT    = 5005;
CHUNK_SIZE  = 60000;
MAGIC_META  = uint32(hex2dec('CAFEBEEF'));
MAGIC_CHUNK = uint32(hex2dec('DEADF00D'));

HEADER_BYTES = 28;   % 4+4+4+8+4
CUBE_BYTES   = numVirt * numADC * numChirps * 8;  % complex64

try
    u = udpport("datagram","IPV4","LocalPort",UDP_PORT, ...
                "Timeout",2,"InputBufferSize",8*1024*1024);
catch
    try
        u = udpport("datagram","IPV4","LocalPort",UDP_PORT,"Timeout",2);
    catch ME
        error('Could not open UDP port %d: %s', UDP_PORT, ME.message);
    end
end
disp('MATLAB FMCW MIMO radar started.  Waiting for CARLA packets ...');

%% ---- Figure 1: TX / RX / Beat waveform -----------------------------
f1 = figure('Name','Fig 1 - Waveform','Color','k','Position',[0 700 620 300]);
ax1 = axes(f1); style_ax(ax1);
hl_tx  = plot(ax1, nan,nan,'c','LineWidth',1.2); hold(ax1,'on');
hl_rx  = plot(ax1, nan,nan,'m','LineWidth',1.0);
hl_mix = plot(ax1, nan,nan,'y','LineWidth',1.4);
xlabel(ax1,'Time (us)'); ylabel(ax1,'Amplitude');
title(ax1,'TX / RX / Beat Signal','Color','w');
legend(ax1,{'TX','RX (target 1)','Beat (virt ch 1)'}, ...
       'TextColor','w','Location','northeast');

%% ---- Figure 2: Range-Doppler Map -----------------------------------
f2 = figure('Name','Fig 2 - Range-Doppler','Color','k','Position',[620 700 620 300]);
ax2 = axes(f2); style_ax(ax2);
hi_rdm = imagesc(ax2, V_ax, R_ax, zeros(numADC,numChirps));
xlabel(ax2,'Velocity (m/s)'); ylabel(ax2,'Range (m)');
title(ax2,'Range-Doppler Map','Color','w');
cb2 = colorbar(ax2); cb2.Color='w'; cb2.Label.String='dB';
caxis(ax2,[-60 0]); colormap(ax2,jet(256));

%% ---- Figure 3: CA-CFAR --------------------------------------------
f3 = figure('Name','Fig 3 - CA-CFAR','Color','k','Position',[1240 700 620 300]);
ax3 = axes(f3); style_ax(ax3);
hi_cfar = imagesc(ax3, V_ax, R_ax, zeros(numADC,numChirps));
xlabel(ax3,'Velocity (m/s)'); ylabel(ax3,'Range (m)');
title(ax3,'CA-CFAR Detections','Color','w');
colormap(ax3,[0 0 0.3; 1 0.2 0.2]);

%% ---- Figure 4: FFT Range-Azimuth Map ------------------------------
f4 = figure('Name','Fig 4 - Range-Azimuth','Color','k','Position',[0 360 620 300]);
ax4 = axes(f4); style_ax(ax4);
hi_ra = imagesc(ax4, ang_ax, R_ax, zeros(numADC,numel(ang_ax)));
xlabel(ax4,'Azimuth (deg)'); ylabel(ax4,'Range (m)');
title(ax4,'FFT Range-Azimuth (virtual aperture beamforming)','Color','w');
cb4 = colorbar(ax4); cb4.Color='w'; cb4.Label.String='dB';
caxis(ax4,[-35 0]); colormap(ax4,jet(256));

%% ---- Figure 5: MUSIC Azimuth Spectrum -----------------------------
f5 = figure('Name','Fig 5 - MUSIC Azimuth','Color','k','Position',[620 360 620 300]);
ax5 = axes(f5); style_ax(ax5);
musicColors = lines(MAX_TAR);
hl_music = gobjects(MAX_TAR,1);
for k = 1:MAX_TAR
    hl_music(k) = plot(ax5, ang_ax, nan(size(ang_ax)), ...
        'Color',musicColors(k,:),'LineWidth',1.5);
    hold(ax5,'on');
end
xlabel(ax5,'Azimuth (deg)'); ylabel(ax5,'Power (dB)');
title(ax5,'MUSIC Azimuth Spectrum (MDL rank)','Color','w');
xlim(ax5,[-90 90]); ylim(ax5,[-40 5]);

%% ---- Figure 6: 3D Point Cloud  -- FROM V5 -------------------------
f6 = figure('Name','Fig 6 - 3D Point Cloud','Color','k','Position',[1240 360 620 300]);
ax6 = axes(f6); style_ax_3d(ax6);
hpc3 = scatter3(ax6, nan,nan,nan, 60, nan, 'filled');
xlabel(ax6,'X (m)'); ylabel(ax6,'Y (m)'); zlabel(ax6,'Z (m)');
title(ax6,'3D Point Cloud - X/Y/Z with velocity color','Color','w');
xlim(ax6,[-70 70]); ylim(ax6,[0 MAX_DEPTH]); zlim(ax6,[-15 15]);
cb6 = colorbar(ax6); cb6.Color='w'; cb6.Label.String='Velocity (m/s)';
colormap(ax6,jet(256));
view(ax6,-30,25);

%% ---- Figure 7: Raw CARLA detections --------------------------------
f7 = figure('Name','Fig 7 - Raw Detections','Color','k','Position',[0 40 620 280]);
ax7 = axes(f7); style_ax(ax7);
hrv = scatter(ax7, nan,nan, 30, nan, 'filled');
xlabel(ax7,'Radial Velocity (m/s)'); ylabel(ax7,'Depth (m)');
title(ax7,'Raw CARLA Detections - R vs V','Color','w');
xlim(ax7,[-MAX_ABS_VEL MAX_ABS_VEL]); ylim(ax7,[0 MAX_DEPTH]);
cb7 = colorbar(ax7); cb7.Color='w'; cb7.Label.String='Azimuth (deg)';
colormap(ax7,jet(256));

%% ---- Figure 8: Info panel -----------------------------------------
f8 = figure('Name','Fig 8 - Info','Color','k','Position',[620 40 620 280]);
ax8 = axes(f8); axis(ax8,'off');
ht_info = text(ax8,0.05,0.95,'Waiting for CARLA FMCW MIMO data...', ...
    'Color','w','FontSize',11,'Units','normalized', ...
    'VerticalAlignment','top','FontName','Courier New');

%% ---- Figure 9: MUSIC Elevation Spectrum  -- FROM V5 ---------------
f9 = figure('Name','Fig 9 - MUSIC Elevation','Color','k','Position',[1240 40 620 280]);
ax9 = axes(f9); style_ax(ax9);
hl_music_el = gobjects(MAX_TAR,1);
for k = 1:MAX_TAR
    hl_music_el(k) = plot(ax9, el_ax, nan(size(el_ax)), ...
        'Color',musicColors(k,:),'LineWidth',1.5);
    hold(ax9,'on');
end
xlabel(ax9,'Elevation (deg)'); ylabel(ax9,'Power (dB)');
title(ax9,'MUSIC Elevation Spectrum (per CFAR target)','Color','w');
xlim(ax9,[-30 30]); ylim(ax9,[-40 5]);

%% ---- Figure 10: FFT CBF Beampattern  -- FROM V5 -------------------
f10 = figure('Name','Fig 10 - FFT Beampattern','Color','k','Position',[1860 360 700 320]);
ax10 = axes(f10); style_ax(ax10);
hl_fft_spec = gobjects(MAX_TAR,1);
hl_fft_gt   = gobjects(MAX_TAR,1);
for k = 1:MAX_TAR
    hl_fft_spec(k) = plot(ax10, ang_ax, nan(size(ang_ax)), ...
        'Color',musicColors(k,:),'LineWidth',1.8);
    hold(ax10,'on');
    hl_fft_gt(k) = xline(ax10, nan, '--', ...
        'Color',musicColors(k,:),'LineWidth',1.0,'Alpha',0.7);
end
xlabel(ax10,'Azimuth (deg)');
ylabel(ax10,'Normalised power (dB)');
title(ax10,'FFT Beampattern  |  dashed = CARLA GT azimuth','Color','w');
xlim(ax10,[-90 90]); ylim(ax10,[-40 5]);
yline(ax10,-13,':','Color',[0.6 0.6 0.6],'LineWidth',0.8, ...
    'Label','-13 dB ULA sidelobe ref', ...
    'LabelHorizontalAlignment','left','FontSize',8);

drawnow;
allFigs    = [f1 f2 f3 f4 f5 f6 f7 f8 f9 f10];
frameCount = 0;
noise_floor_db = -60;
last_wall_t    = tic;   % for FPS tracking  -- FROM V5

%% ---- Pre-compute Hann windows (no Signal Processing Toolbox needed) -
n_range     = (0:numADC-1).';
n_chirp     = (0:numChirps-1).';
win_range   = 0.5 * (1 - cos(2*pi*n_range   / (numADC-1)));
win_doppler = 0.5 * (1 - cos(2*pi*n_chirp   / (numChirps-1))).';
win2d       = win_range * win_doppler;   % [numADC x numChirps]

%% ==================================================================
%% MAIN LOOP
%% ==================================================================
while all(isvalid(allFigs))

    frame_start = tic;   % budget timer  -- FROM V5

    %% Step 1: Wait for META packet (magic = 0xCAFEBEEF) -------------
    %
    %  META  (12 bytes): [uint32 MAGIC_META | uint32 num_chunks | uint32 total_bytes]
    %  CHUNK (>=16 bytes): [uint32 MAGIC_CHUNK | uint32 frame_id |
    %                       uint32 chunk_idx(0-based) | uint32 chunk_len | payload]

    got_meta    = false;
    num_chunks  = 0;
    total_bytes = 0;
    t_meta      = tic;

    while ~got_meta && toc(t_meta) < 3.0
        if u.NumDatagramsAvailable <= 0
            pause(0.005);
            continue;
        end
        dg  = read(u, 1, "uint8");
        dgd = uint8(dg.Data);
        if numel(dgd) < 12
            continue;
        end
        magic = typecast(dgd(1:4), 'uint32');
        if magic ~= MAGIC_META
            continue;
        end
        num_chunks  = double(typecast(dgd(5:8),  'uint32'));
        total_bytes = double(typecast(dgd(9:12), 'uint32'));
        if num_chunks >= 1 && total_bytes >= HEADER_BYTES
            got_meta = true;
        end
    end

    if ~got_meta
        continue;
    end

    %% Step 2: Collect all chunks for this frame ----------------------
    assembled = zeros(total_bytes, 1, 'uint8');
    received  = false(num_chunks, 1);
    t_wait    = tic;

    while ~all(received) && toc(t_wait) < 3.0
        if u.NumDatagramsAvailable <= 0
            pause(0.001);
            continue;
        end
        dg  = read(u, 1, "uint8");
        dgd = uint8(dg.Data);

        if numel(dgd) < 16
            continue;
        end

        magic = typecast(dgd(1:4), 'uint32');
        if magic ~= MAGIC_CHUNK
            continue;
        end

        c_idx = double(typecast(dgd(9:12),  'uint32')) + 1;  % 1-based
        c_len = double(typecast(dgd(13:16), 'uint32'));

        if c_idx < 1 || c_idx > num_chunks || c_len < 1
            continue;
        end
        if numel(dgd) < 16 + c_len
            continue;
        end

        payload    = dgd(17 : 16 + c_len);
        start_byte = (c_idx - 1) * CHUNK_SIZE + 1;
        end_byte   = start_byte + c_len - 1;

        if end_byte <= total_bytes
            assembled(start_byte : end_byte) = payload;
            received(c_idx) = true;
        end
    end

    if ~all(received)
        warning('Incomplete frame (%d/%d chunks). Skipping.', ...
                sum(received), num_chunks);
        continue;
    end

    %% Step 3: Parse fixed header (28 bytes) --------------------------
    if numel(assembled) < HEADER_BYTES
        warning('Frame too short for header. Skipping.');
        continue;
    end

    p = 1;
    pkt_numVirt  = double(typecast(assembled(p:p+3),  'uint32')); p = p+4;
    pkt_numADC   = double(typecast(assembled(p:p+3),  'uint32')); p = p+4;
    pkt_numChirp = double(typecast(assembled(p:p+3),  'uint32')); p = p+4;
    timestamp    =        typecast(assembled(p:p+7),  'double');  p = p+8;
    num_det      = double(typecast(assembled(p:p+3),  'uint32')); p = p+4;

    if pkt_numVirt ~= numVirt || pkt_numADC ~= numADC || pkt_numChirp ~= numChirps
        warning('Dimension mismatch in packet. Skipping.');
        continue;
    end

    %% Step 4: Parse raw CARLA detections -----------------------------
    raw_vals = [];
    if num_det > 0
        det_end = p + num_det*16 - 1;
        if det_end > numel(assembled)
            warning('Packet too short for detections. Skipping.');
            continue;
        end
        raw_float = typecast(assembled(p:det_end), 'single');
        raw_vals  = double(reshape(raw_float, 4, []).');
        % columns: [velocity, azimuth_rad, altitude_rad, depth_m]
        p = det_end + 1;
    end

    %% Step 5: Parse IQ Radar Data Cube [numVirt x numADC x numChirps]-
    cube_end = p + CUBE_BYTES - 1;
    if cube_end > numel(assembled)
        warning('Packet too short for IQ cube. Skipping.');
        continue;
    end

    cube_single = typecast(assembled(p:cube_end), 'single');
    cube_cplx   = complex(cube_single(1:2:end), cube_single(2:2:end));
    % Python C row-major [numVirt, numADC, numChirps] -> permute for MATLAB
    RDC = reshape(cube_cplx, [numChirps, numADC, numVirt]);
    RDC = permute(RDC, [2 1 3]);   % -> [numADC x numChirps x numVirt]

    %% Step 6: Filter raw CARLA detections ----------------------------
    vel_raw = []; dep_raw = []; az_deg = []; alt_deg = [];
    xp = []; yp = []; zp = [];

    if ~isempty(raw_vals)
        vel_r = raw_vals(:,1);
        az_r  = raw_vals(:,2);
        alt_r = raw_vals(:,3);
        dep_r = raw_vals(:,4);
        az_d  = rad2deg(az_r);
        alt_d = rad2deg(alt_r);
        x_r   = dep_r .* sin(az_r);
        y_r   = dep_r .* cos(az_r);
        z_r   = dep_r .* sin(alt_r);   % altitude -> Z  -- FROM V5

        ok = isfinite(x_r) & isfinite(y_r) & isfinite(vel_r) ...
           & y_r  >= 0 ...
           & dep_r >= MIN_DEPTH & dep_r <= MAX_DEPTH ...
           & abs(vel_r) <= MAX_ABS_VEL ...
           & abs(vel_r) > 0.5;

        vel_raw = vel_r(ok);
        dep_raw = dep_r(ok);
        az_deg  = az_d(ok);
        alt_deg = alt_d(ok);
        xp      = x_r(ok);
        yp      = y_r(ok);
        zp      = z_r(ok);
    end

    %% Update Fig 7 ---------------------------------------------------
    if ~isempty(dep_raw)
        set(hrv, 'XData', vel_raw, 'YData', dep_raw, 'CData', az_deg);
        caxis(ax7, [-90 90]);
    end

    %% Step 7: RCS weighting + path loss  -- FROM V5 -----------------
    flush_count = 0;   % not applicable in chunked protocol but kept for info panel
    if ~isempty(dep_raw)
        rcs_weight            = ones(size(dep_raw));
        rcs_weight(abs(vel_raw) < 2.0) = 0.1;
        path_loss  = 1 ./ (dep_raw .^ 4);
        path_loss  = path_loss / max(path_loss);
        amplitude  = rcs_weight .* sqrt(path_loss);
    end

    %% Step 8: Fig 1 - Waveform ---------------------------------------
    sig_tx = real(exp(1j*2*pi*(fc*t_chirp + slope/2*t_chirp.^2)));

    if ~isempty(dep_raw)
        delay1 = 2 * dep_raw(1) / c;
        td     = t_chirp - delay1;
        sig_rx = real(exp(1j*2*pi*(fc*td + slope/2*td.^2)));
    else
        sig_rx = zeros(numADC, 1);
    end

    sig_mix = real(RDC(:, 1, 1));
    mx = max(abs(sig_mix));
    if mx > 0
        sig_mix = sig_mix / mx;
    end

    t_us = t_chirp * 1e6;
    set(hl_tx,  'XData', t_us, 'YData', sig_tx);
    set(hl_rx,  'XData', t_us, 'YData', sig_rx);
    set(hl_mix, 'XData', t_us, 'YData', sig_mix);
    xlim(ax1, [0 t_us(end)]);

    %% Step 9: Apply AWGN from SNR_dB to RDC  -- FROM V5 -------------
    sig_power = mean(abs(RDC(:)).^2);
    if sig_power > 0
        noise_var = sig_power / (10^(SNR_dB/10));
        sigma_n   = sqrt(noise_var / 2);
        RDC = RDC + sigma_n * (randn(size(RDC)) + 1j*randn(size(RDC)));
    end

    %% Step 10: 2D FFT -> Range-Doppler Map (virtual ch 1) -----------
    ch1    = RDC(:,:,1) .* win2d;
    RDM    = fftshift(fft2(ch1, numADC, numChirps), 2);
    RDM_dB = 20*log10(abs(RDM) ./ (max(abs(RDM(:))) + eps));
    set(hi_rdm, 'CData', RDM_dB);
    title(ax2, sprintf('Range-Doppler Map  |  t = %.2f s', timestamp), 'Color','w');

    %% Step 11: CA-CFAR via conv2  -- FROM V5 ------------------------
    RDM_pow   = abs(RDM).^2;
    noise_sum = conv2(RDM_pow, W_train, 'same');
    cfar_mask = RDM_pow > (alpha_cfar / Ncfar) * noise_sum;

    % Zero border cells
    border = numTrain + numGuard;
    cfar_mask(1:border, :)           = 0;
    cfar_mask(end-border+1:end, :)   = 0;
    cfar_mask(:, 1:border)           = 0;
    cfar_mask(:, end-border+1:end)   = 0;
    cfar_mask(R_ax > MAX_DEPTH + 2*dR, :) = 0;   % beyond the sensor range: noise only

    set(hi_cfar, 'CData', double(cfar_mask));
    % Rank cells by power so cluster suppression keeps each cluster's peak
    % and the MAX_TAR strongest targets are passed to CBF / MUSIC.
    det_idx = find(cfar_mask);
    [~, order] = sort(RDM_pow(det_idx), 'descend');
    [cfar_r_all, cfar_d_all] = ind2sub(size(cfar_mask), det_idx(order));
    [cfar_r, cfar_d] = remove_close_detections(cfar_r_all, cfar_d_all, 5);
    K = numel(cfar_r);
    title(ax3, sprintf('CA-CFAR (conv2)  |  %d detections', K), 'Color','w');

    %% Step 12: Range FFT across all virtual channels -----------------
    rangeFFT = zeros(numADC, numChirps, numVirt);
    for vv = 1:numVirt
        chan_win          = RDC(:,:,vv) .* win2d;
        rangeFFT(:,:,vv) = fft(chan_win, numADC, 1);
    end
    dopFFT = fftshift(fft(rangeFFT, numChirps, 2), 2);   % [range x Doppler x virt]

    %% Step 13: Fig 4 - Range-Azimuth using dominant Doppler bin  -- FROM V5
    [~, dom_dop_bin] = max(max(abs(RDM), [], 1));
    virt_snap = squeeze(dopFFT(:, dom_dop_bin, :));    % [numADC x numVirt]
    ra_map = abs(virt_snap * conj(steer_az));          % [numADC x 181]
    ra_dB  = 20*log10(ra_map ./ (max(ra_map(:)) + eps));
    set(hi_ra, 'CData', ra_dB);
    title(ax4, sprintf('FFT Range-Azimuth  |  %d virtual elements', numVirt), ...
          'Color','w');

    %% Step 14: Fig 6 - 3D Point Cloud  -- FROM V5 -------------------
    if ~isempty(xp)
        set(hpc3, 'XData', xp, 'YData', yp, 'ZData', zp, 'CData', vel_raw);
        caxis(ax6, [-MAX_ABS_VEL MAX_ABS_VEL]);
    end
    title(ax6, sprintf('3D Point Cloud  |  %d pts', numel(xp)), 'Color','w');

    %% Step 15: Fig 10 - CBF Beampattern  -- FROM V5 -----------------
    for k = 1:MAX_TAR
        set(hl_fft_spec(k), 'YData', nan(size(ang_ax)));
        set(hl_fft_gt(k),   'Value', nan);
    end

    numTarFFT = min(K, MAX_TAR);
    if numTarFFT > 0
        for k = 1:numTarFFT
            ri       = cfar_r(k);
            x_snap   = squeeze(dopFFT(ri, cfar_d(k), :));   % own range-Doppler cell
            bf_power = abs(steer_az' * x_snap).^2;
            bf_dB    = 10*log10(bf_power ./ (max(bf_power) + eps));
            set(hl_fft_spec(k), 'YData', bf_dB);

            % Ground-truth azimuth from nearest CARLA detection
            if ~isempty(dep_raw)
                cost = abs(dep_raw - R_ax(ri))/dR + abs(vel_raw - V_ax(cfar_d(k)))/dV;
                [~, nearest_det] = min(cost);
                set(hl_fft_gt(k), 'Value', az_deg(nearest_det));
            end
        end
        ylim(ax10, [-40 5]);
    end
    title(ax10, sprintf('FFT Beampattern  |  %d target(s)  |  peak = 0 dB', ...
          numTarFFT), 'Color','w');

    %% Step 16: MUSIC - always runs when targets detected ---------------
    % NOTE: frame budget gate removed because chunked UDP reassembly
    % already consumes >> FRAME_BUDGET_MS, so budget check always blocked
    % MUSIC. MUSIC now runs unconditionally whenever K > 0.

    % Reset all MUSIC lines
    for k = 1:MAX_TAR
        set(hl_music(k),    'YData', nan(size(ang_ax)));
        set(hl_music_el(k), 'YData', nan(size(el_ax)));
    end

    numTar = min(K, MAX_TAR);

    if numTar > 0
        for k = 1:numTar
            ri = cfar_r(k);

            %% Azimuth MUSIC (full virtual aperture)  ----------------
            A_az   = squeeze(rangeFFT(ri, :, :)).';   % [numVirt x numChirps]
            Rxx_az = (A_az * A_az') / numChirps;

            % Forward-backward averaging
            J      = fliplr(eye(numVirt));
            Rxx_fb = 0.5 * (Rxx_az + J * conj(Rxx_az) * J);

            [Q_az, D_az] = eig(Rxx_fb);
            ev_az  = sort(real(diag(D_az)), 'descend');
            d_mdl  = estimate_rank_mdl(ev_az, numChirps);  % MDL  -- FROM V5
            d_mdl  = min(d_mdl, numVirt - 2);   % ensure >=2 noise eigenvectors
            Q_az   = fliplr(Q_az);
            Qn_az  = Q_az(:, d_mdl+1 : end);
            QnQnh  = Qn_az * Qn_az';

            music_az = zeros(numel(ang_ax), 1);
            for ai = 1:numel(ang_ax)
                av = steer_az(:, ai);
                denom = real(av' * QnQnh * av);
                music_az(ai) = real((av' * av) / (denom + eps));
            end
            music_az_dB = 10*log10(abs(music_az) + eps);
            peak_az = max(music_az_dB);
            if isfinite(peak_az)
                music_az_dB = music_az_dB - peak_az;
                set(hl_music(k), 'YData', music_az_dB);
            end

            %% Elevation MUSIC (RX sub-array, TX0 chirps)  -- FROM V5 -
            tx0_chirps = 1:numTX:numChirps;
            rx_snap    = zeros(numRX, numel(tx0_chirps));
            for rx = 1:numRX
                rx_snap(rx,:) = rangeFFT(ri, tx0_chirps, rx);
            end

            Rxx_el   = (rx_snap * rx_snap') / size(rx_snap, 2);
            [Q_el, D_el] = eig(Rxx_el);
            ev_el    = sort(real(diag(D_el)), 'descend');
            d_mdl_el = estimate_rank_mdl(ev_el, size(rx_snap, 2));
            d_mdl_el = min(d_mdl_el, numRX - 2);   % ensure >=2 noise eigenvectors
            Q_el     = fliplr(Q_el);
            Qn_el    = Q_el(:, d_mdl_el+1 : end);
            QnQnh_el = Qn_el * Qn_el';

            music_el = zeros(numel(el_ax), 1);
            for ei = 1:numel(el_ax)
                ev    = steer_el(:, ei);
                denom = real(ev' * QnQnh_el * ev);
                music_el(ei) = real((ev' * ev) / (denom + eps));
            end
            music_el_dB = 10*log10(abs(music_el) + eps);
            peak_el = max(music_el_dB);
            if isfinite(peak_el)
                music_el_dB = music_el_dB - peak_el;
                set(hl_music_el(k), 'YData', music_el_dB);
            end
        end

        ylim(ax5, [-40 5]);
        ylim(ax9, [-40 5]);
    end

    title(ax5,  sprintf('MUSIC Azimuth  |  %d targets (MDL)', numTar),  'Color','w');
    title(ax9,  sprintf('MUSIC Elevation  |  %d targets', numTar),      'Color','w');

    %% Step 17: Fig 8 - Info panel ------------------------------------
    noise_floor_db = 10*log10(mean(abs(RDC(:)).^2) + eps);
    fps = 1 / toc(last_wall_t);   % -- FROM V5
    last_wall_t = tic;

    update_info(ht_info, timestamp, numel(vel_raw), K, numTar, ...
                frameCount, dR, dV, Vmax, numVirt, noise_floor_db, SNR_dB, flush_count);

    frameCount = frameCount + 1;
    fprintf('Frame %3d | t=%7.2f | pts=%3d | CFAR=%2d | MUSIC=%d | %.1f fps\n', ...
            frameCount, timestamp, numel(vel_raw), K, numTar, fps);

    drawnow limitrate;
end

disp('A figure was closed.  Radar stopped.');
clear u;

%% ==================================================================
%% LOCAL HELPER FUNCTIONS
%% ==================================================================

function style_ax(ax)
    set(ax,'Color','k','XColor','w','YColor','w', ...
           'GridColor',[.35 .35 .35],'GridAlpha',.6);
    grid(ax,'on'); hold(ax,'on');
end

function style_ax_3d(ax)
    set(ax,'Color','k','XColor','w','YColor','w','ZColor','w', ...
           'GridColor',[.35 .35 .35],'GridAlpha',.6);
    grid(ax,'on'); hold(ax,'on');
end

function update_info(ht, ts, nPts, nCfar, nMus, frame, dR, dV, Vmax, nVirt, nf, snr, flushed)
    str = sprintf([ ...
        ' Frame          : %d\n', ...
        ' CARLA time     : %.2f s\n', ...
        '\n', ...
        ' Filtered pts   : %d\n', ...
        ' CFAR det       : %d\n', ...
        ' MUSIC targets  : %d\n', ...
        ' Flushed pkts   : %d\n', ...
        '\n', ...
        ' --- TRUE FMCW MIMO ---\n', ...
        ' fc   = 77 GHz\n', ...
        ' BW   = 150 MHz\n', ...
        ' dR   = %.3f m\n', ...
        ' dV   = %.3f m/s\n', ...
        ' Vmax = %.2f m/s\n', ...
        ' Virt = %d elements\n', ...
        ' SNR  = %d dB\n', ...
        ' Noise floor ~ %.1f dB\n'], ...
        frame, ts, nPts, nCfar, nMus, flushed, ...
        dR, dV, Vmax, nVirt, snr, nf);
    set(ht, 'String', str);
end

function [cr, cd] = remove_close_detections(r_in, d_in, thr)
    keep = true(numel(r_in), 1);
    for i = 2:numel(r_in)
        for j = 1:i-1
            if keep(j) ...
               && abs(r_in(i)-r_in(j)) <= thr ...
               && abs(d_in(i)-d_in(j)) <= thr
                keep(i) = false;
                break;
            end
        end
    end
    cr = r_in(keep);
    cd = d_in(keep);
end

function k_hat = estimate_rank_mdl(eigenvalues, N)
%  MDL (Minimum Description Length) rank estimator  -- FROM V5
%  Wax & Kailath 1985.  Returns estimated number of signal sources.
    L   = numel(eigenvalues);
    mdl = zeros(L-1, 1);
    for k = 0:L-2
        noise_eigs = eigenvalues(k+1:end);
        noise_eigs = max(noise_eigs, eps);
        geo_mean   = exp(mean(log(noise_eigs)));
        arith_mean = mean(noise_eigs);
        ratio      = geo_mean / (arith_mean + eps);
        mdl(k+1)   = -N*(L-k)*log(ratio) + 0.5*k*(2*L-k)*log(N);
    end
    [~, idx] = min(mdl);
    k_hat = max(1, min(idx, L-1));
end