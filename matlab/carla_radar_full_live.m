%% carla_radar_full_live.m
% =========================================================================
%  Live CARLA radar  ->  FMCW MIMO signal processing in MATLAB
% =========================================================================
%
%  Receives, per CARLA radar tick, the raw CARLA detections and a
%  synthetic IQ radar data cube [256 ADC x 128 chirps x 8 virtual channels]
%  from python/carla_radar_sender.py and runs the processing chain:
%
%    AWGN injection -> 2-D Hann window -> range FFT -> Doppler FFT
%    -> range-Doppler map -> 2-D CA-CFAR (conv2, limited to the 60 m
%       sensor range) -> power-ranked cluster suppression
%    -> FFT range-azimuth map (dominant Doppler bin)
%    -> per-target CBF beampattern  +  MDL-gated MUSIC (azimuth)
%
%  Array: 2 TX x 4 RX -> 8-element uniform lambda/2 virtual ULA.
%  The sender synthesises all 8 virtual channels on every chirp (ideal
%  MIMO, no TDM time offset), so Vmax = lambda/(4*Tc) = 97.4 m/s.
%
%  10 FIGURE WINDOWS
%   Fig 1  - Baseband TX / RX / beat waveform
%   Fig 2  - Range-Doppler map
%   Fig 3  - CA-CFAR detections
%   Fig 4  - FFT range-azimuth map (virtual aperture, dominant Doppler bin)
%   Fig 5  - MUSIC azimuth spectrum (MDL rank, forward-backward averaged)
%   Fig 6  - 3-D point cloud of the filtered CARLA detections
%   Fig 7  - Raw CARLA detections, depth vs radial velocity
%   Fig 8  - Info panel
%   Fig 9  - RX sub-array MUSIC ("elevation" panel, see note below)
%   Fig 10 - CBF beampattern per CFAR target, dashed = CARLA azimuth
%
%  NOTE on Fig 9: the synthetic cube carries azimuth phase only, and the
%  4 RX elements lie on the same horizontal line as the azimuth array.
%  Fig 9 therefore does NOT measure elevation; it is kept as a
%  placeholder for a future 2-D array. See README "Known limitations".
%
%  HOW TO RUN
%   1. Start the CARLA simulator
%   2. python python/carla_radar_sender.py
%   3. Run this script (MATLAB R2020b or newer, no toolboxes required)

clear; clc; close all;

%% ---- FMCW MIMO radar parameters  (must match the Python sender) -------
c          = 299792458;              % speed of light (m/s); physconst needs a toolbox
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

dR   = c / (2*BW);                      % range bin spacing
dV   = lambda / (2*numChirps*T);        % velocity bin spacing
Vmax = lambda / (4*T);                  % unambiguous velocity (no TDM)

Fs      = numADC / T;
t_chirp = (0:numADC-1).' / Fs;

R_ax = (0:numADC-1) * dR;
V_ax = (-numChirps/2 : numChirps/2-1) * dV;    % matches fftshift bin order

% Steering vectors follow the sender's signal model
%   x_q = s * exp(+j*2*pi/lambda * p_q * sin(theta))
ang_ax   = (-90:1:90).';
steer_az = exp(+1j * 2*pi/lambda * virt_pos * sind(ang_ax.'));   % [numVirt x 181]

el_ax    = (-30:1:30).';
el_idx   = (0:numRX-1).';
steer_el = exp(+1j * 2*pi/lambda * (lambda/2) * el_idx * sind(el_ax.'));  % [numRX x 61]

fprintf('FMCW MIMO Radar Parameters\n');
fprintf('  fc=%g GHz  BW=%g MHz  T=%g us  slope=%.2f MHz/us\n', ...
        fc/1e9, BW/1e6, T*1e6, slope/1e12);
fprintf('  TX=%d  RX=%d  Virtual=%d-element ULA\n', numTX, numRX, numVirt);
fprintf('  ADC=%d  Chirps=%d\n', numADC, numChirps);
fprintf('  dR=%.3f m  Vmax=%.2f m/s  dV=%.3f m/s\n\n', dR, Vmax, dV);

%% ---- Noise / SNR --------------------------------------------------------
SNR_dB = 15;

%% ---- Display filter (same gates as the Python cube filter) -------------
MAX_ABS_VEL = 40;
MIN_ABS_VEL = 0.5;
MIN_DEPTH   = 4.0;
MAX_DEPTH   = 60.0;

%% ---- CA-CFAR parameters  (conv2 implementation) ------------------------
numGuard   = 2;
numTrain   = 4;
P_fa       = 1e-4;
fullW      = 2*(numGuard+numTrain)+1;            % 13
guardW     = 2*numGuard+1;                       % 5
Ncfar      = fullW^2 - guardW^2;                 % 144 training cells
alpha_cfar = Ncfar * (P_fa^(-1/Ncfar) - 1);
MAX_TAR    = 6;

W_train = ones(fullW, fullW);
W_train(numTrain+1:numTrain+guardW, numTrain+1:numTrain+guardW) = 0;

J_virt = fliplr(eye(numVirt));                   % forward-backward exchange matrix

%% ---- Display rate ---------------------------------------------------------
% Redrawing ten windows is the slowest part on a laptop and shares the GPU
% with CARLA. MATLAB processes and shows at most DISPLAY_FPS frames per
% second, stays idle in between and then continues with the newest frame.
% Raise it on a fast PC (10 = every CARLA frame).
DISPLAY_FPS = 2;

%% ---- UDP setup -----------------------------------------------------------
UDP_PORT     = 5005;
CHUNK_SIZE   = 60000;
MAGIC_META   = uint32(hex2dec('CAFEBEEF'));
MAGIC_CHUNK  = uint32(hex2dec('DEADF00D'));
HEADER_BYTES = 24;   % uint32 x3 + double + uint32  (struct "<IIIdI")
CUBE_BYTES   = numVirt * numADC * numChirps * 8;   % complex64

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

% clim (R2022a+) with caxis fallback
if exist('clim','file') == 2 || exist('clim','builtin') == 5
    set_clim = @(ax,lim) clim(ax,lim);
else
    set_clim = @(ax,lim) caxis(ax,lim); %#ok<CAXIS>
end

%% ---- Figure 1: baseband TX / RX / beat waveform ------------------------
f1 = figure('Name','Fig 1 - Waveform','Color','k','Position',[0 700 620 300]);
ax1 = axes(f1); style_ax(ax1);
hl_tx  = plot(ax1, nan,nan,'c','LineWidth',1.2);
hl_rx  = plot(ax1, nan,nan,'m','LineWidth',1.0);
hl_mix = plot(ax1, nan,nan,'y','LineWidth',1.4);
xlabel(ax1,'Time (us)'); ylabel(ax1,'Amplitude');
title(ax1,'Baseband TX / RX / Beat Signal','Color','w');
legend(ax1,{'TX (baseband chirp)','RX (target 1, baseband)','Beat (virt ch 1)'}, ...
       'TextColor','w','Location','northeast');

%% ---- Figure 2: range-Doppler map ----------------------------------------
f2 = figure('Name','Fig 2 - Range-Doppler','Color','k','Position',[620 700 620 300]);
ax2 = axes(f2); style_ax(ax2);
hi_rdm = imagesc(ax2, V_ax, R_ax, zeros(numADC,numChirps));
xlabel(ax2,'Velocity (m/s)'); ylabel(ax2,'Range (m)');
title(ax2,'Range-Doppler Map','Color','w');
cb2 = colorbar(ax2); cb2.Color='w'; cb2.Label.String='dB';
set_clim(ax2,[-60 0]); colormap(ax2,jet(256));

%% ---- Figure 3: CA-CFAR ---------------------------------------------------
f3 = figure('Name','Fig 3 - CA-CFAR','Color','k','Position',[1240 700 620 300]);
ax3 = axes(f3); style_ax(ax3);
hi_cfar = imagesc(ax3, V_ax, R_ax, zeros(numADC,numChirps));
xlabel(ax3,'Velocity (m/s)'); ylabel(ax3,'Range (m)');
title(ax3,'CA-CFAR Detections','Color','w');
colormap(ax3,[0 0 0.3; 1 0.2 0.2]);

%% ---- Figure 4: FFT range-azimuth map ------------------------------------
f4 = figure('Name','Fig 4 - Range-Azimuth','Color','k','Position',[0 360 620 300]);
ax4 = axes(f4); style_ax(ax4);
hi_ra = imagesc(ax4, ang_ax, R_ax, zeros(numADC,numel(ang_ax)));
xlabel(ax4,'Azimuth (deg)'); ylabel(ax4,'Range (m)');
title(ax4,'FFT Range-Azimuth (virtual aperture beamforming)','Color','w');
cb4 = colorbar(ax4); cb4.Color='w'; cb4.Label.String='dB';
set_clim(ax4,[-35 0]); colormap(ax4,jet(256));

%% ---- Figure 5: MUSIC azimuth spectrum -----------------------------------
f5 = figure('Name','Fig 5 - MUSIC Azimuth','Color','k','Position',[620 360 620 300]);
ax5 = axes(f5); style_ax(ax5);
musicColors = lines(MAX_TAR);
hl_music = gobjects(MAX_TAR,1);
for k = 1:MAX_TAR
    hl_music(k) = plot(ax5, ang_ax, nan(size(ang_ax)), ...
        'Color',musicColors(k,:),'LineWidth',1.5);
end
xlabel(ax5,'Azimuth (deg)'); ylabel(ax5,'Power (dB)');
title(ax5,'MUSIC Azimuth Spectrum (MDL rank)','Color','w');
xlim(ax5,[-90 90]); ylim(ax5,[-40 5]);

%% ---- Figure 6: 3-D point cloud ------------------------------------------
f6 = figure('Name','Fig 6 - 3D Point Cloud','Color','k','Position',[1240 360 620 300]);
ax6 = axes(f6); style_ax_3d(ax6);
hpc3 = scatter3(ax6, nan,nan,nan, 60, nan, 'filled');
xlabel(ax6,'X (m)'); ylabel(ax6,'Y (m)'); zlabel(ax6,'Z (m)');
title(ax6,'3D Point Cloud - X/Y/Z with velocity color','Color','w');
xlim(ax6,[-70 70]); ylim(ax6,[0 MAX_DEPTH]); zlim(ax6,[-15 15]);
cb6 = colorbar(ax6); cb6.Color='w'; cb6.Label.String='Velocity (m/s)';
colormap(ax6,jet(256));
view(ax6,-30,25);

%% ---- Figure 7: raw CARLA detections -------------------------------------
f7 = figure('Name','Fig 7 - Raw Detections','Color','k','Position',[0 40 620 280]);
ax7 = axes(f7); style_ax(ax7);
hrv = scatter(ax7, nan,nan, 30, nan, 'filled');
xlabel(ax7,'Radial Velocity (m/s)'); ylabel(ax7,'Depth (m)');
title(ax7,'Raw CARLA Detections - R vs V','Color','w');
xlim(ax7,[-MAX_ABS_VEL MAX_ABS_VEL]); ylim(ax7,[0 MAX_DEPTH]);
cb7 = colorbar(ax7); cb7.Color='w'; cb7.Label.String='Azimuth (deg)';
colormap(ax7,jet(256));

%% ---- Figure 8: info panel ------------------------------------------------
f8 = figure('Name','Fig 8 - Info','Color','k','Position',[620 40 620 280]);
ax8 = axes(f8); axis(ax8,'off');
ht_info = text(ax8,0.05,0.95,'Waiting for CARLA FMCW MIMO data...', ...
    'Color','w','FontSize',11,'Units','normalized', ...
    'VerticalAlignment','top','FontName','Courier New');

%% ---- Figure 9: RX sub-array MUSIC ("elevation" placeholder) -------------
f9 = figure('Name','Fig 9 - RX sub-array MUSIC','Color','k','Position',[1240 40 620 280]);
ax9 = axes(f9); style_ax(ax9);
hl_music_el = gobjects(MAX_TAR,1);
for k = 1:MAX_TAR
    hl_music_el(k) = plot(ax9, el_ax, nan(size(el_ax)), ...
        'Color',musicColors(k,:),'LineWidth',1.5);
end
xlabel(ax9,'Angle (deg)'); ylabel(ax9,'Power (dB)');
title(ax9,'RX sub-array MUSIC - placeholder, cube has no elevation phase','Color','w');
xlim(ax9,[-30 30]); ylim(ax9,[-40 5]);

%% ---- Figure 10: CBF beampattern -----------------------------------------
f10 = figure('Name','Fig 10 - FFT Beampattern','Color','k','Position',[1860 360 700 320]);
ax10 = axes(f10); style_ax(ax10);
hl_fft_spec = gobjects(MAX_TAR,1);
hl_fft_gt   = gobjects(MAX_TAR,1);
for k = 1:MAX_TAR
    hl_fft_spec(k) = plot(ax10, ang_ax, nan(size(ang_ax)), ...
        'Color',musicColors(k,:),'LineWidth',1.8);
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

allFigs    = [f1 f2 f3 f4 f5 f6 f7 f8 f9 f10];

%% ---- Window layout -------------------------------------------------------
% The CARLA window (-ResX=800 -ResY=600 -WinX=0 -WinY=0) sits top-left.
%   rows 1-2, right of CARLA : Fig 1  Fig 4  Fig 3
%                              Fig 2  Fig 5  Fig 7
%   row 3, full width        : Fig 8  Fig 6  Fig 9  Fig 10
% Change CARLA_AREA if the CARLA window on your screen is larger or smaller.
CARLA_AREA = [820 650];          % [width height] reserved for CARLA (pixels)
TASKBAR_H  = 48;                 % Windows taskbar height (pixels)
tile_figures(allFigs, CARLA_AREA, TASKBAR_H);
drawnow;
frameCount = 0;
last_wall_t     = tic;
last_num_chunks = 0;
prebuf          = {};    % datagrams already read but not yet used

%% ---- 2-D Hann window (no Signal Processing Toolbox needed) -------------
win_range   = 0.5 * (1 - cos(2*pi*(0:numADC-1).'    / (numADC-1)));
win_doppler = 0.5 * (1 - cos(2*pi*(0:numChirps-1)   / (numChirps-1)));
win2d       = win_range * win_doppler;     % [numADC x numChirps]

%% =========================================================================
%% MAIN LOOP
%% =========================================================================
while all(isvalid(allFigs))

    t_frame = tic;

    %% Step 1: backlog check ----------------------------------------------
    % If two or more complete frames are already queued, MATLAB has fallen
    % behind: drop the whole input queue and resynchronise on the next
    % META packet. flush_count reports what was discarded for this frame
    % (0 = keeping up in real time).
    flush_count = 0;
    pending = u.NumDatagramsAvailable + numel(prebuf);
    if last_num_chunks > 0 && pending >= 2*(last_num_chunks + 1)
        flush(u, "input");
        prebuf = {};
        flush_count = pending;
    end

    %% Step 2: receive one complete frame (META + all CHUNKs of one ID) ---
    [assembled, num_chunks, dropped, prebuf] = receive_frame(u, prebuf, MAGIC_META, ...
        MAGIC_CHUNK, CHUNK_SIZE, HEADER_BYTES, 1.0);
    t_rx = toc(t_frame);
    flush_count = flush_count + dropped;
    if isempty(assembled)
        drawnow limitrate;
        continue;
    end
    last_num_chunks = num_chunks;

    %% Step 3: parse header (24 bytes) ---------------------------------------
    p = 1;
    pkt_numVirt  = double(typecast(assembled(p:p+3), 'uint32')); p = p+4;
    pkt_numADC   = double(typecast(assembled(p:p+3), 'uint32')); p = p+4;
    pkt_numChirp = double(typecast(assembled(p:p+3), 'uint32')); p = p+4;
    timestamp    =        typecast(assembled(p:p+7), 'double');  p = p+8;
    num_det      = double(typecast(assembled(p:p+3), 'uint32')); p = p+4;

    if pkt_numVirt ~= numVirt || pkt_numADC ~= numADC || pkt_numChirp ~= numChirps
        warning('Dimension mismatch in packet. Skipping.');
        continue;
    end

    %% Step 4: raw CARLA detections ------------------------------------------
    raw_vals = zeros(0,4);
    if num_det > 0
        det_end = p + num_det*16 - 1;
        if det_end > numel(assembled)
            warning('Packet too short for detections. Skipping.');
            continue;
        end
        raw_vals = double(reshape(typecast(assembled(p:det_end), 'single'), 4, []).');
        % columns: [velocity, azimuth_rad, altitude_rad, depth_m]
        p = det_end + 1;
    end

    %% Step 5: IQ radar data cube ---------------------------------------------
    cube_end = p + CUBE_BYTES - 1;
    if cube_end > numel(assembled)
        warning('Packet too short for IQ cube. Skipping.');
        continue;
    end
    cube_single = typecast(assembled(p:cube_end), 'single');
    cube_cplx   = complex(double(cube_single(1:2:end)), double(cube_single(2:2:end)));
    % Python C-order [numVirt, numADC, numChirps] -> MATLAB [numADC x numChirps x numVirt]
    RDC = reshape(cube_cplx, [numChirps, numADC, numVirt]);
    RDC = permute(RDC, [2 1 3]);

    %% Step 6: display filter (same gates as the Python cube filter) --------
    vel_raw = []; dep_raw = []; az_deg = [];
    xp = []; yp = []; zp = [];
    if ~isempty(raw_vals)
        vel_r = raw_vals(:,1);
        az_r  = raw_vals(:,2);
        alt_r = raw_vals(:,3);
        dep_r = raw_vals(:,4);
        x_r   = dep_r .* sin(az_r);
        y_r   = dep_r .* cos(az_r);
        z_r   = dep_r .* sin(alt_r);

        ok = isfinite(x_r) & isfinite(y_r) & isfinite(z_r) & isfinite(vel_r) ...
           & y_r >= 0 ...
           & dep_r >= MIN_DEPTH & dep_r <= MAX_DEPTH ...
           & abs(vel_r) > MIN_ABS_VEL & abs(vel_r) <= MAX_ABS_VEL;

        vel_raw = vel_r(ok);
        dep_raw = dep_r(ok);
        az_deg  = rad2deg(az_r(ok));
        xp = x_r(ok); yp = y_r(ok); zp = z_r(ok);
    end

    if ~isempty(dep_raw)
        set(hrv, 'XData', vel_raw, 'YData', dep_raw, 'CData', az_deg);
        set_clim(ax7, [-90 90]);
    else
        set(hrv, 'XData', nan, 'YData', nan, 'CData', nan);
    end

    %% Step 7: Fig 1 - baseband waveform ------------------------------------
    % The baseband chirp is plotted: at Fs = 25.6 MHz the 77 GHz carrier
    % cannot be represented.
    sig_tx = real(exp(1j*2*pi*(slope/2*t_chirp.^2)));
    if ~isempty(dep_raw)
        td     = t_chirp - 2*dep_raw(1)/c;
        sig_rx = real(exp(1j*2*pi*(slope/2*td.^2)));
    else
        sig_rx = zeros(numADC, 1);
    end
    sig_mix = real(RDC(:,1,1));
    mx = max(abs(sig_mix));
    if mx > 0, sig_mix = sig_mix / mx; end
    t_us = t_chirp * 1e6;
    set(hl_tx,  'XData', t_us, 'YData', sig_tx);
    set(hl_rx,  'XData', t_us, 'YData', sig_rx);
    set(hl_mix, 'XData', t_us, 'YData', sig_mix);
    xlim(ax1, [0 t_us(end)]);

    %% Step 8: AWGN at SNR_dB --------------------------------------------------
    sig_power = mean(abs(RDC(:)).^2);
    if sig_power > 0
        sigma_n = sqrt(sig_power / 10^(SNR_dB/10) / 2);
        RDC = RDC + sigma_n * (randn(size(RDC)) + 1j*randn(size(RDC)));
    end

    %% Step 9: range FFT and Doppler FFT on all virtual channels -------------
    rangeFFT = fft(RDC .* win2d, numADC, 1);             % [ADC x chirp x virt]
    dopFFT   = fftshift(fft(rangeFFT, numChirps, 2), 2); % [range x Doppler x virt]

    %% Step 10: Fig 2 - range-Doppler map (virtual channel 1) ---------------
    RDM    = dopFFT(:,:,1);
    RDM_dB = 20*log10(abs(RDM) ./ (max(abs(RDM(:))) + eps));
    set(hi_rdm, 'CData', RDM_dB);
    title(ax2, sprintf('Range-Doppler Map  |  t = %.2f s', timestamp), 'Color','w');

    %% Step 11: Fig 3 - CA-CFAR ----------------------------------------------
    RDM_pow   = abs(RDM).^2;
    noise_sum = conv2(RDM_pow, W_train, 'same');
    cfar_mask = RDM_pow > (alpha_cfar / Ncfar) * noise_sum;

    border = numTrain + numGuard;            % no complete training window here
    cfar_mask(1:border, :)         = 0;
    cfar_mask(end-border+1:end, :) = 0;
    cfar_mask(:, 1:border)         = 0;
    cfar_mask(:, end-border+1:end) = 0;
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

    %% Step 12: Fig 4 - range-azimuth map at the dominant Doppler bin -------
    [~, dom_dop_bin] = max(max(abs(RDM), [], 1));
    virt_snap = squeeze(dopFFT(:, dom_dop_bin, :));      % [numADC x numVirt]
    ra_map = abs(virt_snap * conj(steer_az));            % [numADC x 181]
    ra_dB  = 20*log10(ra_map ./ (max(ra_map(:)) + eps));
    set(hi_ra, 'CData', ra_dB);
    title(ax4, sprintf('FFT Range-Azimuth  |  %d virtual elements  |  v = %.1f m/s', ...
          numVirt, V_ax(dom_dop_bin)), 'Color','w');

    %% Step 13: Fig 6 - 3-D point cloud ---------------------------------------
    if ~isempty(xp)
        set(hpc3, 'XData', xp, 'YData', yp, 'ZData', zp, 'CData', vel_raw);
        set_clim(ax6, [-MAX_ABS_VEL MAX_ABS_VEL]);
    else
        set(hpc3, 'XData', nan, 'YData', nan, 'ZData', nan, 'CData', nan);
    end
    title(ax6, sprintf('3D Point Cloud  |  %d pts', numel(xp)), 'Color','w');

    %% Step 14: Fig 10 - CBF beampattern per CFAR target --------------------
    for k = 1:MAX_TAR
        set(hl_fft_spec(k), 'YData', nan(size(ang_ax)));
        set(hl_fft_gt(k),   'Value', nan);
    end

    numTar = min(K, MAX_TAR);
    for k = 1:numTar
        x_snap   = squeeze(dopFFT(cfar_r(k), cfar_d(k), :));   % [numVirt x 1]
        bf_power = abs(steer_az' * x_snap).^2;
        bf_dB    = 10*log10(bf_power ./ (max(bf_power) + eps));
        set(hl_fft_spec(k), 'YData', bf_dB);

        % CARLA reference azimuth: filtered detection closest in range and velocity
        if ~isempty(dep_raw)
            cost = abs(dep_raw - R_ax(cfar_r(k)))/dR + abs(vel_raw - V_ax(cfar_d(k)))/dV;
            [~, nearest_det] = min(cost);
            set(hl_fft_gt(k), 'Value', az_deg(nearest_det));
        end
    end
    title(ax10, sprintf('FFT Beampattern  |  %d target(s)  |  peak = 0 dB', numTar), ...
          'Color','w');

    %% Step 15: MUSIC azimuth (8-element virtual array) + RX sub-array -------
    for k = 1:MAX_TAR
        set(hl_music(k),    'YData', nan(size(ang_ax)));
        set(hl_music_el(k), 'YData', nan(size(el_ax)));
    end

    tx0_chirps = 1:numTX:numChirps;
    for k = 1:numTar
        ri = cfar_r(k);

        % Azimuth: snapshots = chirps of range bin ri  [numVirt x numChirps]
        A_az   = squeeze(rangeFFT(ri, :, :)).';
        Rxx_az = (A_az * A_az') / numChirps;
        Rxx_fb = 0.5 * (Rxx_az + J_virt * conj(Rxx_az) * J_virt);

        [Q_az, D_az] = eig(Rxx_fb);
        [ev_az, order] = sort(real(diag(D_az)), 'descend');
        Q_az   = Q_az(:, order);
        d_mdl  = min(estimate_rank_mdl(ev_az, numChirps), numVirt - 2);
        Qn_az  = Q_az(:, d_mdl+1:end);
        music_az = music_spectrum(steer_az, Qn_az);
        set(hl_music(k), 'YData', music_az - max(music_az));

        % RX sub-array (TX0 chirps, virtual elements 1..4) - see note in header
        rx_snap = squeeze(rangeFFT(ri, tx0_chirps, 1:numRX)).';   % [numRX x numChirps/2]
        Rxx_el  = (rx_snap * rx_snap') / size(rx_snap, 2);
        [Q_el, D_el] = eig(Rxx_el);
        [ev_el, order] = sort(real(diag(D_el)), 'descend');
        Q_el     = Q_el(:, order);
        d_mdl_el = min(estimate_rank_mdl(ev_el, size(rx_snap,2)), numRX - 2);
        Qn_el    = Q_el(:, d_mdl_el+1:end);
        music_el = music_spectrum(steer_el, Qn_el);
        set(hl_music_el(k), 'YData', music_el - max(music_el));
    end
    title(ax5, sprintf('MUSIC Azimuth  |  %d targets (MDL)', numTar), 'Color','w');

    %% Step 16: Fig 8 - info panel ------------------------------------------
    noise_floor_db = 10*log10(mean(abs(RDC(:)).^2) + eps);
    fps = 1 / toc(last_wall_t);
    last_wall_t = tic;

    update_info(ht_info, timestamp, numel(vel_raw), K, numTar, frameCount, ...
                dR, dV, Vmax, numVirt, noise_floor_db, SNR_dB, flush_count);

    frameCount = frameCount + 1;
    t_dsp = toc(t_frame) - t_rx;
    t_draw0 = tic;
    drawnow limitrate;
    t_draw = toc(t_draw0);

    fprintf(['Frame %3d | t=%7.2f | pts=%3d | CFAR=%2d | MUSIC=%d | flushed=%d | %.1f fps' ...
             ' | rx %4.0f ms  dsp %4.0f ms  draw %4.0f ms\n'], ...
            frameCount, timestamp, numel(vel_raw), K, numTar, flush_count, fps, ...
            1000*t_rx, 1000*t_dsp, 1000*t_draw);

    % Display cap: idle until the next display slot, then drop the frames
    % that arrived meanwhile (intentional, not counted as backlog).
    idle_s = 1/DISPLAY_FPS - toc(t_frame);
    if idle_s > 0
        pause(idle_s);
        flush(u, "input");
        prebuf = {};
    end
end

disp('A figure was closed.  Radar stopped.');
clear u;

%% =========================================================================
%% LOCAL FUNCTIONS
%% =========================================================================

function [assembled, num_chunks, dropped, prebuf] = receive_frame(u, prebuf, MAGIC_META, ...
        MAGIC_CHUNK, CHUNK_SIZE, HEADER_BYTES, timeout_s)
% Reassemble one frame: a META packet followed by all CHUNK packets that
% carry the same frame_id. Datagrams are read from the port in batches
% (one read call for everything queued); unused ones are returned in
% 'prebuf' for the next frame. Returns [] on timeout. 'dropped' counts the
% datagrams discarded on the way (stale chunks, chunks of another frame,
% duplicates, frames abandoned because a newer META arrived).
    assembled = [];
    num_chunks = 0;
    dropped = 0;
    have_meta = false;
    t0 = tic;

    while toc(t0) < timeout_s
        if isempty(prebuf)
            n_avail = u.NumDatagramsAvailable;
            if n_avail <= 0
                pause(0.001);
                continue;
            end
            dgs    = read(u, n_avail, "uint8");
            prebuf = {dgs.Data};
        end
        dgd = uint8(prebuf{1});
        prebuf(1) = [];
        if numel(dgd) < 12
            dropped = dropped + 1;
            continue;
        end
        magic = typecast(dgd(1:4), 'uint32');

        if magic == MAGIC_META
            n_c = double(typecast(dgd(5:8),  'uint32'));
            tot = double(typecast(dgd(9:12), 'uint32'));
            if have_meta
                dropped = dropped + sum(received);     % abandon the unfinished frame
            end
            if tot < HEADER_BYTES || n_c ~= ceil(tot / CHUNK_SIZE)
                have_meta = false;
                dropped = dropped + 1;
                continue;
            end
            have_meta  = true;
            num_chunks = n_c;
            buf        = zeros(tot, 1, 'uint8');
            received   = false(n_c, 1);
            frame_id   = [];
            t0 = tic;                                  % fresh timeout for this frame

        elseif magic == MAGIC_CHUNK && have_meta && numel(dgd) >= 16
            fid   = typecast(dgd(5:8),   'uint32');
            c_idx = double(typecast(dgd(9:12),  'uint32')) + 1;
            c_len = double(typecast(dgd(13:16), 'uint32'));
            if isempty(frame_id)
                frame_id = fid;                        % lock onto this frame
            end
            s_byte = (c_idx - 1) * CHUNK_SIZE + 1;
            e_byte = s_byte + c_len - 1;
            if fid ~= frame_id || c_idx < 1 || c_idx > num_chunks || c_len < 1 ...
                    || numel(dgd) < 16 + c_len || e_byte > numel(buf) || received(c_idx)
                dropped = dropped + 1;
                continue;
            end
            buf(s_byte:e_byte) = dgd(17:16+c_len);
            received(c_idx) = true;
            if all(received)
                assembled = buf;
                return;
            end

        else
            dropped = dropped + 1;                     % chunk without META, or junk
        end
    end
end

function tile_figures(figs, carla_area, taskbar_h)
% Place the 10 figure windows around the CARLA window on the first monitor.
% figs = [f1 ... f10]; positions are outer positions (incl. title bar).
    mon = get(groot, 'MonitorPositions');
    scr = mon(1, :);                                  % [left bottom width height]
    x0 = scr(1);  y0 = scr(2) + taskbar_h;
    W  = scr(3);  H  = scr(4) - taskbar_h;
    rowH = floor(H / 3);
    yTop = y0 + 2*rowH;  yMid = y0 + rowH;  yBot = y0;

    if W - carla_area(1) >= 900 && 2*rowH >= carla_area(2) - 40
        % Rows 1-2: three columns to the right of CARLA
        xs = x0 + carla_area(1);
        colW = floor((W - carla_area(1)) / 3);
        top = [1 4 3];  mid = [2 5 7];
        for i = 1:3
            set(figs(top(i)), 'OuterPosition', [xs + (i-1)*colW, yTop, colW, rowH]);
            set(figs(mid(i)), 'OuterPosition', [xs + (i-1)*colW, yMid, colW, rowH]);
        end
        % Row 3: four columns across the full width
        bot = [8 6 9 10];  colW4 = floor(W / 4);
        for i = 1:4
            set(figs(bot(i)), 'OuterPosition', [x0 + (i-1)*colW4, yBot, colW4, rowH]);
        end
    else
        % Small screen: plain 4 x 3 grid (CARLA may overlap Fig 1)
        colW4 = floor(W / 4);
        order = [1 4 3 2 5 7 8 6 9 10];
        for k = 1:numel(order)
            r = floor((k-1) / 4);  c = mod(k-1, 4);
            set(figs(order(k)), 'OuterPosition', ...
                [x0 + c*colW4, y0 + (2-r)*rowH, colW4, rowH]);
        end
    end
end

function P_dB = music_spectrum(steer, Qn)
% MUSIC pseudo-spectrum in dB for every column of the steering matrix.
    proj  = Qn' * steer;                               % [noise dim x angles]
    denom = sum(abs(proj).^2, 1);
    P_dB  = 10*log10(real(sum(abs(steer).^2, 1)) ./ (denom + eps));
    P_dB  = P_dB(:);
end

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
        ' CFAR det       : %d (after cluster suppression)\n', ...
        ' MUSIC targets  : %d\n', ...
        ' Flushed pkts   : %d\n', ...
        '\n', ...
        ' --- FMCW MIMO ---\n', ...
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
% Keep one detection per cluster of nearby range-Doppler bins.
    keep = true(numel(r_in), 1);
    for i = 2:numel(r_in)
        for j = 1:i-1
            if keep(j) && abs(r_in(i)-r_in(j)) <= thr && abs(d_in(i)-d_in(j)) <= thr
                keep(i) = false;
                break;
            end
        end
    end
    cr = r_in(keep);
    cd = d_in(keep);
end

function k_hat = estimate_rank_mdl(eigenvalues, N)
% MDL model-order estimator (Wax & Kailath 1985).
% eigenvalues sorted descending, N = number of snapshots.
    L   = numel(eigenvalues);
    mdl = zeros(L-1, 1);
    for k = 0:L-2
        noise_eigs = max(eigenvalues(k+1:end), eps);
        ratio      = exp(mean(log(noise_eigs))) / (mean(noise_eigs) + eps);
        mdl(k+1)   = -N*(L-k)*log(ratio) + 0.5*k*(2*L-k)*log(N);
    end
    [~, idx] = min(mdl);
    k_hat = max(1, min(idx, L-1));
end
