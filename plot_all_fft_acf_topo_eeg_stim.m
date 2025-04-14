
clear all
close all 


%% SETTINGS - Set Output Directory
output_dir = fullfile(pwd, 'MATLAB_exports');
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

fig_dir = '/Users/florarosenzweig/Documents/MATLAB/PROJECTS/XP_Audio_DutyCycles/Global_figures';

%add the necessary library folders to MATLAB path (make 
% sure you change these to your local paths)
addpath(genpath('/Users/florarosenzweig/Documents/MATLAB/Toolboxes/rnb_tools'));
addpath(genpath('/Users/florarosenzweig/Documents/MATLAB/Toolboxes/acf_tools'));
Paths.lw6ElecLoc = 'Users/florarosenzweig/Documents/MATLAB/Toolboxes/letswave6-master/core_functions/Standard-10-20-Cap81.locs';

% add letswave6 for this one 
% https://github.com/NOCIONS/letswave6
addpath(genpath('/Users/florarosenzweig/Documents/MATLAB/Toolboxes/letswave6-master')); 

% small utilities here 
addpath(genpath('./lib')); 

% we will need some EEG data (download it from onedrive)
% https://uclouvain.sharepoint.com/:f:/r/sites/O365G-RBLab/Documents%20partages/RnB-LAB_People/TLenc/acf_workshop?csf=1&web=1&e=0Do0e2
eeg_path =  './EEG_data'; 

cond = cell(1, 4);
cond{1} = '25% Duty Cycle';
cond{2} = '50% Duty Cycle';
cond{3} = '75% Duty Cycle';
cond{4} = '100% Duty Cycle';

% we will need some EEG data (download it from onedrive)
% https://uclouvain.sharepoint.com/:f:/r/sites/O365G-RBLab/Documents%20partages/RnB-LAB_People/TLenc/acf_workshop?csf=1&web=1&e=0Do0e2
eeg_path =  './EEG_data'; 
eeg_fname = cell(1, 4);

eeg_fname{1} = 'chan-select merged_YOUNG_ACF_EEG_25.lw6'; 
eeg_fname{2} = 'chan-select merged_YOUNG_ACF_EEG_50.lw6'; 
eeg_fname{3} = 'chan-select merged_YOUNG_ACF_EEG_75.lw6'; 
eeg_fname{4} = 'chan-select merged_YOUNG_ACF_EEG_100.lw6'; 


%Plot preprocessed data directly not using the script
eeg_fft = cell(1, 4);
eeg_fft{1} = 'grand_average reref_M1M2_EEG_25.lw6';
eeg_fft{2} = 'grand_average reref_M1M2_EEG_50.lw6';
eeg_fft{3} = 'grand_average reref_M1M2_EEG_75.lw6';
eeg_fft{4} = 'grand_average reref_M1M2_EEG_100.lw6';


minall = [];
maxall = [];

% and some cochlear model data (onedrive)
coch_path = eeg_path; 
coch_fname = cell(1, 4);
coch_fname{1} = '50ms_40.80s_300Hz_syncopXP_stimuli_AN.mat'; 
coch_fname{2} = '100ms_40.80s_300Hz_syncopXP_stimuli_AN.mat'; 
coch_fname{3} = '150ms_40.80s_300Hz_syncopXP_stimuli_AN.mat'; 
coch_fname{4} = '200ms_40.80s_300Hz_syncopXP_stimuli_AN.mat'; 


%% Define Colors for Beat-Related & Unrelated Frequencies
col_meter_rel = [0.8706, 0.1765, 0.1490];  % Red
col_meter_unrel = [0.1922, 0.5098, 0.7412];  % Blue
col_neutral = [186, 186, 186] / 255; % Grey

%duty_cycle_labels = {'25', '50', '75', '100'}; % Labels for each row


globalFontSize = 35;  % Set a consistent font size

%% Initialize Min/Max Storage for Topoplot Scaling
minall = [];
maxall = [];

%%

trial_dur = 40.8;  

%% -------------------- FIRST LOOP: Compute Min/Max for Topoplots --------------------
for i = 1:length(coch_fname)

    % Load EEG FFT Data
    [headerFFT, dataFFT] = CLW_load(fullfile(eeg_path, eeg_fft{i}));
    dataFFT = squeeze(dataFFT);
    
    % Get electrode locations
    locs  = readlocs(Paths.lw6ElecLoc);
    for iLocs = 1:length(locs)
        locs(iLocs).topo_enabled = 1;
    end
    [headerFFT] = RLW_edit_electrodes(headerFFT, locs);
    
    % Compute Sample Step & Frequency Vector
    binFFT = headerFFT.xstep;
    maxfreqlim = 5.2;
    freqFFT = [0:binFFT:maxfreqlim];

    % Define Meter-Related and Unrelated Frequencies
    frex = (1/2.4) * (1:maxfreqlim*2.4); % Harmonic Series
    frexidx = round(frex / binFFT) + 1;

    %%%% Vector of frequencies up to limit Hz
    frex    = 1/2.4*(1:maxfreqlim*2.4); % all frequencies
    frexidx = round(frex/headerFFT.xstep)+1;

%%%% Define meter and unrelated frequencies and indices
    nHarm = 1; %1 = 5Hz, 2 = 10Hz, 3 = 15Hz, 4=20Hz, 5=25Hz, 6=30Hz, etc
    
    % Define the range based on a multiple of 12
    maxFreq = nHarm * 12; % This assumes nHarm is the multiple of 12 you're interested in
    allfrex    = 1/2.4*(1:maxFreq); % get all freqs in Hz

    %Define which freqs will not be taken into consideration
    excludeBeg = [1 2]; % If want to exclude 0.416, 0.832 from both sets
    exclude5Hzharms = [12:12:maxFreq];
    excludeAll = [excludeBeg exclude5Hzharms];

    % Define Meter Related frequencies
    allMeterRel = 3:3:maxFreq;
    frexMeterRel= setdiff(allMeterRel, exclude5Hzharms);
    frexMR      = frex(frexMeterRel);
    frexMRidx   = frexidx(frexMeterRel);


    %Define Meter Unrelated frequencies
    frexMeterUnRel = setdiff(1:maxFreq,[frexMeterRel excludeAll]);
    frexMUR        = frex(frexMeterUnRel);
    frexMURidx     = frexidx(frexMeterUnRel);


    % Compute Mean Values for Beat-Related & Unrelated Frequencies
    topoMetRel = mean(dataFFT(:, frexMRidx), 2);
    topoMetUnrel = mean(dataFFT(:, frexMURidx), 2);

    % Store Min/Max Values
    minall = [minall; min(topoMetRel)];
    maxall = [maxall; max(topoMetUnrel)];
end

% Set Global Color Limits for Topoplot
minTopo = min(minall);
maxTopo = max(maxall);

%% -------------------- SECOND LOOP: Plot FFT, ACF, and Topoplots --------------------
figure('Color', 'w', 'Position', [-1254 924 1812 976], 'Units', 'pixels');
%figure('Color', 'w', 'Position', [-1254 924 1812 1320], 'Units', 'pixels');
nrows = 4;
ncols = 4;
ax = gobjects(nrows * ncols, 1);

for k = 1:(nrows * ncols)
    ax(k) = subplot(nrows, ncols, k);  % Create all axes
end

%reposition_subplots(ax, 4, 4, 0.07, 0.1, 0.2, 0.18, 0.03, 0.1);
reposition_subplots(ax, 4, 4, 0.07, 0.07, 0.22, 0.16, 0.03, 0.105);

%(ax, rows, cols, hgap, vgap, top, left, right, bottom)

for i = 1:length(coch_fname)

    % -------------------- Load EEG FFT --------------------
    [headerFFT, dataFFT] = CLW_load(fullfile(eeg_path, eeg_fft{i}));
    dataFFT = squeeze(dataFFT);
    
    %%%% Frequencies of interest
    max_freq = 5; 
    f0_to_excl = 5; 

    [freq_meter_rel, freq_meter_unrel, frex] = get_meter_freq(...
                                                    max_freq, ...
                                                    'f0_to_excl', f0_to_excl);


    % -------------------- Load EEG ACF --------------------
    [header, data] = CLW_load(fullfile(eeg_path, eeg_fname{i}));
    
    fs = 1/header.xstep; 
    t = [0 : header.datasize(end)-1] * header.xstep + header.xstart; 
    
    %%%% Lags of interest (EEG and cochlear model output should have same duration) 
    max_lag = trial_dur/ 2; 

    lag_base_incl_meter_rel = [0.8]; 
    lag_base_excl_meter_rel = [0.6, 1.0, 1.4]; % [0.6, 1.0, 1.4]   [2.4]

    lag_base_incl_meter_unrel = [0.6, 1.0, 1.4]; % [0.6, 1.0, 1.4]   [0.2]
    lag_base_excl_meter_unrel = [0.4, 0.8]; 

    [lags_meter_rel, lags_meter_unrel] = get_meter_lags(...
            max_lag, ...
            lag_base_incl_meter_rel, lag_base_excl_meter_rel, ...
            lag_base_incl_meter_unrel, lag_base_excl_meter_unrel ...
            );

    % -------------------- Load Stimulus --------------------
    coch_output = load(fullfile(coch_path, coch_fname{i}));
    %%%% Extract firing rates from all channels
        cochN   = size(coch_output.AN.an_sout, 2);       
        coch    = squeeze(sum(coch_output.AN.an_sout,1)); %mean rate
        fs_coch = coch_output.AN.Fs_an_sout;   
        N       = size(coch_output.AN.an_sout, 2); 
        t_coch  = [0 : N-1] / fs_coch; 

    %%%% Low-pass filter for nicer acf plot 
        [b,a] = butter(2, 30/(fs_coch/2), 'low'); 
        coch = filtfilt(b, a, coch); 

    %%%%  confirm we have integer number of cycles 
    assert(mod(length(coch) / fs_coch, 2.4) == 0); 

    % -------------------- Compute FFT & ACF from cochlear model --------------------
    %%%%  get ACF (withuout aperiodic subtraction)
    [acf_coch, lags_coch, ~, mX_coch, freq_coch] = get_acf(coch, fs_coch);    

    %%%%  get ACF features
    feat_acf_coch = get_acf_features(...
                                acf_coch, lags_coch, ...
                                lags_meter_rel, lags_meter_unrel);    

    %%%%  get FFT features for the raw spectra                                    
    feat_fft_coch = get_fft_features(mX_coch, freq_coch, ...
                                     freq_meter_rel, freq_meter_unrel); 
    % -------------------- Compute FFT & ACF from EEG --------------------
    %%%%  without 1/f subtraction
    [acf, lags, ~, mX, freq] = get_acf(data, fs);
    
    %%%% with 1/f subtraction
    acf_subtracted_step2 = nan(size(acf));
    acf_subtracted_step3 = nan(size(acf));
    
    
    parfor i_sub=1:size(data, 1)
        
        fprintf('sub-%02d\n', i_sub);
        
        [acf_subtracted_step2(i_sub, :, 1, 1, 1, :)] = ...
            get_acf(data(i_sub, :, 1, 1, 1, :), fs, ...
            'rm_ap', true, ...
            'ap_fit_method', 'irasa', ...
            'plot_diagnostic', false);
        
        [acf_subtracted_step3(i_sub, :, 1, 1, 1, :)] = ...
            get_acf(data(i_sub, :, 1, 1, 1, :), fs, ...
            'rm_ap', true, ...
            'ap_fit_method', 'irasa', ...
            'f0_to_ignore', 1 / 2.4, ...
            'plot_diagnostic', false);
        
    end
    
    all_acf_subtracted_step3 = mean(squeeze(acf_subtracted_step3),1);
    
    mX_subtracted_step2 = subtract_noise_bins(mX, 2, 5);
    
    
    %%%% for illustration, get noise-subtracte FFT and zero out noise bins (step 3)
    frex_signal = [1/2.4 : 1/2.4 : freq(end)];
    frex_signal_idx = dsearchn(freq', frex_signal')';
    frex_noise_idx = setdiff([1:length(freq)], frex_signal_idx);
    
    mX_subtracted_step3 = mX_subtracted_step2;
    %mX_subtracted_step3(:, :, 1, 1, 1, frex_noise_idx) = 0; %%keep this ???
    all_mX_subtracted_step3 = mean(squeeze(mX_subtracted_step3),1);
    
    %%%%  get ACF features
    feat_acf_eeg_raw = get_acf_features(acf, lags, ...
        lags_meter_rel, lags_meter_unrel);
    
    feat_acf_eeg_subtr_step2 = get_acf_features(acf_subtracted_step2, lags, ...
        lags_meter_rel, lags_meter_unrel);
    
    feat_acf_eeg_subtr_step3 = get_acf_features(acf_subtracted_step3, lags, ...
        lags_meter_rel, lags_meter_unrel);
    
    %%%% get FFT features
    feat_fft_eeg_raw = get_fft_features(mX, freq, ...
        freq_meter_rel, freq_meter_unrel);
    
    feat_fft_eeg_subtr_step2 = get_fft_features(mX_subtracted_step2, freq, ...
        freq_meter_rel, freq_meter_unrel);
    
    feat_fft_eeg_subtr_step3 = get_fft_features(mX_subtracted_step3, freq, ...
        freq_meter_rel, freq_meter_unrel);
    
    eeg_z_snr = get_z_snr(mX, freq, frex, 2, 5);

    % -------------------- Compute Topoplots --------------------
    topoMetRel   = mean(dataFFT(:, frexMRidx), 2);
    topoMetUnrel = mean(dataFFT(:, frexMURidx), 2);
    
        %% ----------------- PLOT EEG FFT WITH TOPOPLOT -----------------
    %ax_fft = subplot(4, 4, (i-1)*4 + 1); % Position in grid
    ax_fft      = ax((i-1)*4 + 1);
    cla(ax_fft);  % Clear axis
    hold(ax_fft, 'on');

   s1 = stem(ax_fft, freqFFT, dataFFT(67, 1:length(freqFFT)), 'marker', 'none', 'LineWidth', 3, 'color', col_neutral);
   s2 = stem(ax_fft, frexMR, dataFFT(67, frexMRidx), 'marker', 'none', 'LineWidth', 3, 'color', col_meter_rel);
   s3 = stem(ax_fft, frexMUR, dataFFT(67, frexMURidx), 'marker', 'none', 'LineWidth', 3, 'color', col_meter_unrel);
 
    
    ax_fft.Box = 'off';
    ax_fft.YLim = [-0.001, 0.5];
    ax_fft.YTick = [0, 0.4];
    ax_fft.TickDir = 'out';
    if i == 4
        ax_fft.XLabel.String = 'Frequency (Hz)';
        ax_fft.XTick = [0, 5];
    else
        ax_fft.XTick = [];
    end
    
 % Add Topoplot as Inset
    colormap(parula);
    
% Get current subplot position (relative to figure)
fft_pos = get(ax_fft, 'Position');

% === HARD-CODED INSET SIZE AND SPACING (independent of subplot size)
% Only add insets if we're in column 1

    pos = get(ax((i-1)*4 + 1), 'Position');

    % === Hardcoded dimensions and spacing
    inset_width = 0.1;
    inset_height = 0.1;
    inset_overlap = 0.04;
    inset_x_offset = 0.03;
    inset_y_offset = 0.06;

    x0 = pos(1) + inset_x_offset;
    y0 = pos(2) + pos(4) - inset_height + inset_y_offset;

    inset_pos_rel = [x0, y0, inset_width, inset_height];
    inset_pos_unrel = [x0 + inset_width - inset_overlap, y0, inset_width, inset_height];
    
    % === Topoplot Beat-Related
ax_topoplot_rel = axes('Position', inset_pos_rel);
CLW_topoplot_vector(headerFFT, topoMetRel, ...
    'style', 'map', ...
    'maplimits', [0 maxTopo + 0.1]);
axis off;
axis equal;

% Overlay axis for ring
ax_ring_rel = axes('Position', inset_pos_rel, 'Color', 'none');
axis(ax_ring_rel, 'off');
axis(ax_ring_rel, 'equal');
xlim(ax_ring_rel, [-0.75 0.75]);
ylim(ax_ring_rel, [-0.75 0.75]);
hold(ax_ring_rel, 'on');

theta = linspace(0, 2*pi, 100);
r = 0.7;
plot(ax_ring_rel, r*cos(theta), r*sin(theta), ...
     'Color', col_meter_rel, 'LineWidth', 3);


% === Topoplot Beat-Unrelated
ax_topoplot_unrel = axes('Position', inset_pos_unrel);
CLW_topoplot_vector(headerFFT, topoMetUnrel, ...
    'style', 'map', 'maplimits', [0 maxTopo + 0.1]);
axis off;
axis equal;

% Overlay axis for ring
ax_ring_unrel = axes('Position', inset_pos_unrel, 'Color', 'none');
axis(ax_ring_unrel, 'off');
axis(ax_ring_unrel, 'equal');
xlim(ax_ring_unrel, [-0.75 0.75]);
ylim(ax_ring_unrel, [-0.75 0.75]);
hold(ax_ring_unrel, 'on');

theta = linspace(0, 2*pi, 100);
r = 0.7;
plot(ax_ring_unrel, r*cos(theta), r*sin(theta), ...
     'Color', col_meter_unrel, 'LineWidth', 3);
    
% % Create axes
% ax_topoplot_rel = axes('Position', inset_pos_rel);
% hold on;
% 
% 
% 
% % Plot topoplot
% CLW_topoplot_vector(headerFFT, topoMetRel, ...
%     'style', 'map', ...
%     'maplimits', [0 maxTopo + 0.1]);
%     axis off;
%     axis equal;
%     axis([-0.75 0.75 -0.75 0.75]);  % ensures consistent size across plots
%     
%     % Draw ring
%     hold on;
%     theta = linspace(0, 2*pi, 100);
%     r = 0.53;  % outer head radius from topoplot
%     plot(r * cos(theta), r * sin(theta), 'Color', col_meter_rel, 'LineWidth', 4);
% 
%     
%     % === Beat-Unrelated topoplot (with colored ring)
%     ax_topoplot_unrel = axes('Position', inset_pos_unrel);
%     CLW_topoplot_vector(headerFFT, topoMetUnrel, ...
%         'style', 'map', 'maplimits', [0 maxTopo + 0.1]);
%     axis off;
%     axis equal;
%     axis([-0.75 0.75 -0.75 0.75]);  % ensures consistent size across plots
%     
%     % Draw ring
%     hold on;
%     theta = linspace(0, 2*pi, 100);
%     r = 0.53;  % outer head radius from topoplot
%     plot(r * cos(theta), r * sin(theta), 'Color', col_meter_unrel, 'LineWidth', 4);
% 
%     
% 

    % === LABELS ABOVE INSETS (only in top row)
    if i == 1
        
 

        % === Colorbar only once
        colormap(ax_topoplot_unrel, parula);

        cb_x = inset_pos_unrel(1) + inset_width-0.02;
        cb_y = inset_pos_unrel(2)+0.01;
        cb_width = 0.008;
        cb_height = inset_height * 1.1;

        cb = colorbar(ax_topoplot_unrel, 'Position', [cb_x, cb_y, cb_width, cb_height]); 
        cb.Ticks = [0, maxTopo + 0.1];
        cb.TickLabels = {'0', sprintf('%.1f μV', maxTopo + 0.1)};
        uistack(cb, 'top');
    end

% Dynamically Position Colorbar
drawnow; % Force MATLAB to process everything before adding insets



    %% ----------------- PLOT STIMULUS FFT -----------------
    %ax_fft_stim = subplot(4, 4, (i-1)*4 + 2);
    ax_fft_stim = ax((i-1)*4 + 2);

    cla(ax_fft_stim);  % Clear axis
    hold on;
    
    
    
    % Perform the FFT plot
    plot_fft(freq_coch, mX_coch, ...
        'ax', ax_fft_stim , ...
        'frex_meter_rel', freq_meter_rel, ...
        'frex_meter_unrel', freq_meter_unrel, ...
        'maxfreqlim', max_freq);
    
    ax_fft_stim.YTick = [0, 1];
    ax_fft_stim.YLim = [0, 1];
    if i == 4
        ax_fft_stim.XLabel.String = 'Frequency (Hz)';
        ax_fft_stim.XTick = [0, 5];
    else
        ax_fft_stim.XTick = [];
    end
    
    



    %% ----------------- PLOT EEG ACF -----------------
   
   ax_acf_eeg  = ax((i-1)*4 + 3);

    cla(ax_acf_eeg);

    hold on;
    max_lag_plot = 2.4; 

    plot_acf(ax_acf_eeg, ...
         all_acf_subtracted_step3 (1, :), ...
         lags, ...
         'lags_meter_rel', lags_meter_rel, ...
         'lags_meter_unrel', lags_meter_unrel, ...
         'max_lag', max_lag_plot, ...
         'linew', 2, ...
         'opacity_lagz', 1, ...
         'prec', 1e6); 

    ax_acf_eeg.YTick = [];
    if i == 4
        ax_acf_eeg.XLabel.String = 'Lag (s)';
        ax_acf_eeg.XTick = [0:0.8:2.4];
    else
        ax_acf_eeg.XTick = [];
    end


    %% ----------------- PLOT STIMULUS ACF -----------------
    %ax_acf_stim = subplot(4, 4, (i-1)*4 + 4);
    ax_acf_stim = ax((i-1)*4 + 4);

   
    cla(ax_acf_stim);  % Clear axis
    hold on;
    
    plot_acf(ax_acf_stim, ...
         acf_coch, ...
         lags_coch, ...
         'lags_meter_rel', lags_meter_rel, ...
         'lags_meter_unrel', lags_meter_unrel, ...
         'max_lag', max_lag_plot, ...
         'linew', 2, ...
         'opacity_lagz', 1, ...
         'prec', 1e6); 
ax_acf_stim.YTick = []; 

if i == 4
     ax_acf_stim.XLabel.String = 'Lag (s)';
     ax_acf_stim.XTick = [0:0.8:2.4];
else
        ax_acf_stim.XTick = [];
end

end



%% ----------------- ADD LEGEND -----------------


% % --------- GLOBAL LABELS (OUTSIDE LOOP) --------- %


% Vertical column labels (e.g., next to each column of subplots)
col_labels = {'Magnitude (µV)', 'Magnitude (a.u.)', 'Magnitude (a.u.)', 'Magnitude (a.u.)'};

x_offsets = zeros(1, 4);
y_centers = zeros(1, 4);

for c = 1:4
    col_idxs = c:4:16;  % Get all subplot indices for this column

    pos_top = get(ax(col_idxs(1)), 'Position');        % Top subplot in column
    pos_bottom = get(ax(col_idxs(end)), 'Position');   % Bottom subplot

    % Center Y between top and bottom
    y_top = pos_top(2) + pos_top(4);
    y_bot = pos_bottom(2);
    y_centers(c) = y_bot + (y_top - y_bot) / 2;

    % X offset: more space for column 1, less for others
    switch c
        case 1
            x_offsets(c) = pos_top(1) - 0.05;  % further left
        case 2
            x_offsets(c) = pos_top(1) - 0.04;  % closer
        otherwise
            x_offsets(c) = pos_top(1) - 0.03;  % closer
    end
end

% Create invisible overlay axes
label_axes = axes('Position', [0 0 1 1], 'Units', 'normalized', ...
                  'Visible', 'off', 'HitTest', 'off');

% Add vertical rotated text labels
for c = 1:4
    text('Parent', label_axes, ...
         'Units', 'normalized', ...
         'Position', [x_offsets(c), y_centers(c)], ...
         'String', col_labels{c}, ...
         'FontSize', 24, ...
         'Rotation', 90, ...
         'HorizontalAlignment', 'center', ...
         'VerticalAlignment', 'middle');
end


% % 🔹 Add global horizontal labels

% Get top subplot positions
x_centers = zeros(1, 4);
for c = 1:4
    pos = get(ax(c), 'Position');
    x_centers(c) = pos(1) + pos(3)/2;
end

% Fixed Y for header and sublabels
y_header = 0.94;
y_sub = 0.90;

% Main group labels
group_labels = {'Magnitude Spectrum', '', 'Autocorrelation', ''};
for c = [1, 3]  % Only col groups 1-2 and 3-4
    annotation('textbox', [mean(x_centers(c:c+1)) - 0.15, y_header, 0.3, 0.05], ...
        'String', group_labels{c}, ...
        'FontSize', 24, 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center', ...
        'EdgeColor', 'none');
end

% Sub-column labels
sub_labels = {'EEG', 'Auditory Nerve Model', 'EEG', 'Auditory Nerve Model'};
for c = 1:4
    string{c} = annotation('textbox', [x_centers(c) - 0.05, y_sub, 0.1, 0.04], ...
        'String', sub_labels{c}, ...
        'FontSize', 24, ...
        'HorizontalAlignment', 'center', ...
        'EdgeColor', 'none');
end




% Make all text the same except inset labels
set(findall(gcf, '-property', 'FontSize'), 'FontSize', globalFontSize);
if exist('h1', 'var')
    set([h1, h2], 'FontSize', 20);
end


% Add legend to a separate subplot
% Create an invisible axes to hold the legend (without affecting any subplot)
legend_axes = axes('Position', [0 1 1 1], 'Visible', 'off'); % Full figure normalized space
hold(legend_axes, 'on');

% Create dummy plots for legend entries (these won't appear in the figure)
h1 = plot(NaN, NaN, 'LineWidth', 4, 'Color', col_meter_rel);
h2 = plot(NaN, NaN, 'LineWidth', 4, 'Color', col_meter_unrel);

% Add the legend at the exact position you specified
lgd = legend(legend_axes, [h1, h2], {'Beat-Related', 'Beat-Unrelated'}, ...
    'FontSize', 30, 'Box', 'off');

% Force MATLAB to use the exact positioning you provided
set(lgd, 'Position', [0.013021378572025,0.903000901526686,0.135485651214128,0.079223504721931]);


% --------- CONDITION ROW LABELS  --------- %

row_labels = {'25%', '50%', '75%', '100%'};

for i = 1:4
    ax_idx = (i-1)*4 + 1;  % First column of each row
    pos = get(ax(ax_idx), 'Position');  % [x y w h]
    y_center = pos(2) + pos(4)/2;

    annotation('textbox',...
        [0.01, y_center - 0.025, 0.06, 0.05],... % left x, center y
        'VerticalAlignment','middle',...
        'String',{row_labels{i}, 'Duty', 'Cycle'},...
        'HorizontalAlignment','center',...
        'FontWeight','bold',...
        'FontSize',32,...
        'FitBoxToText','off',...
        'EdgeColor','none');
end


% %% ----------------- SAVE OUTPUT FIGURES -----------------

print(gcf, fullfile(output_dir, 'Overall_paperfig.png'), '-dpng', '-r300');
print(gcf, fullfile(fig_dir, 'Overall_paperfig.png'), '-dpng', '-r300');


% set(gcf, 'PaperPositionMode', 'auto');
% print(gcf, fullfile(fig_dir, 'Overall_paperfig.pdf'), '-dpdf');
% 
% % Export exactly what is on screen, at 300 dpi
% exportgraphics(gcf, fullfile(output_dir, 'Overall_paperfig.png'), ...
%     'Resolution', 300);
% 
% prepare_figure_export('publication');  % small, readable in Word
% % OR
% %prepare_figure_export('presentation');  % big screen display
% 
% print(gcf, fullfile(output_dir, 'figure_export.pdf'), '-dpdf', '-r300');
% print(gcf, fullfile(fig_dir, 'figure_export.pdf'), '-dpdf', '-r300');
% %