%% Get somatosensory evoked potentials (SEPs)
% 
% Author:           Martin Grund, Tilman Stephani
% Last Update:      March 24, 2020

function res = sep_analysis(ID,data_dir)

res.ID = ID;

%% Data files

% Data directory
%data_path = ['/data/pt_02099/ID' ID '/'];
data_path = [data_dir '/ID' ID '/'];

% ECG & respiration recording
wildcard_resp_data = 'respirationCA_*.vhdr';

% Experimental data (mat-files)
if strcmp(ID,'05')

    % txt-file for ID05 because mat-files were overwritten
    trial_data_txt_str = [data_path 'respirationCA_' ID '_trials_0'];

else
    session_file_str = [data_path 'respirationCA_' ID '_data_0'];
end

%% Missed triggers due to late recording start

% Missed first 4 triggers due to late recording start of block #2
if strcmp(ID,'03')
    missed_trials = 4;
    missed_trials_block = 2;
    
% Missed first 2 triggers due to late recording start of block #1    
elseif strcmp(ID,'22')    
    missed_trials = 22;
    missed_trials_block = 1;
    
end

%% SEP data handling

trigger_name = 'R128';

% Select EEG data from 10s before first trigger to 10s after last trigger
time_padding = [10 10];

if str2double(ID)<18 || strcmp(ID,'1001')
   
    line = 2; % SEP
    %line = 3; % Respiration
    
else    
    line = 0;
    %line = 2; % Respiration
    
end

%% Trial filter

% Valid button in maximum time
buttons_valid = [1 2];
resp_t_max = 1.5;


%% Load and cut Brain Products data

% Get all available vhdr-files
hdrfiles = dir([data_path wildcard_resp_data]);

for i = 1:length(hdrfiles)
    
    % Load data with EEGLAB
    EEG_raw(i) = pop_loadbv(data_path, hdrfiles(i).name);
    
    % Get index of first an last trigger
    trigger_first = find(contains({EEG_raw(i).event.type},trigger_name),1,'first');
    trigger_last = find(contains({EEG_raw(i).event.type},trigger_name),1,'last');
    
    % Caculate interval from before first trigger to after last trigger
    % Removes data that reflects not the experiment instructions, talking,
    % etc.
    interval_select = [EEG_raw(i).event(trigger_first).latency - time_padding(1)*EEG_raw(i).srate, ...
                       EEG_raw(i).event(trigger_last).latency + time_padding(1)*EEG_raw(i).srate];
    
    % Select data based on specified interval and SEP channel
    EEG_select(i) = pop_select(EEG_raw(i), 'point', interval_select, 'channel', line);
    
    % Merge blocks continuously 
    if i == 2
        EEG = pop_mergeset(EEG_select(1),EEG_select(2));
    elseif i > 2
        EEG = pop_mergeset(EEG,EEG_select(i));
    end    
    
end

res.EEG_srate = EEG.srate

%% Optional: Interpolate stimulation artifact

% Find event latencies
    % where to cut (triggers) 
    ix_A_Out = find(ismember({EEG.event.type}, trigger_name)); % stimulation events A
    lat_A_Out = [EEG.event(ix_A_Out).latency];

    lat_all_Out = lat_A_Out;
    
        
% pchip = "cubic monotonous hermite spline interpolation" (Gunnar Waterstraat)
        
    % cut out samples (from -2 to 4 ms)
    sr = EEG.srate;
    t_cut = [-2 4]; % interpolated interval around triggers; in ms 
    pt_cut = ceil(t_cut/1000*sr); % in sampling points

    % Piecewise Cubic Hermite Interpolating Polynomial (PCHIP) + replace EEG data
    n_pt_fit = 5; %+1;  number of samples before and after cut used for interpolation fit
    
    x_fit_raw = [pt_cut(1)-n_pt_fit : 1 : pt_cut(1), pt_cut(2) : 1 : pt_cut(2)+n_pt_fit];
    x_sr_raw = [pt_cut(1) : 1 : pt_cut(2)]; % points to be interpolated; in pt

    for i = 1:length(lat_all_Out) % loop through all stimulation events
        x_fit = lat_all_Out(i) + x_fit_raw; % fit point latencies for this event
        x_sr = lat_all_Out(i) + x_sr_raw; % latencies for to-be-interpolated data points

        for c = 1:size(EEG.data,1) % loop through all channels
            y_fit = EEG.data(c, x_fit); % y values to be fitted
            %y_interp = pchip(x_fit, y_fit, EEG.data(c, x_sr-x_sr_raw)); % calculate pchip, obtain values in t_cut interval
            y_interp = pchip(x_fit, y_fit, x_sr); % calculate pchip, obtain values in t_cut interval
            EEG.data(c, x_sr) = y_interp; % replace in EEG data
        end

        if mod(i, 100) == 0 % show message every 100th trial
            fprintf('stimulation event %d \n', i)
        end
    end

%% Inspect power spectrum of data
% figure; pop_spectopo(EEG, 1, [], 'EEG' , 'percent', 100, 'freqrange',[1 100],'electrodes','off'); % spectrogram

%% Filter LOW
%[b,a] = butter(2, [low_pass_filter_resp]/(EEG.srate/2), 'low'); % low-pass filter
% Bach et al. (2016): 5 Hz low-pass
%EEG.data = filtfilt(b, a, double(EEG.data)')'; % apply filter

%% Filter
[b,a] = butter(2, [70]/(EEG.srate/2), 'high'); % high-pass filter
EEG.data = filtfilt(b, a, double(EEG.data)')'; % apply filter

% notch filter 1
% bsFilter = [48 52];
% [b_notch, a_notch] = butter(2, bsFilter/(EEG.srate/2),'stop');
% EEG.data = filtfilt(b_notch, a_notch, double(EEG.data)')';

%% Epoch
trigger_names = {'R128'};
epoch_length_t = [-0.100 0.100];
%epoch_length_t = [-1.000 5.000]; % -5 does not work with ID03
[EEG, ~] = pop_epoch(EEG, trigger_names, epoch_length_t, 'newname', ['Data_epoched'], 'epochinfo', 'yes');
EEG = pop_rmbase(EEG, [-50 -2]); % not really needed when high high-pass filter is used


%% GET EVENT INDICES
% Trial number of conditions to match with triggers

for k = 1:4
    
    if strcmp(ID,'05') % overwrote mat files so we use txt files
        
        block.exp_data = readtable([trial_data_txt_str num2str(k) '.txt']);
        trial_num(k) = size(block.exp_data,1);
    else
        
        block = load([session_file_str num2str(k) '.mat']);
        trial_num(k) = size(block.exp_data.seq,1);
    end
        
    
    % Valid trials (valid button press within maxim response time)
    ind.valid_trials{k} = block.exp_data.resp1_t < resp_t_max ...
                          & block.exp_data.resp2_t < resp_t_max ...
                          & any(block.exp_data.resp1_btn == buttons_valid,2) ...
                          & any(block.exp_data.resp2_btn == buttons_valid,2);

    % Ignore trials when triggers were not recorded due to late start
    if strcmp(ID,'03') || strcmp(ID,'22')

        if k == missed_trials_block
            
            ind.valid_trials{missed_trials_block}(1:missed_trials) = 0;
        end
    end
        
    % Near
    ind.near{k} = find(block.exp_data.intensity & ind.valid_trials{k});
    
    % Near - Hit        
    ind.near_hit{k} = find(block.exp_data.intensity > 0 & block.exp_data.resp1 == 1 & ind.valid_trials{k});
    
    % Near - Miss
    ind.near_miss{k} = find(block.exp_data.intensity > 0 & block.exp_data.resp1 == 0 & ind.valid_trials{k});
    
    % Null
    ind.null{k} = find(block.exp_data.intensity == 0 & ind.valid_trials{k});
end

% Missed first 4 triggers due to late recording start

% Adapt total number of trials to missed triggers

if strcmp(ID,'03') || strcmp(ID,'22')
    
    trial_num(missed_trials_block) = trial_num(missed_trials_block)-missed_trials;

end


% Accumulutating number of trials
trial_num_sum = [trial_num(1) sum(trial_num(1:2)) sum(trial_num(1:3)) sum(trial_num(1:4))];

% Correct indices of following block by trial number of previous block
res.near_i = [ind.near{1}; ind.near{2}+trial_num_sum(1); ind.near{3}+trial_num_sum(2); ind.near{4}+trial_num_sum(3)];
res.near_hit_i = [ind.near_hit{1}; ind.near_hit{1}+trial_num_sum(2); ind.near_hit{2}+trial_num_sum(3); ind.near_hit{4}+trial_num_sum(3)];
res.near_miss_i = [ind.near_miss{1}; ind.near_miss{1}+trial_num_sum(2); ind.near_miss{2}+trial_num_sum(3); ind.near_miss{4}+trial_num_sum(3)];
res.null_i = [ind.null{1}; ind.null{2}+trial_num_sum(1); ind.null{3}+trial_num_sum(2); ind.null{4}+trial_num_sum(3)];

res.trial_num_sum = trial_num_sum;

%% ERP
elec = 1; % channel index
res.ERP_near = mean(EEG.data(elec, :, res.near_i),3);
res.ERP_near_hit = mean(EEG.data(elec, :, res.near_hit_i),3);
res.ERP_near_miss = mean(EEG.data(elec, :, res.near_miss_i),3);
res.ERP_null = mean(EEG.data(elec, :, res.null_i),3);

res.EEG_times = EEG.times;

%% Plot
% elec = 1; % channel index
% figure; hold on;
% plot(EEG.times, mean(EEG.data(elec, :, res.near_i),3), 'r')
% plot(EEG.times, mean(EEG.data(elec, :, res.null_i),3), 'k')
% xlim([-20 60])
% title('Average SSEP NEAR')


end