%% Analysis of oximeter data
%
% Loads all Biopac ACQ files. Cleans them (hard coded), does peak detection
% and locates stimulus onsets relative in the pulse wave cycle.
%
% Furthermore, pulse waves are extracted and resampled to cardiac cycle.
%
% Author:           Martin Grund
% Last update:      November 12, 2020

%% Settings

code_dir = '/Users/martin/ownCloud/promotion/experiment/respirationCA/code/respirationca';

data_dir = '/Users/martin/ownCloud/promotion/experiment/respirationCA/data';


%% Set paths

cd(code_dir)

% Add paths
addpath([code_dir '/oxi']);
%addpath(genpath([code_dir '/assets']));


%% Load and save data initially (ACQ-files)

oxi_data = load_save_data(data_dir);

%% Load data (Structure of load_save_data()

load([data_dir '/oxi_data_all.mat']);

%% Clean data

[oxi_data_c,res_c,acq_num_c] = clean_data(oxi_data);

% Save data
save([data_dir '/oxi_data_all_c.mat'], 'oxi_data_c', 'res_c', 'acq_num_c', '-v7.3');

%% Load clean data (Structure of load_save_data()

load([data_dir '/oxi_data_all_c.mat']);


%% Oxi peak detection

data_rate = 1000; % data sampling rate

srate = 50; % sensor sampling rate

max_bpm = 140; % to calculate min_peak_dist ('MinPeakDistance')
min_peak_dist = (srate*60)/max_bpm; % samples per minute (60 seconds) devided by beats per minute
devide_range = 10; % Fraction of signal range for 'MeanPeakProminence'

% Loop participants
for i = 1:length(oxi_data_c)
        
    % Get index in results structure for current ID
    ID_ind_res = find(strcmp({res_c.ID},oxi_data_c(i).ID));
    
    % Loop blocks
    for j = 1:length(oxi_data_c(i).acq)
        
%         max_bpm = max(oxi_data_c(i).acq(j).data(:,3))*1.1 % to calculate min_peak_dist ('MinPeakDistance')
%         min_peak_dist = (srate*60)/max_bpm; % samples per minute (60 seconds) devided by beats per minute

        % "Resample" data to actual sensor sampling rate
        data_res_tmp = decimate(oxi_data_c(i).acq(j).data(:,2), data_rate/srate);
        
        % "Resample" trigger time to actual sensor sampling rate
        res_c(ID_ind_res).trigger(j).time_res = res_c(ID_ind_res).trigger(j).time / (data_rate/srate);
        
        % Get oxi peaks        
        [res_c(ID_ind_res).pulse(j).pks,...
         res_c(ID_ind_res).pulse(j).locs,...
         res_c(ID_ind_res).pulse(j).w,...
         res_c(ID_ind_res).pulse(j).p] = findpeaks(data_res_tmp,...
                                          'MinPeakDistance',min_peak_dist,...
                                          'MinPeakProminence',range(data_res_tmp)/devide_range);
                                    
    end
end

%% Save oxi peak data

save([data_dir '/oxi_data_peaks.mat'], 'res_c', '-v7.3');

%% Load oxi peak

load([data_dir '/oxi_data_peaks.mat']);


%% Get oxi interval for each trigger

srate = 50; % sensor sampling rate

l = 0;

% Loop participants
for i = 1:length(res_c)
        
    % Loop blocks
    for j = 1:length(res_c(i).trigger)
        
        % Started recording too late for ID22 in block 1
        if strcmp(res_c(i).ID,'22') && j == 1
            
            for k = 1:3
                l = l+1;
                stim_data(l,1:6) = [str2double(res_c(i).ID) j k NaN NaN NaN];
            end
        end
        
        % Loop trials
        for k = 1:length(res_c(i).trigger(j).time)
            
          % Get index of oxi peak before stimulus onset (trigger)                    
          prev_oxi_peak_ind = find(res_c(i).pulse(j).locs < res_c(i).trigger(j).time_res(k),1,'last');

          % If recording started too late for first trigger
          if isempty(prev_oxi_peak_ind)
              i
              j
              k
              prev_oxi_peak_ind
              diff2peak = NaN;
              pulse_t = NaN;
              pulse_t_next = NaN;
              stim_degree = NaN;
          
          else
              
              % Calculate time distance between stimulus onset the previous oxi peak
              diff2peak = (res_c(i).trigger(j).time_res(k) - res_c(i).pulse(j).locs(prev_oxi_peak_ind))/srate;        

              % Calculate interval to subsequent oxi peak in samples
              pulse_t = (res_c(i).pulse(j).locs(prev_oxi_peak_ind+1) - res_c(i).pulse(j).locs(prev_oxi_peak_ind))/srate;
              
              % Calculate interval of following pulse cycle in samples
              pulse_t_next = (res_c(i).pulse(j).locs(prev_oxi_peak_ind+2) - res_c(i).pulse(j).locs(prev_oxi_peak_ind+1))/srate;

              % Calculate relative position of the stimulus onset to the oxi cycle
              stim_degree = (diff2peak/pulse_t) * 360;
          end

          % Store in matrix with 1 row per trial
          l = l+1;
          if strcmp(res_c(i).ID,'22') && j == 1
            
            stim_data(l,1:7) = [str2double(res_c(i).ID) j k+3 diff2peak pulse_t stim_degree pulse_t_next];
              
          else
            stim_data(l,1:7) = [str2double(res_c(i).ID) j k diff2peak pulse_t stim_degree pulse_t_next];          
          end
          
        end % loop triggers/trials        
        
        % Add filter for outlier pulse_t
        pulse_t_median = median(stim_data(stim_data(:,1)==str2double(res_c(i).ID) & stim_data(:,2)==j,5),'omitnan');
        
        stim_data(stim_data(:,1)==str2double(res_c(i).ID) & stim_data(:,2)==j,8) = stim_data(stim_data(:,1)==str2double(res_c(i).ID) & stim_data(:,2)==j,5) < pulse_t_median*1.5;
        
    end
end


%% Write out data for trigger in pulse wave as table

csvwrite([data_dir '/oxi_stim_circ_res.csv'],stim_data);


%% Load data of OXI and ECG matched for only valid trials

dt_oxi_valid = readtable([data_dir '/dt_oxi_valid.csv']);

% Works better with if NaN
% uiopen([data_dir '/dt_oxi_valid.csv'],1)

%% Scale pulse wave to each trial's cardiac cycle

% Settings
sample_t = 0.02; % in s (50 Hz -> 20-ms samples)

data_rate = 1000; % data sampling rate

srate = 50; % sensor sampling rate

% "Degree"-steps for resample
resample_steps = 0:4:356;

% Prepare output matrix
pulse_cc = zeros(height(dt_oxi_valid),3 + length(resample_steps));

% Loop trials
for i = 1:length(dt_oxi_valid.ID)
    
    % Get index in data structure for current ID
    if dt_oxi_valid.ID(i) < 10
        ID_ind_res = find(strcmp({res_c.ID},['0' num2str(dt_oxi_valid.ID(i))]));
        ID_ind_data = find(strcmp({oxi_data_c.ID},['0' num2str(dt_oxi_valid.ID(i))]));
    else
        ID_ind_res = find(strcmp({res_c.ID},num2str(dt_oxi_valid.ID(i))));
        ID_ind_data = find(strcmp({oxi_data_c.ID},num2str(dt_oxi_valid.ID(i))));
    end
    
    % Get block and trial number
    j = dt_oxi_valid.block(i);
    k = dt_oxi_valid.trial(i);
    
    if dt_oxi_valid.ID(i) == 22 && j == 1
        prev_oxi_peak_ind = find(res_c(ID_ind_res).pulse(j).locs < res_c(ID_ind_res).trigger(j).time_res(k-3),1,'last');

    else
          % Get index of oxi peak before stimulus onset (trigger)                    
          prev_oxi_peak_ind = find(res_c(ID_ind_res).pulse(j).locs < res_c(ID_ind_res).trigger(j).time_res(k),1,'last');
    end

    % Get sample index of ECG R peak previous to pulse wave peak
    prev_R_peak_ind = fix(res_c(ID_ind_res).pulse(j).locs(prev_oxi_peak_ind) - (dt_oxi_valid.R2pulse(i)/sample_t));

    % Get sample index of ECG R peak subsequent to pulse wave peak
    next_R_peak_ind = ceil(prev_R_peak_ind + (dt_oxi_valid.ecg_cycle_t(i)/sample_t));

    % "Resample" data to actual sensor sampling rate
    data_res_tmp = decimate(oxi_data_c(ID_ind_data).acq(j).data(:,2), data_rate/srate);

    % Get pulse wave data for length of cardiac cycle
    pulse_in_cc = data_res_tmp(prev_R_peak_ind:next_R_peak_ind);

    %% Resample pulse wave data to fixed number of samples mirroring cardiac
    % cycle in degrees
    
    % Padding based on pop_resample (myresample(),  Andreas Widmann May 5,
    % 2011)
    pnts = length(resample_steps);
    new_pnts = length(pulse_in_cc);
    
    [p, q] = rat(pnts / new_pnts, 1e-12); % Same precision as in resample
    N = 10; % Resample default
    nPad = ceil((max(p, q) * N) / q) * q; % # datapoints to pad, round to integer multiple of q for unpadding
    pulse_in_cc_deg = resample([pulse_in_cc(ones(1, nPad), :); pulse_in_cc; pulse_in_cc(end * ones(1, nPad), :)], pnts, new_pnts);
    nPad = nPad * p / q; % # datapoints to unpad
    pulse_in_cc_deg = pulse_in_cc_deg(nPad + 1:end - nPad, :); % Remove padded data
    

    %% Combine in matrix with ID, block and trial number
    % Resampled pulse wave
    pulse_cc(i,:) = [dt_oxi_valid.ID(i) j k pulse_in_cc_deg'];
    
    % Original pulse wave since R-peak
    pulse_since_R(i,:) = [dt_oxi_valid.ID(i) j k data_res_tmp(prev_R_peak_ind:prev_R_peak_ind+99)'];

end

%% Save data from resampling pulse waves

save([data_dir '/pulse_cc.mat'], 'pulse_cc');

save([data_dir '/pulse_since_R.mat'], 'pulse_since_R');

%% Load data from resampling pulse waves

%save([data_dir '/pulse_cc.mat'], 'pulse_cc');

load([data_dir '/pulse_since_R.mat']);


%% Avergade pulse wave since R-peak for each participant

ID_list = unique(pulse_since_R(:,1));

mean_pulse_since_R = zeros(length(ID_list), size(pulse_since_R,2)-2);
mean_pulse_since_R_miss = zeros(length(ID_list), size(pulse_since_R,2)-2);
mean_pulse_since_R_hit = zeros(length(ID_list), size(pulse_since_R,2)-2);

for i = 1:length(ID_list)    
    
    mean_pulse_since_R(i,:) = [ID_list(i) mean(pulse_since_R(pulse_since_R(:,1)==ID_list(i),4:end))];
    
    mean_pulse_since_R_miss(i,:) = [ID_list(i) mean(pulse_since_R(pulse_since_R(:,1)==ID_list(i) & dt_oxi_valid.stim_type == 1  & dt_oxi_valid.resp1 == 0,4:end))];
    mean_pulse_since_R_hit(i,:) = [ID_list(i) mean(pulse_since_R(pulse_since_R(:,1)==ID_list(i) & dt_oxi_valid.stim_type == 1  & dt_oxi_valid.resp1 == 1,4:end))];

end

% Normalize R-peak locked pulse wave by substracting the minimum and
% dividing by data span
mean_pulse_since_R_norm = (mean_pulse_since_R(:,2:end) - min(mean_pulse_since_R(:,2:end),[],2)) ./ range(mean_pulse_since_R(:,2:end),2);


%% Avergade resampled pulse wave in cardiac cycle for each participant

ID_list = unique(pulse_cc(:,1));

mean_pulse_cc = zeros(length(ID_list), size(pulse_cc,2)-2);

for i = 1:length(ID_list)
    
    mean_pulse_cc(i,:) = [ID_list(i) mean(pulse_cc(pulse_cc(:,1)==ID_list(i),4:end))];        
    
end

% Normalize resample pulse wave
mean_pulse_cc_norm = (mean_pulse_cc(:,2:end) - min(mean_pulse_cc(:,2:end),[],2)) ./ range(mean_pulse_cc(:,2:end),2);


%% Plot mean pulse waves since R-peak or within cardiac cycle

% Add paths
addpath([code_dir '/assets/hline_vline']);

figure

subplot(2,1,1)

% Plot (normalized) mean pulse for each participant
for i = 1:length(mean_pulse_since_R(:,1))

    line_plot1(i) = plot(0:40,mean_pulse_since_R_norm(i,1:41),'LineWidth',1);
    hold on
    
    line_plot1(i).Color(4) = 0.4;

end

% Plot mean across participants
plot(0:40,mean(mean_pulse_since_R_norm(:,1:41)),'LineWidth',4)


xticks(0:2.5:40)
xticklabels((0:2.5:40)*20)

lines = vline([150 300]/20,'k-');
set(lines,'LineWidth',2)

vline(250/20,'k-');

title('Pulse wave since R-peak')
ylabel('Normalized mean pulse wave');
xlabel('Time since R-peak in ms');

subplot(2,1,2)

% Plot (normalized) mean pulse for each participant
for i = 1:length(mean_pulse_since_R(:,1))

    line_plot2(i) = plot(1:40,diff(mean_pulse_since_R_norm(i,1:41))/20,'LineWidth',1);
    hold on
    
    line_plot2(i).Color(4) = 0.4;

end

% Plot mean across participants
line_mean2 = plot(1:40,mean(diff(mean_pulse_since_R_norm(:,1:41),1,2)./20),'Color',[0 0.4470 0.7410],'LineWidth',4);

xticks(0:2.5:40)
xticklabels((0:2.5:40)*20)

lines = vline([150 300]/20,'k-');
set(lines,'LineWidth',2)

vline(250/20,'k-');

title('First derivative of pulse wave since R-peak')
ylabel('d(pulse)/d(time)');
xlabel('Time since R-peak in ms');

hline(0,'k-');

