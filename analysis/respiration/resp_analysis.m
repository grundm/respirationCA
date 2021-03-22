%% Respiratory analyis
%
% Loads respiratory traces and does respiratory cycle detection
%
%
% Author:           Martin Grund
% Last update:      February 8, 2021

%% Settings

code_dir = '/Users/martin/ownCloud/promotion/experiment/respirationCA/code/respirationca';

data_dir = '/Users/martin/ownCloud/promotion/experiment/respirationCA/data';


%% Set paths

cd(code_dir)

% Add paths
addpath([code_dir '/respiration']);
addpath([code_dir '/assets/hline_vline']);
addpath([code_dir '/assets/eeglab14_1_1b']);

% Initialize EEGLab
eeglab


%% Load and save data (vhdr-files)

resp_data = load_save_resp_data(data_dir);
% 35s


%% Select and trim respiratory channel

resp_data_c = cut_resp_data(resp_data);
% 52s

%save([data_dir '/resp_data_c.mat'], 'resp_data_c', '-v7.3');


%% Respiratory cycle detection (peak/troughs)
%load([data_dir '/resp_data_c.mat']);

[resp_cycles] = find_resp_cycle(resp_data_c);
% 220s


%% Trigger onsets relative to respiratory cycle

[stim_data] = locate_trigger_resp(resp_cycles);

% Write out data
csvwrite([data_dir '/resp_stim_circ_20200208_smoothed.csv'],stim_data);


%%
%% Visual inspection of trough detection

% Check ID40 #4 -> very long respiratory cycle

i = find(strcmp({resp_cycles.ID},'29')); % ID
j = 3; % block

figure('Name', ['ID' resp_cycles(i).ID ' #' num2str(j) '']);
%title(['ID' resp_cycles(i).ID ' #' num2str(j) ''])

resp = filloutliers(resp_data_c(i).eeg(j).data,'linear','movmedian',srate_tmp); % remove outliers
smoothresp = smoothdata(resp,'sgolay',srate_tmp); % apply some gentle smoothing
zsmoothresp = zscore(smoothresp); 

plot(zsmoothresp,'g-');
hold on
vline(resp_cycles(i).troughs(j).loc,'r-'); % Inhale onsets
%vline(resp_cycles(i).trigger(j).loc,'k'); % Trigger

%% Optimized quality check of inhale onset

%j = 4;

min_peak_dist = 2;

for i = 8:8
    
    for j = 2:2

        raw_data = resp_data_c(i).eeg(j).data;
        srate_tmp = resp_data_c(i).eeg(j).srate;
        ID_tmp = resp_data_c(i).ID;

        resp = filloutliers(raw_data,'linear','movmedian',srate_tmp);
        %resp = filloutliers(raw_data,'pchip');
        smoothresp = smoothdata(resp,'sgolay',srate_tmp);
        %smoothresp = resp;
        zsmoothresp = zscore(smoothresp); % z-score

        figure('Name',['ID' ID_tmp ' #' num2str(j) ': raw (black), w/o outlier (blue), and smoothed (red)']); 
        %plot(raw_data,'k-'); hold on; plot(resp,'b-'); 
        plot(smoothresp,'r-')

        figure('Name',['ID' ID_tmp ' #' num2str(j) ': Inhale onset detection (smoothed & zscored)']); 
        plot(-zsmoothresp,'r-');
        %findpeaks(-zsmoothresp,'minpeakdistance',min_peak_dist*srate_tmp,'minpeakprominence',0.5)
        findpeaks(-zsmoothresp,'minpeakdistance',min_peak_dist*srate_tmp,'minpeakprominence',iqr(zsmoothresp)*0.9)

        % figure('Name',['ID' ID_tmp ' #' num2str(j) ': Inhale onset detection (raw data)']); 
        % plot(-raw_data,'r-');
        % findpeaks(-raw_data,'minpeakdistance',min_peak_dist*srate_tmp,'minpeakprominence',iqr(raw_data)*0.9)

        %vline(resp_cycles(i).trigger(j).loc,'k'); % Trigger

        % figure('Name',['ID' ID_tmp ' #' num2str(j) ': Exhale onset detection (zscored)']);
        % plot(zsmoothresp,'b-');
        % findpeaks(zsmoothresp,'minpeakdistance',min_peak_dist*srate_tmp,'minpeakprominence',0.5)
    end

end

%% Data quality check of respiratory cycle duration

k = 0;
figure
%for i = 1:24       
for i = 25:41       
    
    for j = 1:length(resp_cycles(i).troughs)
        
        k = k + 1;
        subplot(8,12,k)        
        resp_cycle_t = diff(resp_cycles(i).troughs(j).loc)./resp_cycles(i).srate(j);
        boxplot(resp_cycle_t);
        %boxplot(resp_cycle_t(resp_cycle_t < (median(resp_cycle_t)*3)));
        title(['ID' resp_cycles(i).ID ' #' num2str(j) '']);
        %hline(median(resp_cycle_t)*1.5,'r-')
        hline(median(resp_cycle_t)*2,'k-')
        %disp(['ID' resp_cycles(i).ID ' #' num2str(j) ' - number of outlier: ' num2str(sum(isoutlier(resp_cycle_t)))]);
        disp(['ID' resp_cycles(i).ID ' #' num2str(j) ' - number of outlier: ' num2str(sum(resp_cycle_t > (median(resp_cycle_t)*2)))]);
        disp(['ID' resp_cycles(i).ID ' #' num2str(j) ' - median: ' num2str(median(resp_cycle_t))]);
    end
end

