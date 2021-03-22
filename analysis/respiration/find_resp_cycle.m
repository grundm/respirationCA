function [resp_cycles] = find_resp_cycle(resp_data_c)
% function [resp_cycles] = find_resp_cycle(resp_data_c) does exhale (peak)
% and inhale (trough) detection in respiratory signal after outlier removal
% and smoothing.
%
%
% Author:           Martin Grund
% Last update:      February 4, 2021

%% Settings

trigger_label = 'R128';

% Minimal temporal distance between two peaks/troughs in seconds
min_peak_dist = 2; % in s

% Fraction of z-smoothed data interquartile range for minimum peak
% prominence
frac_iqr = 0.9;
% min_peak_prom = 0.5; -> too low for our data
% min_peak_prom = 1;
% 
% devide_range = 10; % Fraction of signal range for 'MinPeakProminence'

%% Loop participants
for i = 1:length(resp_data_c)
    
    resp_cycles(i).ID = resp_data_c(i).ID;
    
    % Loop files/blocks
    for j = 1:length(resp_data_c(i).eeg)
        
        % Add trigger latencies
        resp_cycles(i).trigger(j).loc = cell2mat({resp_data_c(i).eeg(j).event(find(contains({resp_data_c(i).eeg(j).event.type},trigger_label))).latency});
        
        % Store sampling rate temporarily
        srate_tmp = resp_data_c(i).eeg(j).srate;
        resp_cycles(i).srate(1,j) = srate_tmp;
        
        % Preprocess the respiratory trace (Power et al., 2020)
        resp = filloutliers(resp_data_c(i).eeg(j).data,'linear','movmedian',srate_tmp); % remove outliers
        smoothresp = smoothdata(resp,'sgolay',srate_tmp); % apply some gentle smoothing
        zsmoothresp = zscore(smoothresp); % z-score

        %min_peak_prom = range(zsmoothresp)/devide_range;
        min_peak_prom = iqr(zsmoothresp)*frac_iqr;

        % Find peaks (exhale onsets) and troughs (inhale onsets)
        [~,resp_cycles(i).peaks(j).loc,~,~] = findpeaks(zsmoothresp,'minpeakdistance',min_peak_dist*srate_tmp,'minpeakprominence',min_peak_prom);
        [~,resp_cycles(i).troughs(j).loc,~,~] = findpeaks(-zsmoothresp,'minpeakdistance',min_peak_dist*srate_tmp,'minpeakprominence',min_peak_prom);
        
    end         
    
end
