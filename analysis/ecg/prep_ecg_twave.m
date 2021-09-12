%% Prepare ECG data for t-wave detection
%
% Loads ECG traces, cuts, filters and saves them for further processing
%
%
% Author:           Martin Grund
% Last update:      June 27, 2021

%% Settings

code_dir = '/Users/martin/ownCloud/promotion/experiment/respirationCA/code/respirationca';

data_dir = '/Users/martin/ownCloud/promotion/experiment/respirationCA/data';

output_dir = [data_dir '/kubios/ecg_trace'];

%% Set paths

cd(code_dir)

% Add paths
addpath([code_dir '/respiration']);
addpath([code_dir '/assets/hline_vline']);
addpath([code_dir '/assets/eeglab14_1_1b']);

% Initialize EEGLab
eeglab


%% Load and save data (vhdr-files)

tic
ecg_data = load_save_resp_data(data_dir);
toc
% 35s


%% Select and trim respiratory channel

% ecg_data_c = cut_resp_data(ecg_data,'HR');
% 52s

%% Boil down, filter and save data

mkdir(output_dir)

% Loop participants
for i = 1:length(ecg_data)

    % Loop blocks
    for j = 1:length(ecg_data(i).eeg)
        
        out_tmp.srate = ecg_data(i).eeg(j).srate;
        
        ecg_channel_ind = find(contains({ecg_data(i).eeg(j).chanlocs.labels},'HR'));
        
        out_tmp.ecg = ecg_data(i).eeg(j).data(ecg_channel_ind,:);        
        
        % Low-pass filter 30 Hz, and high-pass filter 0.5 Hz with a 4th order of Butterworth filter
        [b,a] = butter(4, [0.5]/(out_tmp.srate/2), 'high'); % high-pass filter
        out_tmp.ecg = filtfilt(b, a, double(out_tmp.ecg)')'; % apply filter

        [b,a] = butter(4, [30]/(out_tmp.srate/2), 'low'); % low-pass filter
        out_tmp.ecg = filtfilt(b, a, double(out_tmp.ecg)')'; % apply filter
%          
%         [b,a] = butter(4, [0.5 30]/(out_tmp.srate/2), 'bandpass'); % band-pass filter
%         out_tmp.ecg = filtfilt(b, a, double(out_tmp.ecg)')'; % apply filter
        
        if strcmp(ecg_data(i).ID,'21') && j > 1
            save([output_dir '/ecg_trace_' ecg_data(i).ID '_0' num2str(j+1)],'out_tmp')
        else
            save([output_dir '/ecg_trace_' ecg_data(i).ID '_0' num2str(j)],'out_tmp')
        end
    end
    
end

