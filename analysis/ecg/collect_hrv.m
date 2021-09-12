
% 
%
%
% Author:           Martin Grund
% Last update:      June 18, 2021

%% Settings

data_dir = '/data/pt_02099/kubios/';
data_dir = '/Users/martin/ownCloud/promotion/experiment/respirationCA/data/kubios/';

kubios_file_wildcard = 'Kubios_ID*_hrv.mat';

output_file = 'respirationCA_RMSSD.csv';

%% Get RMSSD data for each participant and block

% Get files
kubios_files = dir([data_dir kubios_file_wildcard]);

hrv_data = zeros(length(kubios_files),4);

% Loop files
for i = 1:length(kubios_files)
   
    kubios_data_tmp = load([kubios_files(i).folder '/' kubios_files(i).name]);
    
    % Get ID and block
    filename_split_tmp = split(erase(kubios_data_tmp.Res.f_name, 'Kubios_ID'),'_');
    
    hrv_data(i,1) = str2double(filename_split_tmp{1});
    block_tmp = split(filename_split_tmp{2},'.');
    hrv_data(i,2) = str2double(block_tmp{1});

    hrv_data(i,3) = kubios_data_tmp.Res.HRV.Statistics.RMSSD;
    hrv_data(i,4) = kubios_data_tmp.Res.HRV.Statistics.mean_HRV;
end

%% Save data

writematrix(hrv_data,[data_dir '/' output_file])

