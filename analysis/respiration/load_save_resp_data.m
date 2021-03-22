function resp_data = load_save_resp_data(data_dir)
% function resp_data = load_save_resp_data(data_dir) loads all /ID*/*.acq files
% into one structure resp_data and saves the structure in the data directory
% (data_dir).
%
% Requires EEGLab to be initialized.
%
% Author:           Martin Grund
% Last update:      February 2, 2021

%% Load and save data

% Get ID directories
ID_dir = dir([data_dir '/ID*']);

for i = 1:length(ID_dir)
    
    % Get vhdr-files within each ID directory
    vhdr_files = dir([ID_dir(i).folder '/' ID_dir(i).name '/*.vhdr']);
    
    % Add ID string
    resp_data(i).ID = erase(ID_dir(i).name,'ID');
    
    for j = 1:length(vhdr_files)
        
        % Load EEG data into structure
        resp_data(i).eeg(j) = pop_loadbv(vhdr_files(j).folder, vhdr_files(j).name);
    end
end

% Save respiratory data
%save([data_dir '/oxi_data_all.mat'], 'oxi_data', '-v7.3');