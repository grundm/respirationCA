function oxi_data = load_save_data(data_dir)
% function oxi_data = load_save_data(data_dir) loads all /ID*/*.acq files
% into one structure oxi_data and saves the structure in the data directory
% (data_dir).
%
% Author:           Martin Grund
% Last update:      August 25, 2020

%% Load and save data

% Get ID directories
ID_dir = dir([data_dir '/ID*']);

for i = 1:length(ID_dir)
    
    % Get ACQ-files within each ID directory
    acq_files = dir([ID_dir(i).folder '/' ID_dir(i).name '/*.acq']);
    
    % Add ID string
    oxi_data(i).ID = erase(ID_dir(i).name,'ID');
    
    for j = 1:length(acq_files)
        
        % Load ACQ-file into structure
        oxi_data(i).acq(j) = load_acq([acq_files(j).folder '/' acq_files(j).name]);
    end
end

% Save OXI data
save([data_dir '/oxi_data_all.mat'], 'oxi_data', '-v7.3');