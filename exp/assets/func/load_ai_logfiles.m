function [data,time] = load_ai_logfiles(ai_logfiles,threshold,sampleRate)
% load_ai_logfiles(ai_logfiles) loads all analog input logfiles (*.daq)
% specified with ai_logfiles (e.g., nt_data.ai_logfile) and combines them
% to a single data and time vector.
%
% If threshold is greater than zero, then only the data samples are kept,
% whose absolute value is above the threshold. This considers only the
% first channel.
%
% Input
%   ai_logfiles     - cell with filenames
%   threshold       - vector for each channel
%   sampleRate      - sample rate of data (e.g., ai.SampleRate)
%
%
% Author:           Martin Grund
% Last update:      January 5, 2016

data = [];

for i = 1:length(ai_logfiles)
    
    % Get filename (files were created on Windows)
    [~,name_tmp,ext_tmp] = fileparts(strrep(ai_logfiles{i},'\','/'));        
    
    [data_tmp,time_tmp] = daqread([name_tmp ext_tmp]);
    
    data = [data; data_tmp];
    
    if i == 1
        time = time_tmp;
    else
        time = [time; time_tmp+time(end)+1/sampleRate];
    end
    
end

ind = [];

if sum(threshold)

    for i = 1:numel(threshold)
        
        if threshold(i) > 0            
            ind_tmp = find(abs(data(:,i))>threshold(i));
            
            ind = union(ind,ind_tmp);
        end
        
    end
    
    data = data(ind,:);
    time = time(ind); 
end