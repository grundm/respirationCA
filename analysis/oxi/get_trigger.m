function [res,acq_num] = get_trigger(oxi_data)
% function [res,acq_num] = get_trigger(oxi_data) loops all participants and
% blocks to check number of blocks and triggers per blocks (acq_num) and to
% further process trigger timestamps (res(ID_ind).trigger(block).time).
%
% Author:           Martin Grund
% Last update:      August 25, 2020

for i = 1:length(oxi_data)
    
    % Assign ID string as double
    acq_num(i,1) = str2double(oxi_data(i).ID);
    
    % Check number of files per participant
    acq_num(i,2) = length(oxi_data(i).acq);
    
    % Assign ID string
    res(i).ID = oxi_data(i).ID;
    
    % Get trigger per block/file
    for j = 1:acq_num(i,2)                
        
        % Get trigger timestamps       
        res(i).trigger(j).time = find(oxi_data(i).acq(j).data(:,5));
        
        % Sum up triggers per block/file
        acq_num(i,2+j) = length(res(i).trigger(j).time);
    end
end