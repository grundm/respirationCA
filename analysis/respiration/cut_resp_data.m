function [resp_data_c] = cut_resp_data(resp_data)
% function [resp_data_c] = cut_resp_data(resp_data) selects only channel 
% for respiration data with a specified temporal window before the first 
% trigger and after the last trigger. It also returns the number of events 
% for easy assessment of missing triggers.
%
%
% Author:           Martin Grund
% Last update:      February 8, 2021

%% Settings

trigger_label = 'R128';

% Select EEG data from 10s before first trigger to 10s after last trigger

% Maximal time from onset fixation cross to trigger: 2.5 s
% Maximal time from onset cue to confidence response in last trial: 2.75 ms
% max(0.134 + dt$resp1_t[dt$trial==150] + 0.3 + dt$resp2_t[dt$trial==150])
time_padding = [5 5];

resp_channel_label = 'Resp';

%% Split first file for ID21, because block #1 and #2 recorded together

ID21_ind = find(strcmp({resp_data.ID},'21'));

if length(resp_data(ID21_ind).eeg) < 4

    % Recording of block #1 and #2
    eeg_tmp = resp_data(ID21_ind).eeg(1);

    % Copy block #3 and #4
    resp_data(ID21_ind).eeg(4) = resp_data(ID21_ind).eeg(3);
    resp_data(ID21_ind).eeg(3) = resp_data(ID21_ind).eeg(2);

    % Get trigger indices
    trigger_ind = find(contains({eeg_tmp.event.type},trigger_label));

    % End block #1 = 1st trigger of #2 block (151) - 1s
    block01_end = eeg_tmp.event(trigger_ind(151)).latency - eeg_tmp.srate;

    % Start block #2 = last trigger #1 block (150) + 1s
    block02_start = eeg_tmp.event(trigger_ind(150)).latency + eeg_tmp.srate;

    % Split block #1 and #2
    resp_data(ID21_ind).eeg(1) = pop_select(eeg_tmp, 'point', [1 block01_end]);        
    resp_data(ID21_ind).eeg(2) = pop_select(eeg_tmp, 'point', [block02_start eeg_tmp.pnts]);
    
end
        

%% Loop participants
for i = 1:length(resp_data)
    
    resp_data_c(i).ID = resp_data(i).ID;    
    
    % Number of files
    resp_data_c(i).file_num = length(resp_data(i).eeg);
    
    for j = 1:length(resp_data(i).eeg)

        % Store sampling rate and number of triggers
        resp_data_c(i).srate(1,j) = resp_data(i).eeg(j).srate;
        %resp_data_c(i).event_num(1,j) = length(resp_data(i).eeg(j).event);
        resp_data_c(i).trigger_num(1,j) = sum(contains({resp_data(i).eeg(j).event.type},trigger_label));
        
        % Get index of first an last trigger
        trigger_first_ind = find(contains({resp_data(i).eeg(j).event.type},trigger_label),1,'first');
        trigger_last_ind = find(contains({resp_data(i).eeg(j).event.type},trigger_label),1,'last');

        % Caculate interval from before first trigger to after last trigger
        % Removes data that reflects not the experiment instructions, talking,
        % etc.
        interval_select = [resp_data(i).eeg(j).event(trigger_first_ind).latency - time_padding(1)*resp_data(i).eeg(j).srate, ...
                           resp_data(i).eeg(j).event(trigger_last_ind).latency + time_padding(2)*resp_data(i).eeg(j).srate];

        resp_channel_ind = find(contains({resp_data(i).eeg(j).chanlocs.labels},resp_channel_label));                
                       
        % Select data based on specified interval and respiration channel
        resp_data_c(i).eeg(j) = pop_select(resp_data(i).eeg(j), 'point', interval_select, 'channel', resp_channel_ind);
    
    end
    
    % Sum up number of events across files/blocks    
    resp_data_c(i).trigger_num(1,end+1) = sum(resp_data_c(i).trigger_num(1,:),2);
    
end
    
