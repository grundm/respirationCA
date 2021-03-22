function [stim_data] = locate_trigger_resp(resp_cycles)
% function [stim_data] = locate_trigger_resp(resp_cycles) locates
% trigger/stimulus onset in respiratory cycle (starting with inhale onset)
% and return absolute and relative position for each trial/trigger.
%
%
% Author:           Martin Grund
% Last update:      February 8, 2021

%% Settings

factor_median_outlier = 2;

%% Loop participants

l = 0; % continious row counter

for i = 1:length(resp_cycles)
    
    % Loop blocks/files
    for j = 1:length(resp_cycles(i).trigger)        
        
        srate = resp_cycles(i).srate(1,j);
        
        % Started recording to late for ID03 in block #2 (missed 4 trials)
        if strcmp(resp_cycles(i).ID,'03') && j == 2
            
            for k = 1:4
                l = l+1;
                stim_data(l,1:6) = [str2double(resp_cycles(i).ID) j k NaN NaN NaN];
            end
        end
        
        % Started recording too late for ID22 in block 1 (missed 22 trials)        
        if strcmp(resp_cycles(i).ID,'22') && j == 1
            
            for k = 1:22
                l = l+1;
                stim_data(l,1:6) = [str2double(resp_cycles(i).ID) j k NaN NaN NaN];
            end
        end
        
        % Loop trials
        for k = 1:length(resp_cycles(i).trigger(j).loc)
            
          % Get index of inhale onsetbefore stimulus onset (trigger)                    
          prev_inhale_onset_ind = find(resp_cycles(i).troughs(j).loc < resp_cycles(i).trigger(j).loc(k),1,'last');

          % If recording started too late for first trigger
          if isempty(prev_inhale_onset_ind)
              i
              j
              k
              prev_inhale_onset_ind
              diff2inhale_onset = NaN;
              resp_cycle_t = NaN;
              resp_cycle_t_next = NaN;
              stim_degree = NaN;
          
          else
              
              % Calculate temporal distance from stimulus onset to
              % preceding inspiration onset
              diff2inhale_onset = (resp_cycles(i).trigger(j).loc(k) - resp_cycles(i).troughs(j).loc(prev_inhale_onset_ind))/srate;        

              % Calculate interval to subsequent trough (inhale onset)
              if (prev_inhale_onset_ind+1) <= length(resp_cycles(i).troughs(j).loc)
                resp_cycle_t = (resp_cycles(i).troughs(j).loc(prev_inhale_onset_ind+1) - resp_cycles(i).troughs(j).loc(prev_inhale_onset_ind))/srate;
              else
                i
                j
                k
                resp_cycle_t = NaN;
              end
              
              % Calculate interval of following respiratory cycle
              if (prev_inhale_onset_ind+2) <= length(resp_cycles(i).troughs(j).loc)
                resp_cycle_t_next = (resp_cycles(i).troughs(j).loc(prev_inhale_onset_ind+2) - resp_cycles(i).troughs(j).loc(prev_inhale_onset_ind+1))/srate;
              else
                i
                j
                k
                resp_cycle_t_next = NaN;
              end

              % Calculate relative position of the stimulus onset to the oxi cycle
              stim_degree = (diff2inhale_onset/resp_cycle_t) * 360;
          end

          % Store in matrix with 1 row per trial
          l = l+1;
          if strcmp(resp_cycles(i).ID,'03') && j == 2 % missed first 4 triggers in #2 block
            
            stim_data(l,1:7) = [str2double(resp_cycles(i).ID) j k+4 diff2inhale_onset resp_cycle_t stim_degree resp_cycle_t_next];         
            
          elseif strcmp(resp_cycles(i).ID,'22') && j == 1 % missed first 22 triggers in #1 block
            
            stim_data(l,1:7) = [str2double(resp_cycles(i).ID) j k+22 diff2inhale_onset resp_cycle_t stim_degree resp_cycle_t_next];         
              
          else
              
            stim_data(l,1:7) = [str2double(resp_cycles(i).ID) j k diff2inhale_onset resp_cycle_t stim_degree resp_cycle_t_next];          
          end
          
        end % loop triggers/trials        
        
        
    end % loop files/blocks            
    
    % Add median to filter outlier pulse_t
    % participant-wise
    resp_cyle_t_median = median(stim_data(stim_data(:,1)==str2double(resp_cycles(i).ID),5),'omitnan');

    stim_data(stim_data(:,1)==str2double(resp_cycles(i).ID),8) = resp_cyle_t_median;
    
end % loop participants