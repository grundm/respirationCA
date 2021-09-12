function [stim_data] = locate_trigger_resp(resp_cycles)
% function [stim_data] = locate_trigger_resp(resp_cycles) locates
% trigger/stimulus onset in respiration cycle (starting with inhale and 
% exhale onset) and returns absolute/relative position for each trial/trigger.
%
%
% Author:           Martin Grund
% Last update:      June 22, 2021


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
                stim_data(l,1:11) = [str2double(resp_cycles(i).ID) j k NaN NaN NaN NaN NaN NaN NaN NaN];
            end
        end
        
        % Started recording too late for ID22 in block 1 (missed 22 trials)        
        if strcmp(resp_cycles(i).ID,'22') && j == 1
            
            for k = 1:22
                l = l+1;
                stim_data(l,1:11) = [str2double(resp_cycles(i).ID) j k NaN NaN NaN NaN NaN NaN NaN NaN];
            end
        end
        
        % Loop trials
        for k = 1:length(resp_cycles(i).trigger(j).loc)
          
          % RESPIRATION CYCLE LOCKED TO INSPIRATION %  
            
          % Get index of inhale onset before stimulus onset (trigger)
          prev_inhale_onset_ind = find(resp_cycles(i).inhale_onsets(j).loc < resp_cycles(i).trigger(j).loc(k),1,'last');

          % If recording started too late for first trigger
          if isempty(prev_inhale_onset_ind)
              i
              j
              k
              prev_inhale_onset_ind
              diff2inhale_onset = NaN;
              resp_cycle_t_inhale = NaN;
              resp_cycle_t_next_inhale = NaN;
              stim_degree_inhale = NaN;
              inhale_duration_inhale = NaN;
              exhale_duration_inhale = NaN;
          
          else
              
              % Calculate temporal distance from stimulus onset to
              % preceding inspiration onset
              diff2inhale_onset = (resp_cycles(i).trigger(j).loc(k) - resp_cycles(i).inhale_onsets(j).loc(prev_inhale_onset_ind))/srate;

              % Calculate interval to subsequent trough (inhale onset)
              if (prev_inhale_onset_ind+1) <= length(resp_cycles(i).inhale_onsets(j).loc)
                resp_cycle_t_inhale = (resp_cycles(i).inhale_onsets(j).loc(prev_inhale_onset_ind+1) - resp_cycles(i).inhale_onsets(j).loc(prev_inhale_onset_ind))/srate;
              else
                i
                j
                k
                resp_cycle_t_inhale = NaN;
              end
              
              % Calculate interval of following respiratory cycle (inhale)
              if (prev_inhale_onset_ind+2) <= length(resp_cycles(i).inhale_onsets(j).loc)
                resp_cycle_t_next_inhale = (resp_cycles(i).inhale_onsets(j).loc(prev_inhale_onset_ind+2) - resp_cycles(i).inhale_onsets(j).loc(prev_inhale_onset_ind+1))/srate;
              else
                i
                j
                k
                resp_cycle_t_next_inhale = NaN;
              end

              % Calculate relative position of the stimulus onset to the respiration cycle
              stim_degree_inhale = (diff2inhale_onset/resp_cycle_t_inhale) * 360;
              
              % Determine inspiration and expiration duration,
              % only if respiration cylce (inhale) has been identified
              if ~isnan(resp_cycle_t_inhale)
                  
                % Check which exhale onset are inbetween the two inhale onsets defining the respiration cycle (inhale)  
                exhale_within_resp_cycle_inhale = intersect(resp_cycles(i).inhale_onsets(j).loc(prev_inhale_onset_ind):resp_cycles(i).inhale_onsets(j).loc(prev_inhale_onset_ind+1),...
                                                            resp_cycles(i).exhale_onsets(j).loc);

                % Only if one exhale onset lies inbetween, calculate
                % inhale/exhale durations
                if length(exhale_within_resp_cycle_inhale) == 1
                    inhale_duration_inhale = (exhale_within_resp_cycle_inhale - resp_cycles(i).inhale_onsets(j).loc(prev_inhale_onset_ind))/srate;
                    exhale_duration_inhale = (resp_cycles(i).inhale_onsets(j).loc(prev_inhale_onset_ind+1) - exhale_within_resp_cycle_inhale)/srate;
                else
                    disp('More inhale onsets within resp cycle (exhale)')
                    i
                    j
                    k
                    inhale_duration_inhale = NaN;
                    exhale_duration_inhale = NaN;
                end

              else
                  inhale_duration_inhale = NaN;
                  exhale_duration_inhale = NaN;
              end
                            
          end
          
          
          % RESPIRATION CYCLE LOCKED TO EXPIRATION %
          
          % Get index of exhale onset before stimulus onset (trigger)
          prev_exhale_onset_ind = find(resp_cycles(i).exhale_onsets(j).loc < resp_cycles(i).trigger(j).loc(k),1,'last');
          
          % If recording started too late for first trigger
          if isempty(prev_exhale_onset_ind)
              i
              j
              k
              prev_exhale_onset_ind
              diff2exhale_onset = NaN;
              resp_cycle_t_exhale = NaN;
              resp_cycle_t_next_exhale = NaN;
              stim_degree_exhale = NaN;
              inhale_duration_exhale = NaN;
              exhale_duration_exhale = NaN;
          
          else
              
              % Calculate temporal distance from stimulus onset to
              % preceding expiration onset
              diff2exhale_onset = (resp_cycles(i).trigger(j).loc(k) - resp_cycles(i).exhale_onsets(j).loc(prev_exhale_onset_ind))/srate;

              % Calculate interval to subsequent peak (exhale onset)
              if (prev_exhale_onset_ind+1) <= length(resp_cycles(i).exhale_onsets(j).loc)
                resp_cycle_t_exhale = (resp_cycles(i).exhale_onsets(j).loc(prev_exhale_onset_ind+1) - resp_cycles(i).exhale_onsets(j).loc(prev_exhale_onset_ind))/srate;
              else
                i
                j
                k
                resp_cycle_t_exhale = NaN;
              end
              
              % Calculate interval of following respiratory cycle (exhale)
              if (prev_exhale_onset_ind+2) <= length(resp_cycles(i).exhale_onsets(j).loc)
                resp_cycle_t_next_exhale = (resp_cycles(i).exhale_onsets(j).loc(prev_exhale_onset_ind+2) - resp_cycles(i).exhale_onsets(j).loc(prev_exhale_onset_ind+1))/srate;
              else
                i
                j
                k
                resp_cycle_t_next_exhale = NaN;
              end

              % Calculate relative position of the stimulus onset to the respiration cycle
              stim_degree_exhale = (diff2exhale_onset/resp_cycle_t_exhale) * 360;
              
              % Determine inspiration and expiration duration,
              % only if respiration cylce (exhale) has been identified
              if ~isnan(resp_cycle_t_exhale)
                  
                % Check which inhale onsets are inbetween the two exhale onsets defining the respiration cycle (exhale)  
                inhale_within_resp_cycle_inhale = intersect(resp_cycles(i).exhale_onsets(j).loc(prev_exhale_onset_ind):resp_cycles(i).exhale_onsets(j).loc(prev_exhale_onset_ind+1),...
                                                            resp_cycles(i).inhale_onsets(j).loc);

                % Only if one exhale onset lies inbetween, calculate
                % inhale/exhale durations
                if length(inhale_within_resp_cycle_inhale) == 1
                    exhale_duration_exhale = (inhale_within_resp_cycle_inhale - resp_cycles(i).exhale_onsets(j).loc(prev_exhale_onset_ind))/srate;
                    inhale_duration_exhale = (resp_cycles(i).exhale_onsets(j).loc(prev_exhale_onset_ind+1) - inhale_within_resp_cycle_inhale)/srate;
                else
                    disp('More exhale onsets within resp cycle (inhale)')
                    i
                    j
                    k
                    inhale_duration_exhale = NaN;
                    exhale_duration_exhale = NaN;
                end

              else
                  inhale_duration_exhale = NaN;
                  exhale_duration_exhale = NaN;
              end
          end

          % Store in matrix with 1 row per trial
          l = l+1;
          if strcmp(resp_cycles(i).ID,'03') && j == 2 % missed first 4 triggers in #2 block
            
            stim_data(l,1:15) = [str2double(resp_cycles(i).ID) j k+4 diff2inhale_onset resp_cycle_t_inhale stim_degree_inhale resp_cycle_t_next_inhale diff2exhale_onset resp_cycle_t_exhale stim_degree_exhale resp_cycle_t_next_exhale inhale_duration_inhale exhale_duration_inhale inhale_duration_exhale exhale_duration_exhale];          
            
          elseif strcmp(resp_cycles(i).ID,'22') && j == 1 % missed first 22 triggers in #1 block
            
            stim_data(l,1:15) = [str2double(resp_cycles(i).ID) j k+22 diff2inhale_onset resp_cycle_t_inhale stim_degree_inhale resp_cycle_t_next_inhale diff2exhale_onset resp_cycle_t_exhale stim_degree_exhale resp_cycle_t_next_exhale inhale_duration_inhale exhale_duration_inhale inhale_duration_exhale exhale_duration_exhale];         
              
          else
              
            stim_data(l,1:15) = [str2double(resp_cycles(i).ID) j k diff2inhale_onset resp_cycle_t_inhale stim_degree_inhale resp_cycle_t_next_inhale diff2exhale_onset resp_cycle_t_exhale stim_degree_exhale resp_cycle_t_next_exhale inhale_duration_inhale exhale_duration_inhale inhale_duration_exhale exhale_duration_exhale];          
          end
          
        end % loop triggers/trials
        
    end % loop files/blocks            
    
    % Add median to filter outlier respiration cycle participant-wise
    resp_cyle_t_inhale_median = median(stim_data(stim_data(:,1)==str2double(resp_cycles(i).ID),5),'omitnan');
    resp_cyle_t_exhale_median = median(stim_data(stim_data(:,1)==str2double(resp_cycles(i).ID),9),'omitnan');

    stim_data(stim_data(:,1)==str2double(resp_cycles(i).ID),16) = resp_cyle_t_inhale_median;
    stim_data(stim_data(:,1)==str2double(resp_cycles(i).ID),17) = resp_cyle_t_exhale_median;    
    
end % loop participants
