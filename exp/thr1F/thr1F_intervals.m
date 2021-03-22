function data = thr1F_intervals(data)
% data = thr1F_intervals(data) computes the actual time intervals
% between screen onsets (e.g., fixation to cue) and analog output trigger 
% in thr1F_run.
%
% Input:
%   data     - thr1F_run output data structure (e.g., thr1F_data)
%   
% Author:           Martin Grund
% Last update:      December 18, 2018

for i = 1:size(data.seq,1)

%% SCREEN ONSET INTERVALS

    data.t_fix_cue(i,1) = (data.onset_cue{i,1}-data.onset_fix{i,1})*1000; % fix to cue interval
    data.t_cue_resp(i,1) = (data.onset_resp{i,1}-data.onset_cue{i,1})*1000; % cue to response screen    
    if i < size(data.seq,1)
        data.t_trial_fix(i,1) = (data.onset_fix{i+1,1}-data.onset_fix{i,1})*1000; % Trial
    else
        data.t_trial_fix(i,1) = (data.onset_resp{i,1} + data.resp_t(i,1) - data.onset_fix{i,1})*1000; % Trial
    end
    
%% AO TRIGGER    
    
    % (Cue onset to pre AO trigger) - stimulus delay [pre AO trigger is MRI trigger locked]
    data.t_cue_ao_trigger_pre_diff_stim_delay(i,1) = (data.ao_trigger_pre(i,1)-data.onset_cue{i,1}-data.stim_delay(i))*1000;    
end

data.t_trigger_ao = data.ao_trigger_post - data.ao_trigger_pre;