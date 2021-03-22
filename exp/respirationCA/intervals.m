function data = intervals(s,data)
% data = intervals(s,data) computes the actual time intervals between
% screen onsets (e.g., fixation to cue) and analog output trigger in run_exp.
%
% Input:
%   s        - run_exp settings structure (doc setup)
%   data     - run_exp output data structure (e.g., exp_data1)
%   
% Author:           Martin Grund
% Last update:      December 18, 2018

%%
for i = 1:size(data.seq,1)

%% SCREEN ONSET INTERVALS

    data.t_fix_cue(i,1) = (data.onset_cue{i,1}-data.onset_fix{i,1})*1000; % fix to cue interval
    data.t_cue_resp1(i,1) = (data.onset_resp1{i,1}-data.onset_cue{i,1})*1000; % fix to cue interval
    data.t_resp1_resp2(i,1) = (data.onset_resp2{i,1}-data.onset_resp1{i,1})*1000; % response 1 to response 2 screen interval

    % Trial
    if i < size(data.seq,1)
        data.t_trial_fix(i,1) = (data.onset_fix{i+1,1}-data.onset_fix{i,1})*1000;
        data.t_resp2_fix(i,1) = (data.onset_fix{i+1,1}-data.onset_resp2{i,1})*1000; % response 2 screen to fix interval
    else
        data.t_trial_fix(i,1) = (data.wait_block_end-data.onset_fix{i,1})*1000;
    end

  
%% CUE TO AO TRIGGER

    % (Cue onset to pre AO trigger) - stimulus delay
    data.t_cue_ao_trigger_pre_diff_stim_delay(i,1) = (data.ao_trigger_pre(i,1)-data.onset_cue{i,1}-data.seq(i,4))*1000;    
end

% STIMULUS ONSET LOCKED TO FIRST FIXATION ONSET
data.t_fix_stim_onset = data.ao_trigger_pre(:,1) + data.stim_offset/1000 - data.onset_fix{1,1};

data.t_trigger_ao = data.ao_trigger_post - data.ao_trigger_pre;