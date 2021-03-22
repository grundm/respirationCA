function log_detection(thr1F_data,exp_data)
% detection_log(thr1F_data,exp_data) displays the applied intensities and 
% their detection rates in the experimental block, as well as the expected
% detection rates based on the estimated psychometric function in thr1F_run.
%
% Input:
%   thr1F_data      - output of threshold assessment thr1F_run
%   exp_data        - output of experimental block (run_exp)
%
% Author:           Martin Grund
% Last update:      December 18, 2018

% Display block number
disp(['Block #' num2str(exp_data.seq(1,1))]);

% Calculate detection rates for all intensities in experiment
nt_detection = count_resp([exp_data.intensity exp_data.resp1]);
    
% Display expected detection rates
arrayfun(@(intensity) disp(['PF(' num2str(intensity) ' mA) = ' num2str(PAL_Quick(thr1F_data.PF_params_PM,intensity))]), nt_detection(:,1));
    
% Display actual detection rates
disp(nt_detection);