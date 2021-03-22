function print_UD_res(UD_res,PF_params_PM,finger_label,block,run)
% print_UD_res(UD_res,PF_params_PM,finger_label,block,run) prints the 
% range, mean, 50% threshold estimate and detection rates for the up/down
% method, as well as the 50% threshold estimate of the psi method
%
% Input:
%   UD_res          - results of up/down method in thr2F_run.m
%   PF_param_PM     - psyometric function parameter based on psi method
%   finger_label    - text string to indicate the finger
%
% Author:           Martin Grund
% Last update:      July 11, 2018

disp(sprintf(['\nFinger ' finger_label '\nThA #' num2str(block) '-' num2str(run) '\n']));
disp(['UD range: ' num2str(UD_res.UD_range(1)) '-' num2str(UD_res.UD_range(2))]);
disp(['UD mean: ' num2str(UD_res.UD_mean)]);
disp(['UD PF 50%: ' num2str(UD_res.PF_params_UD(1))]);
disp(UD_res.x_resp_freq_UD);

disp(['PM PF 50%: ' num2str(PF_params_PM(1))]);