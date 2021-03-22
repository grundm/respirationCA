function thr1F_data = thr1F_loop(thr1F,p_data,aio_s,block,run)
% thr1F_data = thr1F_loop(thr1F,p_data,aio_s,block,run) runs threshold 
% assessment with thr1F_run, saves its output (thr1F_data) to the 
% participants data directory and plots the results of the assessment on 
% the modelled psychometric function.
%
% Input:
%   thr1F           - settings structure (doc thr1F_setup)
%   p_data          - output of participant_data
%   aio_s           - daq acquisition session object
%   block           - block number for filename end
%   run             - run number for filename end
%
% Author:           Martin Grund
% Last update:      December 18, 2018


%%

    % Run threshold assessment
    thr1F_data = thr1F_run(thr1F,aio_s,p_data.ID);
    
    % Compute intervals
    thr1F_data = thr1F_intervals(thr1F_data);

    % Save threshold assessment data
    thr1F_save(p_data,thr1F_data,thr1F,['0' num2str(block) '_0' num2str(run)]);
    
    % Plot up/down trials
    plot_UD_run(thr1F_data.UD,'1');
    
    % Print up/down & psi results
    print_UD_res(thr1F_data,thr1F_data.PF_params_PM,'1',block,run)    
    
    % Plot psi method trials
    plot_PM_run(thr1F_data.PM,'1');
    
    % Plot estimated psychometric function with test results
    plot_PM_PF(thr1F_data.PM,thr1F_data,1,1,'1');