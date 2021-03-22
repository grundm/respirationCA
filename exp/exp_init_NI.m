function [p_data,aio_s,ao,ai] = exp_init_NI(exp_dir,thr_dir,ao_num,ai_num,s_rate)
% [p_data,ao,ai] = exp_init(exp_dir) runs initial procedures for 
% experiment:
%   - sets paths
%   - starts diary
%   - participant_data
%	- aio_setup
%
% Author:           Martin Grund
% Last update:      May 9, 2018

%%
% Make all assets available (e.g., Palamedes toolbox)
addpath(genpath([pwd, '/', exp_dir]))
addpath(genpath([pwd, '/', thr_dir]))
addpath(genpath([pwd, '/assets']))

%% Particpant data
p_data = participant_data(['data/', exp_dir, '/ID']);

%% Diary logfile   
diary([p_data.dir 'exp_' p_data.ID '_log.txt']);

%% Setup analog output (ao) and input (ai)
[aio_s,ao,ai] = aio_setup_NI(ao_num,ai_num,s_rate);
