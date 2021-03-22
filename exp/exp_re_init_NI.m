function [aio_s,ao,ai] = exp_re_init_NI(exp_dir,thr_dir,p_data,ao_num,ai_num,s_rate)
% [aio_s,ao,ai] = exp_re_init(exp_dir,p_data) runs initial procedures for 
% experiment:
%   - sets paths
%   - starts diary
%	- aio_setup
%
% Author:           Martin Grund
% Last update:      December 13, 2018

%%
% Make all assets available (e.g., Palamedes toolbox)
addpath(genpath([pwd, '/', exp_dir]))
addpath(genpath([pwd, '/', thr_dir]))
addpath(genpath([pwd, '/assets']))

%% Diary logfile   
diary([p_data.dir 'exp_' p_data.ID '_log.txt']);

%% Setup analog output (ao) and input (ai)
[aio_s,ao,ai] = aio_setup_NI(ao_num,ai_num,s_rate);
