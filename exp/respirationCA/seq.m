function [seq,s] = seq(s)
% [seq,nt] = seq(s) returns a sequence matrix (seq) that is shuffled 
% per block as well as the random number generator state that was set.
%
%   seq: block x stimulus type x stimulus step x stimulus delay
%
% The number of stimulus delays is equal the number of trials per block,
% defined as sum(nt.stim_types_num) and lineraly spaced between the minimum
% and maximim of the defined stimulus delay interval (s.stim_delay).
%
% Relevant settings in nt_setup (doc setup):
%   s.rng_state
%   s.blocks
%   s.stim_types
%   s.stim_types_num
%   s.stim_steps
%   s.stim_delay
%   s.fix_t
%
% Author:           Martin Grund
% Last update:      December 18, 2018

%%

% Set random number generator state
s.rng_state = set_rng_state(s);

seq = [];

for i = 1:s.blocks
    
    seq_type_tmp = [];
    seq_step_tmp = [];
    
    % Loop stimulus types
    for j = 1:numel(s.stim_types)
        seq_type_tmp = [seq_type_tmp; ones(s.stim_types_num(j),1)*s.stim_types(j);];   
    
        % Vector stimulus steps
        if s.stim_types(j) == 0
            steps = zeros(s.stim_types_num(j),1);
        else
            if i == 1 && mod(s.stim_types_num(j),numel(s.stim_steps)) ~= 0
                warning('nt_exp:trialNum',...
                        ['Sequence generation - No equal stimulus step distribution ',...
                         'for stimulus type "' num2str(s.stim_types(j)) '", because the ',...
                         'frequency of this stimulus type (' num2str(s.stim_types_num(j)) ') ',...
                         'is not a multiple of the number of stimulus steps (',...
                         num2str(numel(s.stim_steps)) ').']);
            end            
            steps = repmat(s.stim_steps',ceil(s.stim_types_num(j)/numel(s.stim_steps)),1);
        end

        seq_step_tmp = [seq_step_tmp; steps(1:s.stim_types_num(j))];

    end
    
    % Shuffle stimulus types
    seq_ind = Shuffle(1:sum(s.stim_types_num));
    
    % Shuffle stimlus delays (rounded on ms)
    seq_stim_delay = Shuffle(round_dec(linspace(s.stim_delay(1),s.stim_delay(2),sum(s.stim_types_num)),4))';
    
    % Shuffle fix_t (rounded on ms)
    seq_fix_t = Shuffle(round_dec(linspace(s.fix_t(1),s.fix_t(2),sum(s.stim_types_num)),4))';
    
    % Sequence matrix: block - type - step - delay - fix_t
    seq = [seq; i*ones(sum(s.stim_types_num),1) seq_type_tmp(seq_ind) seq_step_tmp(seq_ind) seq_stim_delay seq_fix_t];        
end