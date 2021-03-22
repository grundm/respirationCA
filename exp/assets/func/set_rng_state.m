function rng_state = set_rng_state(exp_settings)
% rng_state = set_rng_state(exp_settings) sets the random number generator
% state (rng state). It checks if a rng state is defined in a field 
% 'rng_state' of the experiment settings structure (e.g. "nt.rng_state")
% and uses this one. If not, it generates an unique one based on the
% current time.
%
% Author:           Martin Grund
% Last update:      November 19, 2015


%% Check if a state is defined. If not, generate one.
if isfield(exp_settings, 'rng_state')
    rng_state = exp_settings.rng_state;
    disp(['Random number generator state set to ' inputname(1) '.rng_state.']);
else
    rng_state = sum(100*clock);
    disp('Random number generator state set to unique time point.');
end

%% Initialize random generator (Shuffle uses rand)
rand('state',rng_state);