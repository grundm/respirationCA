function [stim,offset] = rectpulse(t,I,sample_rate,wave_t)
% [stim,offset] = rectpulse(t,I,sample_rate,wave_t) generates a single 
% rectangular pulse with the intensity I for the duration t starting after 
% half the wave form duration wave_t
%
% Input:
%   - t in ms
%   - sampling rate in Hz
%   - wave_t in ms
%
% Author:           Martin Grund
% Last update:      December 8, 2015

    sample_t = 1000/sample_rate; % sample duration in ms
    pulse_samples = t/sample_t; % samples per pulse
    half_wave_samples = (0.5*wave_t)/sample_t; % half samples per wave
    offset = 0.5*wave_t;
    
    % Check if durations are positive integers and can be sampled.
    if pulse_samples < 0 || mod(pulse_samples,1)~=0
        error('Pulse duration does not match sampling rate (cycle) or is negative');
    end
    
    if half_wave_samples < 0 || mod(half_wave_samples,1)~=0
        error('Wave duration does not match sampling rate (cycle) or is negative');
    end
    
    stim = [zeros(half_wave_samples,1); ones(pulse_samples,1)*I; zeros(half_wave_samples-pulse_samples,1);];