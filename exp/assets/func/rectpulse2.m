function [stim,stim_offset] = rectpulse2(pulse_t,I,sample_rate,pre_pulse_t,wave_t)
% [stim,stim_offset] = rectpulse2(pulse_t,I,sample_rate,pre_pulse_t,wave_t)
% generates a single rectangular pulse with the intensity I for the duration 
% pulse_t starting after pre_pulse_t and ending after wave_t
%
% Input:
%   pulse_t         - duration of rectangular pulse in ms
%   I               - intensity/height of rectangular pulse in mA
%   sample_rate     - sampling rate of analog output channel in Hz
%   pre_pulse_t     - duration before rectangular pulse in ms
%   wave_t          - duration of waveform in ms
%
% Output:
%   stim            - waveform vector with rectangular pulse
%   stim_offset     - time in ms when rectangular pulse ends in waveform
%
% Author:           Martin Grund
% Last update:      July 12, 2018

%% Calculate number of sample pulse and pre/post pulse interval
sample_t = 1000/sample_rate; % sample duration in ms

pulse_samples = pulse_t/sample_t; % samples per pulse
pre_pulse_samples = pre_pulse_t/sample_t;
wave_samples = wave_t/sample_t;

%% Check if durations are positive integers and can be sampled.
if pulse_samples < 0 || mod(pulse_samples,1)~=0
    error('Pulse duration does not match sampling rate (cycle) or is negative');
end

if pre_pulse_samples < 0 || mod(pre_pulse_samples,1)~=0
    error('Pre-pulse duration does not match sampling rate (cycle) or is negative');
end

if wave_samples < 0 || mod(wave_samples,1)~=0
    error('Waveform duration does not match sampling rate (cycle) or is negative');
end

%% Construct stimulus waveform

stim = [zeros(pre_pulse_samples,1); ones(pulse_samples,1)*I; zeros(wave_samples - pre_pulse_samples - pulse_samples,1);];

stim_offset = pre_pulse_t + pulse_t;
