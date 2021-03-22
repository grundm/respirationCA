function [ao,ai] = aio_setup
% [ao,ai] = aio_setup returns an analog output and input object for the data
% acquisition card from Data Translation (DT-9812)
% 
% The settings (e.g., sample rate) have to be edited within this function.
% 
% Author:           Martin Grund
% Last update:      November 5, 2015

%% Helpful commands
% - daqhwinfo(ai)
% - propinfo(ai)
% - daqreset

%% Stop any data acquisition objects

if (~isempty(daqfind))
    stop(daqfind)
end

%% Analog output setup

% create output object
ao = analogoutput('dtol',0); % -> Warning: External clock not suported

% Settings
ao.SampleRate = 50000;
ao.TriggerType = 'Immediate';%'Manual';
ao.BufferingMode = 'Auto'; % match MATLAB buffer to DT buffer and the buffer in the PutData method (?)

% DAQ Adaptor shows bug for DT912
% Add channels after setting SampleRate, because
% - error for ao.SampleRate = 50000; if 2 channels are added,
% - likely ao.SampleRate is doubled for each channel after adding channels

% set channels
addchannel(ao,0:1); % -> Warning: Random channel assignment not supported

% zero analog output
% BUG: Sets ao.SampleRate = 25000
% putsample(ao,repmat(0,1,size(ao.channel,1)));

% THESE SETTINGS CAUSE MATLAB TO CRASH
% Reset output value to zero after output is done
% ao.OutOfDataMode = 'DefaultValue'; % default: 'Hold'
% ao.Channel.DefaultChannelValue = 0;

% Further possible seetings

% number of samples in each waveform used in PutData
% set(ao,'BufferingConfig',[samplesPerBuffer,10]);


%% Input settings

% create input object
ai = analoginput('dtol',0);

% set channel
addchannel(ai,0:1);

% Settings
ai.SampleRate = 25000;
ai.TriggerType = 'Immediate';
ai.SamplesPerTrigger = Inf;
ai.BufferingMode = 'Auto';

% Data logging
ai.LogFileName = [tempdir 'ai_rec_' datestr(now,'yyyy-mm-dd-HH-MM') '_01'];
ai.LogToDiskMode = 'Index';
ai.LoggingMode = 'Disk'; %Memory';
