function [aio_s,ao,ai] = aio_setup_NI(ao_num,ai_num,s_rate)
% s = aio_setup returns a session (aio_s), analog output (ao) and input 
% object (ai) for the data acquisition card from National Instruments 
% (e.g., USB-6343)
% 
% Input
%   ao_num         - number of analog output channels
%   ai_num         - number of analog input channels
%   s_rate         - sampling rate (NI USB-6343: max. 500 kHz / ai_num)
% 
% Author:           Martin Grund
% Last update:      July 4, 2018

%% Helpful commands
% - daqhwinfo(ai)
% - propinfo(ai)
% - daqreset

%% Stop any data acquisition objects

if (~isempty(daqfind))
    stop(daqfind)
end

%% Create session

aio_s = daq.createSession('ni');

% Settings
aio_s.Rate = s_rate;

% Open questions:
% How do we set input/output rate if it is only possible for the session?
% -> Only for session

%% Analog output setup

% Add channels
ao = addAnalogOutputChannel(aio_s,'Dev1',0:ao_num-1,'Voltage');

%% Analog input setup

% Add channels
ai = addAnalogInputChannel(aio_s,'Dev1',0:ai_num-1,'Voltage');

% Settings
for i = 1:ai_num
    ai(i).Range = [-10.0 10.0];
end
%ai.TerminalConfig = 'SingleEnded';

% ai.SampleRate = 25000;
% ai.TriggerType = 'Immediate';
% ai.SamplesPerTrigger = Inf;
% ai.BufferingMode = 'Auto';

% Work-around to log data continously to file
% https://de.mathworks.com/help/daq/examples/log-analog-input-data-to-a-file-using-ni-devices.html
%
% For now I will save analog input data only while analog output data is
% send

% % Data logging
% ai.LogFileName = [tempdir 'ai_rec_' datestr(now,'yyyy-mm-dd-HH-MM') '_01'];
% ai.LogToDiskMode = 'Index';
% ai.LoggingMode = 'Disk'; %Memory';
