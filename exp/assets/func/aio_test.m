function ai_rec = aio_test(ao,ai,stim,scale_out)
% ai_rec = aio_test(ao,ai,stim) sends a rectangular pulse (data vector 
% "stim") to all defined analog output channels and records on all defined 
% analog input channels.
%
% It plots the acquired data and displays the number of acquired samples 
% that are about the threshold (mean value of channel + 25% rectangular 
% pulse height). The recorded analog input data is returned in ai_rec.
%
% scale_out - allows to scale the recorded data for plotting and pulse
%             detection, e.g. DS5 output is 1V = 10mA
%
% After reading the analog input logfile it is deleted.
%
% Author:           Martin Grund
% Last update:      November 16, 2017

data = repmat(stim,1,length(ao.channel));

stop([ao ai]);

putdata(ao,data);

start(ai);

f=ai.LogFileName;

pause(1);

start(ao);
wait(ao,.5);

pause(1)

stop(ai);

[ai_rec.data, ai_rec.time, ai_rec.abstime, ai_rec.events, ai_rec.daqinfo] = daqread(f);

delete(f);

plot(ai_rec.time,ai_rec.data.*scale_out);

disp(['stim: ' num2str(numel(find(stim>.5*max(stim)))) ' sent -> ' num2str((ai.SampleRate/ao.SampleRate)*numel(find(stim>.5*max(stim)))) ' expected']);

for i = 1:length(ai.channel)
    disp(['ch' num2str(i) ': ' num2str( numel( find((ai_rec.data(:,i).*scale_out)>(mean(ai_rec.data(:,i).*scale_out)+.25*max(stim))) ) ) ]);
end