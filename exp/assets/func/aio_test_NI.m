function ai_rec = aio_test_NI(aio_s,ao,ai,stim,scale_out)
% ai_rec = aio_test_NI(aio_s,ao,ai,stim,scale_out) sends a rectangular
% pulse (data vector "stim") to all defined analog output channels and 
% records on all defined analog input channels.
%
% It plots the acquired data and displays the number of acquired samples 
% that are about the threshold (mean value of channel + 25% rectangular 
% pulse height). The recorded analog input data is returned in ai_rec.
%
% scale_out - allows to scale the recorded data for plotting and pulse
%             detection, e.g. DS5 output is 1V = 10mA
%
% Author:           Martin Grund
% Last update:      July 4, 2018

ao_data = repmat(stim,1,length(ao));

stop(aio_s);

queueOutputData(aio_s,ao_data);

[ai_rec.data,ai_rec.time] = startForeground(aio_s);

stop(aio_s);

plot(ai_rec.time,ai_rec.data.*scale_out);

disp(['stim: ' num2str(numel(find(stim>.5*max(stim)))) ' sent -> ' num2str(numel(find(stim>.5*max(stim)))) ' expected']);

for i = 1:length(ai)
    disp(['AI' num2str(i) ': ' num2str( numel( find((ai_rec.data(:,i).*scale_out)>(mean(ai_rec.data(:,i).*scale_out)+.25*max(stim))) ) ) ]);
end