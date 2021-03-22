function peaks = eval_DS5_output(time,data,intensities)
% eval_DS5_output(data,time,intensities) compares the DS5 current output 
% with administered intensities
%
% Requires Signal Processing Toolbox
%
% Input
%   time        - time of recorded analog data
%   data        - recorded analog data of two channels
%   intensities - administered intensities (e.g., nt_data.intensity)
%
% Author:           Martin Grund
% Last update:      January 5, 2016


% Find peaks in DS5 current output and analog output channel 2
[peaks_DS5,locs_DS5] = findpeaks(data(:,1)*10,time,'MinPeakHeight',1.0,'MinPeakDistance',.5);
[peaks_ao2,locs_ao2] = findpeaks(data(:,2),time,'MinPeakHeight',1.0,'MinPeakDistance',.5);

% Get administered intensities and attach measured peaks
peaks = [intensities(intensities>0) peaks_DS5 peaks_ao2];

% Compute peak difference for the different recordings
% I-DS5 & I-AO2 & DS5-AO2
peaks(:,4:6) = [peaks(:,1)-peaks(:,2) peaks(:,1)-peaks(:,3) peaks(:,2)-peaks(:,3)];

% Plot administered intensities, measured peaks from DS5 and their
% difference
fig = figure;
set(fig,'Name','Difference of DS5 output and intensities');
hold on;

% Plot administered intensities
plot(1:length(peaks(:,1)),peaks(:,1),'k');

% Plot measured DS5 current peaks
plot(1:length(peaks(:,1)),peaks(:,2),'b');

% Plot measured analog output channel 2 peaks
plot(1:length(peaks(:,1)),peaks(:,3),'m--');

% Plot difference between intensities and DS5 current peaks
plot(1:length(peaks(:,1)),peaks(:,4),'r');

% Plot difference between intensities and AO2 current peaks
plot(1:length(peaks(:,1)),peaks(:,5),'c--');

hold off;

xlim([0 size(peaks,1)+1]);

xlabel('Trials');
ylabel('Current in mA');

legend('Intensities',...
       'DS5 current peaks',...
       'AO2 peaks',...
       'I-DS5 difference',...
       'I-AO2 difference',...
       'Location','EastOutside');

header = {'Intensities','DS5','AO2','I-DS5','I-AO2','DS5-AO2'};

disp(header);
disp(peaks);
disp(['Int-DS5 mean (sd): ' num2str(mean(abs(peaks(:,4)))) ' (' num2str(std(peaks(:,4))) ')']);
disp(['Int-AO2 mean (sd): ' num2str(mean(abs(peaks(:,5)))) ' (' num2str(std(peaks(:,5))) ')']);
disp(['DS5-AO2 mean (sd): ' num2str(mean(abs(peaks(:,6)))) ' (' num2str(std(peaks(:,6))) ')']);