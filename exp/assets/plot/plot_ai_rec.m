function plot_ai_rec(time,data)
% plot_ai_rec(time,data) creates a figure with the recording of two analog
% input channels (e.g., current and voltage output by Digitimer DS5).
%
% Author:           Martin Grund
% Last update:      November 10, 2015

%% Plot analog input recording
% voltage - DS5 voltage signal of 1V = 10V
% current - DS5 signal of 1V = 10 mA

fig = figure;
set(fig,'Name','Analog input recording');

plot(time+1,data(:,1)*10);  % DS5 current output
hold on;
plot(time,data(:,2)*1,'g'); % analog output channel 2

xlim([0 max(time)+5]);
xlabel('Time in s');
ylabel('Intensity in mA');

legend('DS5 current',...
       'Analog ouput channel 2',...
       'Location','SouthEast');

hold off;