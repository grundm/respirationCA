function plot_PM_run(PM,finger_label)
% plot_PM_run(PM,finger_label) creates a figure of the applied intensities 
% and threshold estimations over the run of the psi method block.
%
% Additionally, the figure indicates 'yes' responses and the stimulus 
% range, as well as the prior alpha range.
%
% Input:
%   PM              - output structure by psi method (PAL_AMPM_updatePM.m)
%   finger_label    - text string to indicate the finger
%
% Author:           Martin Grund
% Last update:      July 11, 2018

Fig_PM = figure;
set(Fig_PM,'Name',['Finger ' finger_label ': Psi method estimation over trials']);

hold on;

plot([0 length(PM.x)],[min(PM.stimRange) min(PM.stimRange)],'Color',[.3 .3 .3],'LineStyle',':');
plot([0 length(PM.x)],[max(PM.stimRange) max(PM.stimRange)],'Color',[.3 .3 .3],'LineStyle',':');
               
rectangle('Position',[.1 min(PM.priorAlphaRange) length(PM.x) max(PM.priorAlphaRange)-min(PM.priorAlphaRange)],...
          'FaceColor',[.9 .9 .9],...
          'LineStyle','none');
      
text(.5,max(PM.stimRange)+.15,'PM.stimRange');
text(.5,max(PM.priorAlphaRange)-.15,'PM.priorAlphaRange');

p_x = plot(1:length(PM.x),PM.x,'k-');
p_t = plot(1:length(PM.threshold),PM.threshold,'m-');
% p_s = plot(1:length(PM.slope),10.^PM.slope);
b_r = bar(1:length(PM.response),PM.response/2,'FaceColor',[.6 .6 .6],'LineStyle','none');

% ylim([0 max(10.^PM.slope)+1]);
ylim([0 max(PM.stimRange)+1]);
xlim([0 length(PM.x)]);

xlabel('Trial');
ylabel('Slope / Intensity in mA');

legend([p_x p_t b_r],...
       {'PM.x',...
       'PM.threshold',...
       'PM.response (yes)'},...
       'Location','NorthEast');

% legend([p_x p_t p_s b_r],...
%        {'PM.x',...
%        'PM.threshold',...
%        'PM.slope',...
%        'PM.response (yes)'},...
%        'Location','NorthWest');

hold off;