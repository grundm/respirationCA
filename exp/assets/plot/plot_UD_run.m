function plot_UD_run(UD,finger_label)
% plot_UD_run(UD,finger_label) creates a figure of the applied intensities
% over the run of the up/down method block.
%
% Additionally, the figure indicates 'yes' responses.
%
% Input:
%   UD              - output structure by up/down method (PAL_AMUD_updateUD.m)
%   finger_label    - text string to indicate the finger
%
% Author:           Martin Grund
% Last update:      July 11, 2018

Fig_UD = figure;
hold on;
set(Fig_UD,'Name',['Finger ' finger_label ': Up/down method trials']);

plot(1:length(UD.x),UD.x);

bar(1:length(UD.x),UD.response/2,'FaceColor',[.6 .6 .6],'LineStyle','none');

hold off;