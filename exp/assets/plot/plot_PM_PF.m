function PM_PF_params = plot_PM_PF(PM,thr_data,finger_num,finger_index,finger_label)
% plot_PM_PF(PM,thr_data,finger_index,finger_label) estimates and plots the
% psychometric function (PF) based on the psi method (PM) used in psi_1AFC.m
%
% It creates a figure with the PF and maps the actual applied intensities
% while the PM and test block on the PF, plus the measured performance of
% the test intensities.
%
% Input:
%   PM              - output structure by psi method (PAL_AMPM_updatePM.m)
%   thr_data        - output structure by psi_1AFC.m or thr2F_run.m
%   finger_num      - number of fingers
%   finger_index    - integer to indicate the finger (e.g., 1 or 2)
%   finger_label    - text string to indicate the finger
%
% Output:
%   PM_PF_params    - vector of PF threshold, slope, guess and lapse rate
%
% Author:           Martin Grund
% Last update:      December 17, 2018

%% Psychometric function

PM_PF_params = [PM.threshold(end); 10.^PM.slope(end); PM.guess(end); PM.lapse(end);]; % [alpha beta gamma lambda]
PF_stims = 0:.01:max([thr_data.supra(1,finger_index); thr_data.intensity(:,finger_index)])+1; % PM.stimRange

% Create intensities x responses matrices
if finger_num == 1
    x_resp = [thr_data.intensity thr_data.resp];
    x_resp_test = x_resp(end-size(thr_data.seq_test,1)+1:end,:);
elseif finger_num == 2
    x_resp = [thr_data.intensity(thr_data.seq(:,2)==finger_index,finger_index) thr_data.resp(thr_data.seq(:,2)==finger_index)];
    x_resp_test = x_resp(end-size(thr_data.seq_test(thr_data.seq_test(:,1)==finger_index,:),1)+1:end,:);    
end

x_resp_psi = [PM.x(1:end-1)' PM.response'];

% Compute response frequencies
x_resp_freq_psi = count_resp(x_resp_psi);
x_resp_freq_test = count_resp(x_resp_test)

% Prepare figure
Fig_PF = figure;
set(Fig_PF,'Name',['Finger ' finger_label ': Psychometric function']);
ylabel('Proportion of "yes"');
xlabel('Intensity in mA');
ylim([0 1]);
xlim([min(PF_stims) max(PF_stims)]);
hold on;

% Plot psychometric function (PF)
plot(PF_stims,PAL_Quick(PM_PF_params,PF_stims));

% Plot applied psi method intensities on psychometric function
plot(PM.x,PAL_Quick(PM_PF_params,PM.x),'k+');

% Plot test intensities on psychometric function
plot(x_resp_freq_test(:,1),PAL_Quick(PM_PF_params,x_resp_freq_test(:,1)),'r+');

% Plot real performance of test intensities
plot(x_resp_freq_test(:,1),x_resp_freq_test(:,4),'ro');

% Plot constant for p = .5
plot(get(gca,'xlim'),[.5 .5],'Color',[.7 .7 .7]);

% Plot constant for estimated threshold
plot([PM_PF_params(1) PM_PF_params(1)],get(gca,'ylim'),'g-');

% Plot constant for supra-threshold intensity
plot([thr_data.supra(1,finger_index) thr_data.supra(1,finger_index)],get(gca,'ylim'),'m-');


legend('Psychometric function (Quick)',...
       'Applied intensities',...
       'Test intensities',...
       'Performance test intensities',...
       '50%',...
       'Threshold estimation',...
       ['Supra-threshold (~' num2str(PAL_Quick(PM_PF_params,thr_data.supra(1,finger_index))*100) '%)'],...
       'Location','NorthWest');

hold off;