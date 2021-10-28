% This script calculates hierachichal meta d-prime (type 2 SDT) according to
% Fleming et al. (2019)

% Author:           Martin Grund, Carina Forster
% Last update:      October 25, 2021

%% Settings

code_dir = '/Users/martin/ownCloud/promotion/experiment/respirationCA/code/respirationca';

data_dir = '/Users/martin/ownCloud/promotion/experiment/respirationCA/data';


%% Set paths
cd(code_dir)

% Add paths which store HMeta-d scripts (for this script we need 
% trails2counts.m and fit_mcmc_group)
addpath(genpath([code_dir '/assets/HMeta-d']));

%% Prepare data for Mratio calculation

% Load and prepare data
data = readtable([data_dir '/trials_ecg_bin_20210804']);

% IMPORTANT:
% Make sure that confidence is coded as 1,2 and not 0,1
data.conf(data.resp2==1)=2;
data.conf(data.resp2==0)=1;


%% Count trials for cardiac cycle bins
 
IDs = unique(data.ID);

% Loop cardiac phases
for k = 1:max(data.dist)

    trial_counts(k).nR_S1 = zeros(length(IDs),4);
    trial_counts(k).nR_S2 = zeros(length(IDs),4);

    % Loop all participants to get confidence vector for meta d prime
    for i = 1:length(IDs)
    
        data_tmp = data(data.ID == IDs(i) & data.dist == k,:);
         
        % Input: you need columns that contain stimulus,
        % response, confidence (converted to 1 and 2),
        % rating, number of confidence levels
        % S1 = catch trials, no response
        % S2 = signal trials, yes response
        
        [trial_counts(k).nR_S1(i,:), trial_counts(k).nR_S2(i,:)] = trials2counts(data_tmp.stim_type, data_tmp.resp1, data_tmp.conf, 2, 0);
         
    end

    % Convert to cell array for meta d prime function
    trial_counts(k).nR_S1_cell = mat2cell(trial_counts(k).nR_S1, ones(1, size(trial_counts(k).nR_S1, 1)), size(trial_counts(k).nR_S1, 2))';
    trial_counts(k).nR_S2_cell = mat2cell(trial_counts(k).nR_S2, ones(1, size(trial_counts(k).nR_S2, 1)), size(trial_counts(k).nR_S2, 2))';

end

%% Report number of trials per cardicac cycle bin

% No responses
for k = 1:4
    no_resp(:,k) = sum(trial_counts(k).nR_S1(:,1:2) + trial_counts(k).nR_S2(:,1:2),2);
    yes_resp(:,k) = sum(trial_counts(k).nR_S1(:,3:4) + trial_counts(k).nR_S2(:,3:4),2);
end

mean(no_resp)
std(no_resp)
mean(yes_resp)
std(yes_resp)


%% Fit M-ratio for each cardiac cycle bin separately

% Get default parameters
mcmc_params = fit_meta_d_params;

% Change default to estimate d-prime
mcmc_params.estimate_dprime = 1;

% Change defaults to make response-conditional
mcmc_params.response_conditional = 1;

% Fit group data all at once for each phase

for k = 1:length(trial_counts)

    mratio.fit(k) = fit_meta_d_mcmc_group(trial_counts(k).nR_S1_cell, trial_counts(k).nR_S2_cell, mcmc_params);

end

% save([data_dir, '/fit_mratio_rs.mat'],'mratio');


%%
figure;
title_str = {'0-200 ms', '200-400 ms', '400-600 ms', '600-800 ms'};

i = 1;

%%% Posterior distribution of group-level M-ratios for yes/no-responses %%%
% + boxplots of participant-level M-ratios
% + posterior distribution of difference between yes/no

% Loop cardiac intervals
for k = 1:4
    
    % Posterior distributions yes/no + boxplots (columns 1:2 of 3)
    subplot(5,3,i:(i+1));
    i = i + 2;
        
    % Posterior distribution of group-level M-ratio for no/yes-responses
    yyaxis left
    histogram(exp(mratio.fit(k).mcmc.samples.mu_logMratio_rS1(:,:,1)),'FaceColor','#e78ac3','EdgeColor','none')
    hold on
    histogram(exp(mratio.fit(k).mcmc.samples.mu_logMratio_rS2(:,:,1)),'FaceColor','#66c2a5','EdgeColor','none')
    
    ylabel('Samples')
    %xlabel('M-ratio')
    xlim([0.1,1.9])    
    ylim([0,1800])

    if k == 1; legend({'no','yes'}); end

    % Boxplots of individual-level M-ratio
    yyaxis right
    boxplot([mratio.fit(k).Mratio_rS1', mratio.fit(k).Mratio_rS2'],'Orientation','horizontal','Colors','k','Symbol','.k')
    set(findobj(gcf,'LineStyle','--'),'LineStyle','-')    

    % Hard-coded to shift boxplots a bit down
    ylim([0.5, 3.5])

    % Title for cardiac cycle bin
    title(['M-ratio ' title_str{k}])
    
    % High-density interval (HDI) of posterior difference (column 3 of 3)
    subplot(5,3,i)
    i = i + 1;
    
    % Calculate HDI of posterior difference
    diff_tmp = mratio.fit(k).mcmc.samples.mu_logMratio_rS2(:,:,1)-mratio.fit(k).mcmc.samples.mu_logMratio_rS1(:,:,1);
    hdi_tmp = calc_HDI(diff_tmp(:));

    % Plot posterior difference
    histogram(diff_tmp,'BinLimits',hdi_tmp,'BinWidth',0.03,'EdgeColor','none','FaceColor','#8da0cb')
    hold on
    
    % Vertical line at zero (critical value for HDI)
    xline(0,'Color','#A2142F');
    
    xlim([-0.5, 1.1])
    ylim([0, 1800])
    

    title('Diff(yes-no) 95%-HDI')

end

%%% 95-% HDI of posterior difference between cardiac intervals %%%

% Define colors for the three time interval contrasts
edge_cols = {'k','#EDB120','#7E2F8E'}; % black, orange, purple


%%% HDIs for no-responses %%%
subplot(5,3,i)
i = i + 1;

% Loop time interval contrasts
for k = 1:3       

    diff_tmp = mratio.fit(k+1).mcmc.samples.mu_logMratio_rS1(:,:,1)-mratio.fit(k).mcmc.samples.mu_logMratio_rS1(:,:,1);
    hdi_tmp = calc_HDI(diff_tmp(:));

    histogram(diff_tmp,'BinLimits',hdi_tmp,'BinWidth',0.03,'EdgeColor',edge_cols{k},'DisplayStyle','stairs','LineWidth',0.8)
    hold on

    xlim([-0.7, 1.0])
    ylim([0, 2000])
end

% Vertical line at zero
xline(0,'Color','#A2142F');

title('No: Diff(bin_{i+1} - bin_{i})')


%%% HDIs for no-responses %%%
subplot(5,3,i)
i = i + 1;

% Loop time interval contrasts
for k = 1:3       

    diff_tmp = mratio.fit(k+1).mcmc.samples.mu_logMratio_rS2(:,:,1)-mratio.fit(k).mcmc.samples.mu_logMratio_rS2(:,:,1);
    hdi_tmp = calc_HDI(diff_tmp(:));    

    histogram(diff_tmp,'BinLimits',hdi_tmp,'BinWidth',0.03,'EdgeColor',edge_cols{k},'DisplayStyle','stairs','LineWidth',0.8)
    hold on

    xlim([-0.7, 1.0])
    ylim([0, 2000])

end

% Vertical line at zero
xline(0,'Color','#A2142F');

title('Yes: Diff(bin_{i+1} - bin_{i})')


legend({'0-200ms vs. 200-400ms','200-400ms vs. 400-600ms','400-600ms vs. 600-800ms',' '})


%% Display range of MCMC samples (min/max)
for k = 1:4; disp(min(min(exp(mratio.fit(k).mcmc.samples.mu_logMratio_rS1(:,:,1))))); end
for k = 1:4; disp(min(min(exp(mratio.fit(k).mcmc.samples.mu_logMratio_rS2(:,:,1))))); end
for k = 1:4; disp(max(max(exp(mratio.fit(k).mcmc.samples.mu_logMratio_rS1(:,:,1))))); end
for k = 1:4; disp(max(max(exp(mratio.fit(k).mcmc.samples.mu_logMratio_rS2(:,:,1))))); end
