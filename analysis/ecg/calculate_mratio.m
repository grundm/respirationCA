% This script calculate hierachichal meta d prime (type 2 SDT) according to
% Fleming et al. 2019
%there is no implementation of the response conditional version in R
%therefore this matlab script (I prefer R)

%IMPORTANT:

%make sure that confidence is coded as 1,2 and not 0,1 


%
% written by Carina Forster, 2021
%

clear all
close all


addpath('/data/hu_forster/Documents/HMeta-d-master/Matlab/');
addpath('/data/hu_forster/Documents/');


%%

%load and prepare data
data = readtable('/data/hu_forster/Documents/trials_cc_bin.csv');


%how many subjects ?
%Nsub = max(data(:,1)); %subject 1001? 

Nsub = 40;

%VariableNames = {'isyes', 'sayyes', 'conf', 'ID', 'RT1', 'RT2', 'block', 'subblock', 'intensity'};

%delete unnecessary columns

data(:,1) = [];
data(:,1) = [];
data(:,3) = [];
data(:,3) = [];


%change confidence from 0 to 1

data.conf(data.resp2==1)=2;
data.conf(data.resp2==0)=1;

%convert to array


columnnames = data(1,:);

data = table2array(data);

%% create vector for Mratio
%loop over all subjects to get confidence vector for meta d prime



%vector = zeros(40,4);

for i = 1:Nsub
      
        idx=ismember(data(:,1),i)
        
        d =data(idx,:)
         
        %input: stimulus, response,confidence rating
        
        [nR_S1(i,:), nR_S2(i,:)] = trials2counts(d(:,5), d(:,8), d(:,36),2);
         
end

%%
%S1 = catch trials, no response
%S2 = signal trials, yes response

%convert to cell array for meta d prime function

nR_S1 = mat2cell(nR_S1, ones(1, size(nR_S1, 1)), size(nR_S1, 2))';
nR_S2 = mat2cell(nR_S2, ones(1, size(nR_S2, 1)), size(nR_S2, 2))';

%% fit M ratio

%Get default parameters

mcmc_params = fit_meta_d_params;

% Change defaults to make response-conditional (recommended for detection
% taks)

mcmc_params.response_conditional = 1;

% Fit group data all at once

fit = fit_meta_d_mcmc_group(nR_S1, nR_S2, mcmc_params);

%DIC Martins data: 877,72 RC (response conditional)
%without rc: 928,25

%better fit for response conditional model
%%
%those scripts calculate meta d per subject using MLE or SSE 

% for i = 1:Nsub
%  
%     fit = fit_rs_meta_d_MLE(nR_S1(i),nR_S2(i));
% 
%     %fit = type2_SDT_SSE(nR_S1, nR_S2);
%     
% end

%% Plot output

%lower M ratio for noise trials

%catch trials = S1

plotSamples(exp(fit.mcmc.samples.mu_logMratio_rS1))

%signal trials = S2

plotSamples(exp(fit.mcmc.samples.mu_logMratio_rS2))

%% for Martins data
%calculate Meta d for different phases in cardiac cycle
%which phase in cardiac cycle?

idx_first = ismember(data(:,35),1);
idx_second = ismember(data(:,35),2);
idx_third = ismember(data(:,35),3);
idx_fourth = ismember(data(:,35),4);

data_first = data(idx_first,:);
data_second = data(idx_second,:);
data_third = data(idx_third,:);
data_fourth = data(idx_fourth,:);


%calculate M ratio for each phase 

nR_S1_first = zeros(40,4);
nR_S2_first = zeros(40,4);

for i = 1:Nsub
%         
        idx=ismember(data_first(:,1),i);
        
        d =data_first(idx,:);
        
        [nR_S1_first(i,:), nR_S2_first(i,:)] = trials2counts(d(:,5), d(:,8), d(:,36),2);
        
end


for i = 1:Nsub
%         
        idx=ismember(data_second(:,1),i);
        
        d =data_second(idx,:);
        
        [nR_S1_second(i,:), nR_S2_second(i,:)] = trials2counts(d(:,5), d(:,8), d(:,36),2);
        
end



for i = 1:Nsub
%         
        idx=ismember(data_third(:,1),i);
        
        d =data_third(idx,:);
        
        [nR_S1_third(i,:), nR_S2_third(i,:)] = trials2counts(d(:,5), d(:,8), d(:,36),2);
        
end



for i = 1:Nsub
%         
        idx=ismember(data_fourth(:,1),i);
        
        d =data_fourth(idx,:);
        
        [nR_S1_fourth(i,:), nR_S2_fourth(i,:)] = trials2counts(d(:,5), d(:,8), d(:,36),2);
        
end


%S1 = catch trials, no response
%S2 = signal trials, yes response

%convert to cell array for meta d prime function

nR_S1_first = mat2cell(nR_S1_first, ones(1, size(nR_S1_first, 1)), size(nR_S1_first, 2))';
nR_S2_first = mat2cell(nR_S2_first, ones(1, size(nR_S2_first, 1)), size(nR_S2_first, 2))';

nR_S1_second = mat2cell(nR_S1_second, ones(1, size(nR_S1_second, 1)), size(nR_S1_second, 2))';
nR_S2_second = mat2cell(nR_S2_second, ones(1, size(nR_S2_second, 1)), size(nR_S2_second, 2))';

nR_S1_third = mat2cell(nR_S1_third, ones(1, size(nR_S1_third, 1)), size(nR_S1_third, 2))';
nR_S2_third = mat2cell(nR_S2_third, ones(1, size(nR_S2_third, 1)), size(nR_S2_third, 2))';

nR_S1_fourth = mat2cell(nR_S1_fourth, ones(1, size(nR_S1_fourth, 1)), size(nR_S1_fourth, 2))';
nR_S2_fourth = mat2cell(nR_S2_fourth, ones(1, size(nR_S2_fourth, 1)), size(nR_S2_fourth, 2))';

%%
%fit M ratio

%Get default parameters

mcmc_params = fit_meta_d_params;

% Change defaults to make response-conditional
mcmc_params.response_conditional = 0;

% Fit group data all at once for each phase

fit_first = fit_meta_d_mcmc_group(nR_S1_first, nR_S2_first, mcmc_params);
fit_second = fit_meta_d_mcmc_group(nR_S1_second, nR_S2_second, mcmc_params);
fit_third = fit_meta_d_mcmc_group(nR_S1_third, nR_S2_third, mcmc_params);
fit_fourth = fit_meta_d_mcmc_group(nR_S1_fourth, nR_S2_fourth, mcmc_params);

%%
%save outputs

%Mratio rc

Mratios_first = [fit_first.Mratio_rS1, fit_first.Mratio_rS2];
Mratios_second = [fit_second.Mratio_rS1, fit_second.Mratio_rS2];
Mratios_third = [fit_third.Mratio_rS1, fit_third.Mratio_rS2];
Mratios_fourth = [fit_fourth.Mratio_rS1, fit_fourth.Mratio_rS2];

writematrix(Mratios_first,'Mratios_firstint.csv'); 
writematrix(Mratios_second,'Mratios_secondint.csv');
writematrix(Mratios_third,'Mratios_thirdint.csv');
writematrix(Mratios_fourth,'Mratios_fourthint.csv');

%Mratio non rc

Mratios_first = [fit_first.Mratio];
Mratios_second = [fit_second.Mratio];
Mratios_third = [fit_third.Mratio];
Mratios_fourth = [fit_fourth.Mratio];

writematrix(Mratios_first,'Mratiosnonrc_firstint.csv'); 
writematrix(Mratios_second,'Mratiosnonrc_secondint.csv');
writematrix(Mratios_third,'Mratiosnonrc_thirdint.csv');
writematrix(Mratios_fourth,'Mratiosnonrc_fourthint.csv');

%Meta d prime

Mdprime_first = [fit_first.meta_d];
Mdprime_second = [fit_second.meta_d];
Mdprime_third = [fit_third.meta_d];
Mdprime_fourth = [fit_fourth.meta_d];

writematrix(Mdprime_first,'Md_firstint.csv'); 
writematrix(Mdprime_second,'Md_secondint.csv');
writematrix(Mdprime_third,'Md_thirdint.csv');
writematrix(Mdprime_fourth,'Md_fourthint.csv');

%criteria

c_first = [fit_first.t2ca_rS1, fit_first.t2ca_rS2];
c_second = [fit_second.t2ca_rS1, fit_second.t2ca_rS2];
c_third = [fit_third.t2ca_rS1, fit_third.t2ca_rS2];
c_fourth = [fit_fourth.t2ca_rS1, fit_fourth.t2ca_rS2];



writematrix(c_first,'c_firstint.csv'); 
writematrix(c_second,'c_secondint.csv');
writematrix(c_third,'c_thirdint.csv');
writematrix(c_fourth,'c_fourthint.csv');


%I prefer R for data visualization
%%
%Plot M-ratio for the 4 heart phases and compare

figure
boxplot(fit_fourth.Mratio_rS1)


plotSamples(fit_second.mcmc.samples.Mratio_rS1(:,:,1));
plotSamples(fit_third.mcmc.samples.mu_logMratio(:,:,1))
plotSamples(fit_fourth.mcmc.samples.mu_logMratio(:,:,1))

[R,P] = corrcoef(fit_first.Mratio_rS1, fit_second.Mratio_rS1)

% Compute HDI of difference between tasks 
sampleDiff = fit_first.mcmc.samples.mu_logMratio(:,:,1) - fit_fourth.mcmc.samples.mu_logMratio(:,:,1);
hdi = calc_HDI(sampleDiff(:));
fprintf(['\n HDI on difference in log(meta-d''/d''): ', num2str(hdi) '\n\n'])

% % Plot difference in meta-d/d ratio between two tasks 
plotSamples(sampleDiff)
% 
% Plot estimate of correlation
h1 = figure;
set(gcf, 'Position', [200 200 400 300])
h= histogram('Normalization', 'probability');
xlabel('\rho');
ylabel('Posterior density');
%line([rho rho],[0 max(h.Values)+0.015], 'LineWidth', 2, 'Color', 'k', 'LineStyle', '--')
ci = calc_CI(fit_group_no.mcmc.samples.rho(:));
line([ci(1) ci(2)],[0.002 0.002], 'LineWidth', 3, 'Color', [1 1 1])
box off
set(gca, 'FontSize', 14, 'XLim', [-1 1])