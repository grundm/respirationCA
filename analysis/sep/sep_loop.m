%% Analysis of somatsensory evoked potential (SEP) data
%
% Author:           Martin Grund
% Last Update:      March 24, 2020


%% Settings

code_dir = '/Users/martin/ownCloud/promotion/experiment/respirationCA/code/respirationca';

data_dir = '/Users/martin/ownCloud/promotion/experiment/respirationCA/data';


%% Set paths

cd(code_dir)

% Add paths
addpath([code_dir '/sep']);
addpath([code_dir '/assets/eeglab14_1_1b']);
addpath([code_dir '/assets/hline_vline']);
%addpath(genpath([code_dir '/assets']));

eeglab


%% Define data

exclude_IDs = 0; % 16?

%IDs = [1001, 1:40];

IDs = [1001, 1:17];

for i = 1:length(IDs)
    
    if IDs(i) < 10
        ID = ['0' num2str(IDs(i))];
    else
        ID = num2str(IDs(i));
    end    

    if ~any(exclude_IDs == IDs(i))
        
        ID
        result(i) = sep_analysis(ID,data_dir);
    
    end        
end

%%

% ID1001, ID02 -> 
% ID1001, ID01-05 -> wrong filter for SEP recording (< 250 Hz)

figure; 

for i = 1:length(result)
%for i = 7:length(result)
    
    subplot(4,5,i)
    %subplot(3,4,i-6)
    hold on;
    plot(result(i).EEG_times,result(i).ERP_near_miss, 'r')
    plot(result(i).EEG_times,result(i).ERP_near_hit, 'k')
    xlim([-20 60])
    title(result(i).ID)

end

%% 

for i = 7:length(result)
    
    ERP_mean.near(i-6,:) = result(i).ERP_near;
    ERP_mean.null(i-6,:) = result(i).ERP_null;
    ERP_mean.near_hit(i-6,:) = result(i).ERP_near_hit;
    ERP_mean.near_miss(i-6,:) = result(i).ERP_near_miss;
    
end

%% Mean across all participants

for i = 1:length(result)
    
    ERP_mean.near(i,:) = result(i).ERP_near;
    ERP_mean.null(i,:) = result(i).ERP_null;
    ERP_mean.near_hit(i,:) = result(i).ERP_near_hit;
    ERP_mean.near_miss(i,:) = result(i).ERP_near_miss;
    
end

%%

figure; hold on;

plot(result(10).EEG_times,mean(ERP_mean.null,1), 'k')
plot(result(10).EEG_times,mean(ERP_mean.near_hit,1), 'g')
plot(result(10).EEG_times,mean(ERP_mean.near_miss,1), 'r')

plot(result(10).EEG_times,mean(ERP_mean.near_hit-ERP_mean.near_miss,1), 'g-')
hline(0)

%xlim([-20 60])
%title(result(i).ID)

