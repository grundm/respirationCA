function thr1F_save(p_data,thr1F_data,thr1F,file_name_end)
% thr1F_save(p_data,thr1F_data,thr1F,file_name_end) saves the output of the 
% threshold assesment with thr1F_run (thr1F_data, incl. Palamedes up/down 
% method and psi method structures), as well as the settings for thr1F_run 
% (thr1F).
%
% Additionally, it creates a table with the single trial data in each line.
%
% % Input variables %
%   p_data          - output of participant_data
%   thr1F           - output of thr1F_setup (setting structure)
%   thr1F_data      - output of thr1F_run (threshold assesment data)
%   file_name_end   - string that defines end of filename
%
% Author:           Martin Grund
% Last update:      January 18, 2019

% Setup data logging
file_name = [thr1F.file_prefix p_data.ID];

% Create participant data directory
if ~exist(p_data.dir,'dir');
    mkdir('.',p_data.dir);
end

% Save Matlab variables
mat_file_tmp = [p_data.dir file_name '_data_' file_name_end '.mat'];

if exist(mat_file_tmp, 'file')
    disp('MAT-file of experiment exists. Generated random file name to prevent overwritting.')
    save([p_data.dir file_name '_data_' file_name_end '_' num2str(round(sum(100*clock))) '.mat'],'p_data','thr1F_data','thr1F');
else
    save([p_data.dir file_name '_data_' file_name_end '.mat'],'p_data','thr1F_data','thr1F');
end


%% Save trial data

% Make thr1F_data easily accessbile
d = thr1F_data;

% Get current date
date_str = datestr(now,'yyyy/mm/dd');

% Open file
data_file = fopen([p_data.dir file_name '_trials_' file_name_end '.txt'],'a');

% Write header
fprintf(data_file,'ID\tage\tgender\tdate\tblock\ttrial\tstim_type\tintensity\tresp\tresp_t\tresp_btn\tstim_delay\tfix_t\tao_error\n');

for i = 1:length(d.seq)
   fprintf(data_file,'%s\t%s\t%s\t%s\t%.0f\t%.0f\t%.0f\t%.6f\t%.0f\t%.4f\t%.0f\t%.4f\t%.6f\t%.0f\n',p_data.ID,p_data.age,p_data.gender,date_str,d.block(i),i,d.seq(i,1),d.intensity(i),d.resp(i),d.resp_t(i),d.resp_btn(i),d.stim_delay(i),d.fix_t(i),d.ao_error(i,1));
end

fclose(data_file);