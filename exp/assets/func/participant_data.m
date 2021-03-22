function p_data = participant_data(data_path)
% p_data = participant_data(data_path) opens a dialog that allows to enter the 
% participant's ID, age and gender. It also creates a data directory under the specified 
% data_path and ID (p_data.dir = [data_path p_data.ID '/'];), if no such data directory 
% exists.
%
% It returns the structure p_data with the cells "ID", "age", "gender" and "dir".
%
% Do not use "data_" as directory name: "exist('data_01','dir')" or 
% "exist('data','dir')" return everywhere "7", indicating an existing directory even 
% though it does not exist.
%
%
% Author:           Martin Grund
% Last update:      December 17, 2018

%% Settings
dialogTitle = 'Participant data';
prompt = {'Participant ID:' 
          'Age:'
          'Gender:'
          };
defaultAnswer = {'01','99','x'};
numLines = 1;

%% Dialog
answer = inputdlg(prompt,dialogTitle,numLines,defaultAnswer);

p_data.ID = answer{1};
p_data.age = answer{2};
p_data.gender = answer{3};

%% Check if data directory exists
if exist([data_path p_data.ID],'dir')
    e = errordlg(sprintf(['Data directory exists already. This indicates that ID ' p_data.ID ' was used.']));                                  
    uiwait(e);
    
    % Call input dialog again
    p_data = participant_data(data_path);
        
else

    % Create data directory
    p_data.dir = [data_path p_data.ID '/'];
    mkdir('.',p_data.dir);

end