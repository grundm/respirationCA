function [oxi_data_c,res_c,acq_num_c] = clean_data(oxi_data)
% function [oxi_data_c,res_c,acq_num_c] = clean_data(oxi_data) - hard coded
% data separation for blocks with more than 150 expected triggers. Returns
% "cleaned" data, triggers (res_c) and overview triggers per block
% (acq_num).
%
% Author:           Martin Grund
% Last update:      August 25, 2020


%% Check data based on triggers

[res,acq_num] = get_trigger(oxi_data);


%% Participants with blocks unequal 150 trials

IDs2check = acq_num(any(acq_num(:,3:6) ~= 150,2),:)

% 1           4         242         392         542         692
% 2           4         150         300         300         150
% 3           4         250         400         550         700
% 4           4         150         300         450         600
% 5           4         150         323         473         681         
% 1001        4         150         300         450         600
% 19          4         150         150         150         300        
% 22          4         147         150         150         150
% 39          4         336         150         150         150

% ID22 -> missed first 3 trials in block #1

%% Clean ACQ data

oxi_data_c = oxi_data;

t_padding = 4000; % 4 seconds at 1000 Hz


%% Cut recording before and after last trigger if 150 triggers

% Loop participants
for i = 1:length(oxi_data_c)
    
    ID_ind_res = find(strcmp({res.ID},oxi_data_c(i).ID));
    
    % Loop blocks
    for j = 1:length(oxi_data_c(i).acq)
        
        % Check if 150 triggers
        if acq_num(acq_num(:,1)==str2double(oxi_data_c(i).ID),2+j) == 150
            
            oxi_data_c(i).acq(j).data = oxi_data_c(i).acq(j).data(res(ID_ind_res).trigger(j).time(1)-t_padding:res(ID_ind_res).trigger(j).time(end)+t_padding,:);        
            
        end
    end
    
end

%% ID01
% - 92 additional triggers before block #1 in 1 recording of block #1
% - Recordings of subsequent blocks include all previous blocks.
% -> Block #1: res(1).trigger(1).time(93:end)

ID_ind_data = find(strcmp({oxi_data_c.ID},'01'));
ID_ind_res = find(strcmp({res.ID},'01'));

% One line for each block .acq(1:4) & .trigger(1:4) 
% -> cut previous recordings
oxi_data_c(ID_ind_data).acq(1).data = oxi_data_c(ID_ind_data).acq(1).data(res(ID_ind_res).trigger(1).time(93)-t_padding:res(ID_ind_res).trigger(1).time(end)+t_padding,:);
oxi_data_c(ID_ind_data).acq(2).data = oxi_data_c(ID_ind_data).acq(2).data(res(ID_ind_res).trigger(2).time(93+150)-t_padding:res(ID_ind_res).trigger(2).time(end)+t_padding,:);
oxi_data_c(ID_ind_data).acq(3).data = oxi_data_c(ID_ind_data).acq(3).data(res(ID_ind_res).trigger(3).time(93+150*2)-t_padding:res(ID_ind_res).trigger(3).time(end)+t_padding,:);
oxi_data_c(ID_ind_data).acq(4).data = oxi_data_c(ID_ind_data).acq(4).data(res(ID_ind_res).trigger(4).time(93+150*3)-t_padding:res(ID_ind_res).trigger(4).time(end)+t_padding,:);


%% ID02
% - Recording of block #2 includes block #1.
% - Recording of block #3 includes first block #4 and then block #3.
% - Based on the MATLAB log file: "successful" blocks: #1, #2, #3, #4, #3
% - First "successful" run of block #3 stored seperately, because no ExG
% recording.
% - Last run of block #3 considered for behavioral analysis.
% -> Delete block #1 in recording of block #2.
% -> Delete block #4 in recording of block #3.

ID_ind_data = find(strcmp({oxi_data_c.ID},'02'));
ID_ind_res = find(strcmp({res.ID},'02'));

% Check if data for first 10 triggers in recording of block #1 and #2 are
% indentical:
figure('Name','ID02 - Block #1 vs. #2');
plot(oxi_data_c(ID_ind_data).acq(1).data(res(ID_ind_res).trigger(1).time(1):res(ID_ind_res).trigger(1).time(10),2),'k-')
hold on
plot(oxi_data_c(ID_ind_data).acq(2).data(res(ID_ind_res).trigger(2).time(1):res(ID_ind_res).trigger(2).time(10),2)+1,'b-')

% Check if data for first 10 triggers in recording of block #4 and #3 are
% indentical:
figure('Name','ID02 - Block #4 vs. #3');
plot(oxi_data_c(ID_ind_data).acq(4).data(res(ID_ind_res).trigger(4).time(1):res(ID_ind_res).trigger(4).time(10),2),'k-')
hold on
plot(oxi_data_c(ID_ind_data).acq(3).data(res(ID_ind_res).trigger(3).time(1):res(ID_ind_res).trigger(3).time(10),2)+1,'b-')

% Cut first block #1 from recording of block #2
oxi_data_c(ID_ind_data).acq(2).data = oxi_data_c(ID_ind_data).acq(2).data(res(ID_ind_res).trigger(2).time(151)-t_padding:res(ID_ind_res).trigger(2).time(end)+t_padding,:);

% Cut first block #4 from recording of block #3
oxi_data_c(ID_ind_data).acq(3).data = oxi_data_c(ID_ind_data).acq(3).data(res(ID_ind_res).trigger(3).time(151)-t_padding:res(ID_ind_res).trigger(3).time(end)+t_padding,:);


%% ID03 
% - Recording of block #1 100 triggers before block #1.
% - Recordings of subsequent blocks include all previous blocks.
% -> res(3).trigger(4).time(101:end)

ID_ind_data = find(strcmp({oxi_data_c.ID},'03'));
ID_ind_res = find(strcmp({res.ID},'03'));

% One line for each block .acq(1:4) & .trigger(1:4) 
% -> cut previous recordings
oxi_data_c(ID_ind_data).acq(1).data = oxi_data_c(ID_ind_data).acq(1).data(res(ID_ind_res).trigger(1).time(101)-t_padding:res(ID_ind_res).trigger(1).time(end)+t_padding,:);
oxi_data_c(ID_ind_data).acq(2).data = oxi_data_c(ID_ind_data).acq(2).data(res(ID_ind_res).trigger(2).time(101+150)-t_padding:res(ID_ind_res).trigger(2).time(end)+t_padding,:);
oxi_data_c(ID_ind_data).acq(3).data = oxi_data_c(ID_ind_data).acq(3).data(res(ID_ind_res).trigger(3).time(101+150*2)-t_padding:res(ID_ind_res).trigger(3).time(end)+t_padding,:);
oxi_data_c(ID_ind_data).acq(4).data = oxi_data_c(ID_ind_data).acq(4).data(res(ID_ind_res).trigger(4).time(101+150*3)-t_padding:res(ID_ind_res).trigger(4).time(end)+t_padding,:);

%% ID04
% - Recordings of subsequent blocks include all previous blocks.

ID_ind_data = find(strcmp({oxi_data_c.ID},'04'));
ID_ind_res = find(strcmp({res.ID},'04'));

% One line for each block .acq(2:4) & .trigger(2:4) 
% -> cut previous recordings
oxi_data_c(ID_ind_data).acq(2).data = oxi_data_c(ID_ind_data).acq(2).data(res(ID_ind_res).trigger(2).time(1+150)-t_padding:res(ID_ind_res).trigger(2).time(end)+t_padding,:);
oxi_data_c(ID_ind_data).acq(3).data = oxi_data_c(ID_ind_data).acq(3).data(res(ID_ind_res).trigger(3).time(1+150*2)-t_padding:res(ID_ind_res).trigger(3).time(end)+t_padding,:);
oxi_data_c(ID_ind_data).acq(4).data = oxi_data_c(ID_ind_data).acq(4).data(res(ID_ind_res).trigger(4).time(1+150*3)-t_padding:res(ID_ind_res).trigger(4).time(end)+t_padding,:);


%% ID05
% - Recordings of subsequent blocks include all previous blocks.
% - Recording of block #2 includes 23 additional triggers at the end of block #1.
% - Recording of block #4 includes 58 additional triggers at the end of block #3.
% -> Block #2: res(5).trigger(2).time(174:end)
% -> Block #3: res(5).trigger(3).time(324:end)
% -> Block #4: res(5).trigger(4).time(532:end)

ID_ind_data = find(strcmp({oxi_data_c.ID},'05'));
ID_ind_res = find(strcmp({res.ID},'05'));

% One line for each block .acq(2:4) & .trigger(2:4) 
% -> cut previous recordings
oxi_data_c(ID_ind_data).acq(2).data = oxi_data_c(ID_ind_data).acq(2).data(res(ID_ind_res).trigger(2).time(24+150)-t_padding:res(ID_ind_res).trigger(2).time(end)+t_padding,:);
oxi_data_c(ID_ind_data).acq(3).data = oxi_data_c(ID_ind_data).acq(3).data(res(ID_ind_res).trigger(3).time(24+150*2)-t_padding:res(ID_ind_res).trigger(3).time(end)+t_padding,:);
oxi_data_c(ID_ind_data).acq(4).data = oxi_data_c(ID_ind_data).acq(4).data(res(ID_ind_res).trigger(4).time(24+58+150*3)-t_padding:res(ID_ind_res).trigger(4).time(end)+t_padding,:);

%% ID1001
% - Recordings of subsequent blocks include all previous blocks.

ID_ind_data = find(strcmp({oxi_data_c.ID},'1001'));
ID_ind_res = find(strcmp({res.ID},'1001'));

% One line for each block .acq(2:4) & .trigger(2:4) 
% -> cut previous recordings
oxi_data_c(ID_ind_data).acq(2).data = oxi_data_c(ID_ind_data).acq(2).data(res(ID_ind_res).trigger(2).time(1+150)-t_padding:res(ID_ind_res).trigger(2).time(end)+t_padding,:);
oxi_data_c(ID_ind_data).acq(3).data = oxi_data_c(ID_ind_data).acq(3).data(res(ID_ind_res).trigger(3).time(1+150*2)-t_padding:res(ID_ind_res).trigger(3).time(end)+t_padding,:);
oxi_data_c(ID_ind_data).acq(4).data = oxi_data_c(ID_ind_data).acq(4).data(res(ID_ind_res).trigger(4).time(1+150*3)-t_padding:res(ID_ind_res).trigger(4).time(end)+t_padding,:);


%% ID19 
% Recording of block #4 includes block #3.

ID_ind_data = find(strcmp({oxi_data_c.ID},'19'));
ID_ind_res = find(strcmp({res.ID},'19'));

% Cut first block #3 from recording of block #4
oxi_data_c(ID_ind_data).acq(4).data = oxi_data_c(ID_ind_data).acq(4).data(res(ID_ind_res).trigger(4).time(151)-t_padding:res(ID_ind_res).trigger(4).time(end)+t_padding,:);


%% ID39 
% - Recording of block #1 entails additonal 186 triggers before block #1.
% -> Block #1: res(40).trigger(1).time(187:end)

ID_ind_data = find(strcmp({oxi_data_c.ID},'39'));
ID_ind_res = find(strcmp({res.ID},'39'));

% Cut first block #3 from recording of block #4
oxi_data_c(ID_ind_data).acq(1).data = oxi_data_c(ID_ind_data).acq(1).data(res(ID_ind_res).trigger(1).time(187)-t_padding:res(ID_ind_res).trigger(1).time(end)+t_padding,:);


%% Check data based on triggers

[res_c,acq_num_c] = get_trigger(oxi_data_c);

IDs2check_c = acq_num_c(any(acq_num_c(:,3:6) ~= 150,2),:)

% 22          4         147         150         150         150
