function data = run_exp(s,aio_s,ID,PF_params,seq)
% run_exp(s,ao,ai,ID,PF_params,seq) starts the experiment.
% 
% % Input variables %
%   s               - settings structure (doc setup)
%   aio_s           - daq acquisition session object
%   ID              - participant ID as string (e.g., p_data.ID)
%   PF_params       - psychometric function parameters (e.g., psi_data.PF_params_PM)
% 
% % Output variables %
%   data            - data output structure
% 
% Author:           Martin Grund
% Last update:      December 18, 2018


%% Intensity input dialog %%
data.near = dlg_1intensity(PF_params,s.near_d,s.near_name,s.stim_steps,s.stim_max,s.stim_dec);


%% SETUP %%    
try


%% Random number generator
data.rng_state = set_rng_state(s);


%% Response-button mapping

% respirationCA: 2*2 responses -> 4 conditions (2^2)
% 1: JN & SU
% 2: JN & US
% 3: NJ & SU
% 4: NJ & US

data.resp_btn_map_num = 4;

data.resp_btn_map = mod(str2double(ID)-1,data.resp_btn_map_num)+1;

% 1st response
switch data.resp_btn_map
    case {1,2} % JN
        data.resp1_txt = s.resp1_txt;
    case {3,4} % NJ
        data.resp1_txt = flipud(s.resp1_txt);
end

% 2nd response
switch data.resp_btn_map
    case {1,3} % SU
        data.resp2_txt = s.resp2_txt;
    case {2,4} % US
        data.resp2_txt = flipud(s.resp2_txt);
end

%% Sequence (doc seq)

data.seq = seq;


%% Stimulus

% Generate default waveform vector
[data.stim_wave, data.stim_offset] = rectpulse2(s.pulse_t,1,aio_s.Rate,s.pre_pulse_t,s.wave_t);

% Generate TTL pulse waveform vector
[data.TTL_wave, data.TTL_offset] = rectpulse2(s.TTL_t,s.TTL_V,aio_s.Rate,s.pre_pulse_t,s.wave_t);


%% Timing

data.trial_t = s.fix_t + s.cue_t + s.resp1_max_t + s.resp2_max_t;

% Trial screen flips
[data.onset_fix,...
 data.onset_cue,...
 data.onset_resp1,...
 data.onset_resp1_p,...
 data.onset_resp2] = deal(cell(size(data.seq,1),5));

%  data.onset_stim_p,...
%  data.onset_ITI

% % End of wait for trial end
% data.wait_trial_end = zeros(size(data.seq,1),1);

% Before block screen flips
[data.btn_instr,...
data.onset_instr] = deal(cell(1,5));

% Pause screen flips
% data.onset_pause = cell(s.blocks-1,5);

% Analog input data
data.ai_data = cell(size(data.seq,1),1);

% Stimulus trigger
[data.ao_trigger_pre,...
 data.ao_trigger_post,...
 data.ao_error] = deal(zeros(size(data.seq,1),1));

% Response tracking
[data.resp1_btn,...
data.resp1_t,...
data.resp1,...
data.resp2_btn,...
data.resp2_t,...
data.resp2] = deal(zeros(size(data.seq,1),1));

%% Parallel port
lpt = dio_setup(s.lpt_adr1,s.lpt_adr2,s.lpt_dir);


%% Screen

window = Screen('OpenWindow',0,s.window_color);
HideCursor;

Screen('TextFont',window,s.txt_font);               
             
% Get screen frame rate
Priority(1); % recommended by Mario Kleiner's Perfomance & Timing How-to
flip_t = Screen('GetFlipInterval',window,s.get_flip_i);
Priority(0);
data.flip_t = flip_t;

% Update sequence - make fix_t divisible by flip time
data.seq(:,5) = data.seq(:,5) - mod(data.seq(:,5), flip_t);

% Compute response text location
[data.window(1),data.window(2)] = Screen('WindowSize',window);

% left top right bottom
resp_rect1 = [data.window(1)*0.5 - s.txt_size,...
              data.window(2)*0.5 - s.resp1_offset - s.txt_size,...
              data.window(1)*0.5 + s.txt_size,...
              data.window(2)*0.5 - s.resp1_offset];
resp_rect2 = [data.window(1)*0.5 - s.txt_size,...
              data.window(2)*0.5 + s.resp1_offset + s.txt_size,...
              data.window(1)*0.5 + s.txt_size,...
              data.window(2)*0.5 + s.resp1_offset];

    
%% EXPERIMENTAL PROCEDURE %%
% Priority(1) seems to cause DataMissed events for analog input


%% Instructions

% Get directory (based on response-button mapping)
instr_dir = [fileparts(mfilename('fullpath')) s.instr_dir];
instr_subdir = dir([instr_dir s.instr_subdir_wildcard num2str(data.resp_btn_map) '*']);

% Load image data
instr_images = load_images([instr_dir instr_subdir.name '/'],s.instr_img_wildcard);

% Show images
[data.btn_instr, data.onset_instr] = show_instr_img(instr_images,window,lpt);

% Delete image data from memory
clear instr_images img_texture


%% Trial loop

for i = 1:size(data.seq,1)    
    
    %%% FIX %%%
    
    % Set font size for symbols
    Screen('TextSize',window,s.cue_size);
    DrawFormattedText(window,s.fix,'center','center',s.txt_color);
    
    [data.onset_fix{i,:}] = Screen('Flip',window);
    
    
    %%% CUE %%%
    
    DrawFormattedText(window,s.cue,'center','center',s.cue_color);
    
    [data.onset_cue{i,:}] = Screen('Flip',window,data.onset_fix{i,1} + data.seq(i,5) - flip_t);
    
    
    %%% STIMULUS %%%
    
    % Select intensity
    
    switch data.seq(i,2)
        case 0 % null
            data.intensity(i,1) = 0;
        case 1 % near
            data.intensity(i,1) = round_dec(data.near*data.seq(i,3),s.stim_dec);
        %case 2 % supra
        %    data.intensity(i,1) = round_dec(data.supra*data.seq(i,3),s.stim_dec);
    end
    
    % Buffer waveform
    stop(aio_s);
    
    queueOutputData(aio_s,[data.stim_wave*data.intensity(i,1) data.TTL_wave]);       
    
    % Stimulus delay (locked to scanner trigger)
    data.ao_trigger_pre(i,1) = WaitSecs('UntilTime',data.onset_cue{i,1} + data.seq(i,4) - (s.pre_pulse_t/1000));
    
    % Start analog output (triggers waveform immediately)
    try
        [ai_rec.data,ai_rec.time] = startForeground(aio_s);
        stop(aio_s);
    catch lasterr
        disp(['Trial ', num2str(i), ': ', lasterr.message]);
        data.ao_error(i,1) = 1;
        stop(aio_s);
    end
    
    data.ao_trigger_post(i,1) = GetSecs;
    data.ai_data{i,1} = ai_rec;
    
    
    %%% RESPONSE 1 - DETECTION %%%
    
    % Response options
    Screen('TextSize',window,s.txt_size);
    DrawFormattedText(window,data.resp1_txt(1,:),'center','center',s.txt_color,[],[],[],[],[],resp_rect1);
    DrawFormattedText(window,data.resp1_txt(2,:),'center','center',s.txt_color,[],[],[],[],[],resp_rect2);

    [data.onset_resp1{i,:}] = Screen('Flip',window,data.onset_cue{i,1} + s.cue_t - flip_t);
        
    % Wait for key press
    [data.resp1_btn(i,1),data.resp1_t(i,1),data.resp1_port(i,:)] = parallel_button(s.resp1_max_t,data.onset_resp1{i,1},s.resp_window,s.debounceDelay,lpt);        

    
    % RT dependent fix between responses
%     if s.resp1_max_t-data.resp1_t(i,1) > s.resp_p_min_t
        Screen('TextSize',window,s.cue_size);
        DrawFormattedText(window,s.fix,'center','center',s.txt_color);
        [data.onset_resp1_p{i,:}] = Screen('Flip',window);
%     end        
    
    
    %%% RESPONSE 2 - CONFIDENCE %%%
    
    % Response options
    Screen('TextSize',window,s.txt_size);
    DrawFormattedText(window,data.resp2_txt(1,:),'center','center',s.txt_color,[],[],[],[],[],resp_rect1);
    DrawFormattedText(window,data.resp2_txt(2,:),'center','center',s.txt_color,[],[],[],[],[],resp_rect2);
    
    [data.onset_resp2{i,:}] = Screen('Flip',window,data.onset_resp1_p{i,1} + s.resp_p_min_t - flip_t);
        
    % Wait for key press
    [data.resp2_btn(i,1),data.resp2_t(i,1),data.resp2_port(i,:)] = parallel_button(s.resp2_max_t,data.onset_resp2{i,1},s.resp_window,s.debounceDelay,lpt);   
    
    
    %%% RESPONSE EVALUATION %%%
    
    % Response 1
    switch data.resp1_btn(i,1)
        case num2cell(s.btn_resp1)
            switch data.resp1_txt(s.btn_resp1==data.resp1_btn(i,1))
                case s.resp1_txt(1)
                    data.resp1(i,1) = 1; % yes
                case s.resp1_txt(2)
                    data.resp1(i,1) = 0; % no
            end
        case s.btn_esc
            break
        otherwise
            data.resp1(i,1) = 0;
    end             

    % Response 2
    switch data.resp2_btn(i,1)
        case num2cell(s.btn_resp2)
            switch data.resp2_txt(s.btn_resp2==data.resp2_btn(i,1))
                case s.resp2_txt(1)
                    data.resp2(i,1) = 1; % confident
                case s.resp2_txt(2)
                    data.resp2(i,1) = 0; % unconfident
            end
        case s.btn_esc
            break
        otherwise
            data.resp2(i,1) = 0;
    end
    
    
    %%% WAIT UNTIL TRIAL END %%%

    if i == size(data.seq,1)
        [data.onset_fix{i+1,:}] = Screen('Flip',window);
        data.wait_block_end = WaitSecs('UntilTime',data.onset_resp2{i,1} + s.resp2_max_t);
    end
    
end


%% End procedures

% Close all screens
sca;


%% Error handling

catch lasterr
    sca;
    rethrow(lasterr);
end