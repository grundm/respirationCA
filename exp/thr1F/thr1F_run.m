function data = thr1F_run(thr1F,aio_s,ID)
% thr1F_run(thr1F,aio_s,ID) starts the threshold estimation for a 
% one-alternative forced choice (1AFC) paradigm (yes/no task).
%   
% It relies on the Palamedes toolbox by Prins & Kingdom (2009), see
% http://www.palamedestoolbox.org
%
% % Input variables %
%   thr1F            - settings structure (doc thr1F_setup)
%   aio_s           - daq acquisition session object
%   ID              - participant ID as string (e.g., p_data.ID)
%
% % Output variables %
%   data            - output structure (e.g., sequence, responses, ...)
%
% Author:           Martin Grund
% Last update:      December 18, 2018

try

%% SETUP %%    
    
%% Input device
lpt = dio_setup(thr1F.lpt_adr1,thr1F.lpt_adr2,thr1F.lpt_dir);


%% Screen

% Open window
window = Screen('OpenWindow',0,thr1F.window_color);
HideCursor;

% Get screen frame rate
Priority(1); % recommended by Mario Kleiner's Perfomance & Timing How-to
flip_t = Screen('GetFlipInterval',window,200);
Priority(0);
data.flip_t = flip_t;

% Set font
Screen('TextFont',window,thr1F.txt_font);

% Compute response text location
[data.window(1),data.window(2)] = Screen('WindowSize',window);

% left top right bottom
resp_rect1 = [data.window(1)*0.5 - thr1F.txt_size,...
              data.window(2)*0.5 - thr1F.resp_offset - thr1F.txt_size,...
              data.window(1)*0.5 + thr1F.txt_size,...
              data.window(2)*0.5 - thr1F.resp_offset];
resp_rect2 = [data.window(1)*0.5 - thr1F.txt_size,...
              data.window(2)*0.5 + thr1F.resp_offset + thr1F.txt_size,...
              data.window(1)*0.5 + thr1F.txt_size,...
              data.window(2)*0.5 + thr1F.resp_offset];


%% Random number generator
data.rng_state = set_rng_state(thr1F);


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
        data.resp_txt = thr1F.resp_txt;
    case {3,4} % NJ
        data.resp_txt = flipud(thr1F.resp_txt);
end


%% Up/down method

UD = PAL_AMUD_setupUD('up', thr1F.UD_up, ...
                      'down', thr1F.UD_down, ...
                      'stepSizeUp', thr1F.UD_stepSizeUp, ...
                      'stepSizeDown', thr1F.UD_stepSizeDown, ...
                      'stopCriterion', thr1F.UD_stopCriterion, ...
                      'stopRule', thr1F.UD_stopRule, ...
                      'startValue', thr1F.UD_startValue, ...
                      'xMax', thr1F.UD_xMax, ...
                      'xMin', thr1F.UD_xMin);
                  

%% Sequence

% Compute number of null trials based on specified rate
data.trials_psi_null = round(thr1F.trials_psi * (thr1F.psi_null_rate/(1-thr1F.psi_null_rate)) );

% Shuffle sequence of psi method target and non-target (null) trials
data.seq_psi = Shuffle([ones(thr1F.trials_psi,1); zeros(data.trials_psi_null,1);]);

% Test block
data.seq_test = [];

for i = 1:length(thr1F.trials_test)
	data.seq_test = [data.seq_test; thr1F.stim_test(i)*ones(thr1F.trials_test(i),1)];
end

data.seq_test = Shuffle(data.seq_test);

% Create complete sequence block vector
% 1 - UD; 2 - PM; 3 - Test
data.block = [ones(UD.stopRule,1); 2*ones(length(data.seq_psi),1); 3*ones(length(data.seq_test),1)];
data.seq = [-1*ones(UD.stopRule,1); data.seq_psi; data.seq_test;];
 
% Random stimlus delays
% data.stim_delay = round( thr1F.stim_delay(1) +
data.stim_delay = Shuffle(round_dec(linspace(thr1F.stim_delay(1),thr1F.stim_delay(2),size(data.seq,1)),4))';

% Random fixation cross presentaion (fix_t)
data.fix_t = Shuffle(round_dec(linspace(thr1F.fix_t(1),thr1F.fix_t(2),size(data.seq,1)),4))';
% Make fix_t divisible by flip time
data.fix_t = data.fix_t - mod(data.fix_t, flip_t);


%% Stimulus

% Generate default waveform vector
[data.stim_wave, data.stim_offset] = rectpulse2(thr1F.pulse_t,1,aio_s.Rate,thr1F.pre_pulse_t,thr1F.wave_t);

% Generate TTL pulse waveform vector
[data.TTL_wave, data.TTL_offset] = rectpulse2(thr1F.TTL_t,thr1F.TTL_V,aio_s.Rate,thr1F.pre_pulse_t,thr1F.wave_t);



%% Timing

% Screen flips
[data.onset_fix,...
 data.onset_cue,...
 data.onset_resp,...
 data.onset_ITI] = deal(cell(size(data.seq,1),5));

% Analog input data
data.ai_data = cell(size(data.seq,1),1);

% Stimulus trigger
[data.ao_trigger_pre,...
 data.ao_trigger_post,...
 data.ao_error] = deal(zeros(size(data.seq,1),1));

% Response tracking
[data.resp_btn,...
data.resp_t,...
data.resp] = deal(zeros(size(data.seq,1),1));


%% EXPERIMENTAL PROCEDURE %%
% Priority(1) seems to cause DataMissed events for analog input

%% Instructions

% Get directory of instruction images for response-button mapping condition
instr_dir = [fileparts(mfilename('fullpath')) thr1F.instr_dir];
instr_subdir = dir([instr_dir thr1F.instr_subdir_wildcard num2str(data.resp_btn_map) '*']);

% Load image data
instr_images = load_images([instr_dir instr_subdir.name '/'],thr1F.instr_img_wildcard);

% Show images
[data.btn_instr, data.onset_instr] = show_instr_img(instr_images,window,lpt);

% Delete image data from memory
clear instr_images img_texture


%% Threshold assesment

for i = 1:size(data.seq,1);
    
    % Set font size for symbols
    Screen('TextSize',window,thr1F.cue_size);

    %%% FIX %%%
    
    DrawFormattedText(window,thr1F.fix,'center','center',thr1F.txt_color);

    [data.onset_fix{i,:}] = Screen('Flip',window);
    
    
    %%% CUE %%%
    
    DrawFormattedText(window,thr1F.cue,'center','center',thr1F.cue_color);
    [data.onset_cue{i,:}] = Screen('Flip',window,data.onset_fix{i,1} + data.fix_t(i,1) - flip_t);   
    
    
    %%% STIMULUS %%%
    
    % Select intensity
    
    if data.seq(i) == 0;                     % Null trial
      
        data.intensity(i,1) = 0;
                        
    elseif data.seq(i) == -1;                % Up/down method trial
        
        data.intensity(i,1) = UD.xCurrent;
            
    elseif (data.seq(i) > 0) && ~PM.stop     % Psi method trial
        
        data.intensity(i,1) = PM.xCurrent;
        
    elseif (data.seq(i) > 0) && PM.stop      % Test trials        
                
        data.intensity(i,1) = round_dec(PAL_Quick(data.PF_params_PM,data.seq(i),'Inverse'),thr1F.stim_dec);        
        
        % Check if test intensity would exceed stimulus range
        if abs(data.intensity(i,1)) > max(thr1F.stim_range)
            data.intensity(i,1) = max(thr1F.stim_range);            
        end                
    end
    
    % Buffer waveform
    stop(aio_s);
    
    queueOutputData(aio_s,[data.stim_wave*data.intensity(i,1) data.TTL_wave]);       
                
    % Stimulus delay
    data.ao_trigger_pre(i,1) = WaitSecs('UntilTime',data.onset_cue{i,1}+data.stim_delay(i)-(thr1F.pre_pulse_t/1000));
       
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
    
    %%% RESPONSE %%%
    
    % Response options    
    Screen('TextSize',window,thr1F.txt_size);
    DrawFormattedText(window,data.resp_txt(1,:),'center','center',thr1F.txt_color,[],[],[],[],[],resp_rect1);
    DrawFormattedText(window,data.resp_txt(2,:),'center','center',thr1F.txt_color,[],[],[],[],[],resp_rect2);
    
    [data.onset_resp{i,:}] = Screen('Flip',window,data.onset_cue{i,1} + thr1F.cue_t - flip_t);
        
    % Wait for key press
    [data.resp_btn(i,1),data.resp_t(i,1),data.resp_port(i,:)] = parallel_button(thr1F.maxRespTime,data.onset_resp{i,1},thr1F.resp_window,thr1F.debounceDelay,lpt);
           
    % Evaluate response
    switch data.resp_btn(i,1)
        case num2cell(thr1F.btn)
            switch data.resp_txt(thr1F.btn==data.resp_btn(i,1))
                case thr1F.resp_txt(1)
                    data.resp(i,1) = 1; % yes
                case thr1F.resp_txt(2)
                    data.resp(i,1) = 0; % no
            end
%         case thr1F.btn_esc
%             break
        otherwise
            data.resp(i,1) = 0;
    end   
    
    % Update UD or PM structure with response (1: correct, 0: incorrect)
    if data.seq(i) == -1;
        
        UD = PAL_AMUD_updateUD(UD, data.resp(i));
        
        if UD.stop
            
            %% Analyze up/down method %%
            
            % Compute mean from up/down method
            data.UD_mean = PAL_AMUD_analyzeUD(UD,thr1F.UD_stopCriterion,thr1F.UD_meanNumber);
            
            % Compute intensity range of trials used for mean
            data.UD_range = [min(UD.x(end-thr1F.UD_meanNumber-1:end)) max(UD.x(end-thr1F.UD_meanNumber-1:end))];            
            
            % Define prior alpha range ± thr1F.UD_range_factor of up/down mean range
            tmp_priorAlphaRange = data.UD_range(1)*thr1F.UD_range_factor(1):thr1F.priorAlphaSteps:data.UD_range(2)*thr1F.UD_range_factor(2);

            % Define priors for psychometric function fitting
            grid.alpha = tmp_priorAlphaRange;
            grid.beta = thr1F.priorBetaRange;
            grid.gamma = thr1F.priorGammaRange;
            grid.lambda = thr1F.priorLambdaRange;
            
            % Compute intensity - response frequency matrix of up/down
            % method
            data.x_resp_freq_UD = count_resp([UD.x' UD.response']);

            % Fit psychometric function to up/down method data
            [data.PF_params_UD, data.posterior_UD] = PAL_PFBA_Fit(data.x_resp_freq_UD(:,1),data.x_resp_freq_UD(:,2),data.x_resp_freq_UD(:,3),grid,@PAL_Quick);
                       
            % ??? MAYBE USE THE THRESHOLD ESTIMATION            
            
            %% Prepare psi method %%
            
            % Shift prior alpha range
            tmp_priorAlphaRange = tmp_priorAlphaRange + (data.PF_params_UD(1,1)-mean(tmp_priorAlphaRange));
            
            % Restrict stimulus range to prior alpha
%             tmp_stim_range = thr1F.stim_range(thr1F.stim_range>=min(tmp_priorAlphaRange));
%             tmp_stim_range = tmp_stim_range(tmp_stim_range<=max(tmp_priorAlphaRange));
            tmp_stim_range = thr1F.stim_range;
            
            % The current warning can be ignored, because the posterior 
            % output by PAL_PFBA_Fit is in the format "alpha x beta"
            
            PM = PAL_AMPM_setupPM('numtrials', thr1F.trials_psi, ...
                      'stimRange', tmp_stim_range, ...
                      'PF', thr1F.PF, ...
                      'priorAlphaRange', tmp_priorAlphaRange, ...
                      'priorBetaRange', thr1F.priorBetaRange, ...
                      'priorGammaRange', thr1F.priorGammaRange, ...
                      'priorLambdaRange', thr1F.priorLambdaRange, ...
                      'prior', data.posterior_UD);            
        end
    
    elseif (data.seq(i) > 0) && ~PM.stop
        PM = PAL_AMPM_updatePM(PM, data.resp(i));
        if PM.stop
            data.near = PM.threshold(end);
            data.PF_params_PM = [PM.threshold(end); 10.^PM.slope(end); PM.guess(end); PM.lapse(end);];
            data.supra = PAL_Quick(data.PF_params_PM, thr1F.supra_p, 'Inverse');
        end
    end

end


%% End procedures

% Save up/down method and psi method structures
data.UD = UD;
data.PM = PM;

% Close window
Screen('Close',window);
ShowCursor;


%% Error handling

catch lasterr
    sca;
    rethrow(lasterr);
end

end