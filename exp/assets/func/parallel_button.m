function [button,respTime,port] = parallel_button(maxRespTime,startTime,timeWindow,debounceDelay,lpt)
% [button,respTime,port] = parallel_button(maxRespTime,startTime,timeWindow,lpt)
% monitors the parallel data port for a specified time frame (maxRespTime) 
% since the specified start time (startTime).
%
% It returns the data port pin that changed from 1 to 0, response time in
% s and the status of the parallel data, status, control and extended
% control register port. For special cases the variable "button" is assigned:
% 0  - no button press
% -1 - nonstop button press starting before with parallel part monitoring
% 9  - multiple button presses
%
% Optional: Set timeWindow = 'fixed' if you want no break after the first 
% button press and want to consider only the last pressed button.
%
% Input:
%   maxRespTime     - response window in seconds (e.g., 1.500)
%   startTime       - real system time in seconds (e.g., GetSecs)
%   timeWindow      - 'fixed' or 'variable' (see above for details)
%   lpt             - structure (output by dio_setup)
%   debounceDelay   - necessary button press duration in seconds
%
% Author:           Martin Grund
% Last update:      November 20, 2017

%% Prepare inline functions

bibit_on = @() lpt.set_bibit(lpt,1);
byte_mode_on = @() lpt.set_mode(lpt,lpt.mode_byte);

get_data = @() lpt.get(lpt.dio,lpt.data_adr);

%% Set byte mode and bidirectional bit
byte_mode_on();
bibit_on();

button = 0;

%% Loop until no buttons are pressed
data_port = get_data();

while data_port ~= 255 && (GetSecs-startTime <= maxRespTime)
    % Set bidirectional bit
    bibit_on();
    
    data_port = get_data();
    button = -1;
end

%% Wait for button press

last_data_port = data_port;
lastDebounceTime = 0;

while GetSecs-startTime <= maxRespTime

    % Set bidirectional bit
    bibit_on();
    
    data_port = get_data();
    
    % Debouncing code idea by https://www.arduino.cc/en/Tutorial/Debounce#toc5
    % Change in data_port_state?
    if data_port ~= last_data_port
       % Reset debounce time
       lastDebounceTime = GetSecs; 
    end
    
    if GetSecs-lastDebounceTime > debounceDelay
    
        if data_port ~= 255
            % all data pins low suggests unidirectional mode
            if data_port ~= 0            

                respTime = lastDebounceTime-startTime;

                % Get zero data bit
                button = find(bitget(data_port,1:8)==0);

                if strcmp(timeWindow,'variable') && length(button) == 1
                    break
                end
            end
        end
    
    end    
        
    last_data_port = data_port;
end

% Save parallel ports: data - status - control - ecr
port = [data_port lpt.get(lpt.dio,lpt.status_adr) lpt.get(lpt.dio,lpt.ctrl_adr) lpt.get(lpt.dio,lpt.ecr)];

if length(button) > 1
    disp(datestr(now));
    disp(['button = ' num2str(button)]);
    disp(['respTime = ' num2str(respTime)]);
    disp(['port = ' num2str(port)]);
    
    button = 9;
end

% If no button (0) was pressed or signal with request onset that stayed nonstop or no button press followed (-1)
% assign maximum response time
if button == 0 || button == -1
    respTime = maxRespTime;
    disp(datestr(now));
    disp(['button = ' num2str(button)]);
    disp(['port = ' num2str(port)]);
end