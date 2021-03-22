function [trigger_t,trigger_date,start_wait,stop_wait] = wait_for_scanner(lpt,TR,trigger_bit,trigger_max)
% [trigger_t,trigger_date,start_wait,stop_wait] = wait_for_scanner(lpt,TR,trigger_bit,trigger_max)
% waits for a specified number of triggers send by the scanner.
%
% For this, it enables interrupt requests (IRQ) by setting the control port
% bit 5 on 1, and disables IRQs at the end again. For details see:
% http://retired.beyondlogic.org/spp/parallel.htm#5
%
% Input:
%   lpt             - structure (output by dio_setup)
%   TR              - MRI repetition time in s
%   trigger_bit     - status port bit reflecting hardware pin with MRI trigger input
%   trigger_max     - number of triggers to wait
%
% Author:           Martin Grund
% Last update:      November 11, 2018

%% Subfunctions
get_status = @() lpt.get(lpt.dio,lpt.status_adr);

%% Setup
start_wait = GetSecs;

trigger_count = 0;
trigger_t = zeros(trigger_max,1);
trigger_date = trigger_t;

% Enable interrupt request (IRQ)
%io32(dio,ctrl_adr,bitset(get_lpt(dio,ctrl_adr),5,1));
lpt.set_IRQ(lpt,1);

%% Monitor

while 1
    
    if bitget(get_status(),trigger_bit) == 0
        
        trigger_count = trigger_count + 1;
        
        trigger_t(trigger_count) = GetSecs;
        
        trigger_date(trigger_count) = now;
        
        if trigger_count == trigger_max
           break 
        end
        
        WaitSecs('UntilTime',trigger_t(trigger_count)+TR/2);
    end
end

%% Close
% Disable interrupt request (IRQ)
lpt.set_IRQ(lpt,0);

stop_wait = GetSecs;