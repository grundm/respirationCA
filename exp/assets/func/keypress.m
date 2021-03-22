function [key_code,respTime] = keypress(maxRespTime,startTime,timeWindow)
%   Function monitors key presses for a specified time frame (maxRespTime),
%   since the set start time and returns the key code and response time.
%
%   In-/output times are milliseconds.
%
%   Optional: Set timeWindow = 'fixed' if you want no break after first key
%   pres.
%
%   Author:           Martin Grund
%   Last update:      July 22, 2015

if nargin < 3
    timeWindow = 'variable';
end

keyHit = 0;
maxRespTime = maxRespTime/1000; % ms -> s

while GetSecs-startTime <= maxRespTime
    
    [keyIsDown, secs, keyCode] = KbCheck;
%     disp(GetSecs-startTime)
    if keyIsDown && keyHit == 0
        respTime = (GetSecs-startTime)*1000;
        
        key_code = find(keyCode==1,1,'first');
        
        keyHit = 1;
        
        if strcmp(timeWindow,'variable')
            break
        end
    end    
end

if keyHit == 0;
    key_code = 0;
    respTime = maxRespTime*1000;
end