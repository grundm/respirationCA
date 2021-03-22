function [btn,onset] = show_instr_img(instr_images,window,lpt)
% [btn,onset] = show_instr_img(instr_images,window,lpt) displays 
% sequentially all instruction images in the structure instr_images. It
% returns the flip onset times (onset) and the pressed buttons (btn.n) as
% well as the response time (btn.t).
%
% It requires an open window.
%
% Author:           Martin Grund
% Last update:      November 20, 2017

% Prepare output variables
onset = cell(length(instr_images),5);
btn.n = zeros(length(instr_images),1);
btn.t = btn.n;

for i = 1:length(instr_images)
    img_texture = Screen('MakeTexture',window,instr_images{i});
    Screen('DrawTexture',window,img_texture);
    [onset{i,:}] = Screen('Flip',window);    
    Screen('Close',img_texture);
    
    % Wait for button press
    [btn.n(i),btn.t(i),btn.port(i,:)] = parallel_button(inf,onset{i,1},'variable',0.025,lpt);    
end

% Blank screen after last button press
[onset{i,:}] = Screen('Flip',window);