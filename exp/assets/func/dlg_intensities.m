function [near,supra] = dlg_intensities(PF_params,near_p,supra_p,stim_steps,stim_max,dec_digits)
% [near,supra] = dlg_intensities(PF_params,near_p,supra_p,stim_steps,stim_max,dec_digits)
% creates an input dialog with a near- and supra-threshold intensity. The
% default values are based on the psychometric function parameters (PF_params)
% and the performance levels (detection probability: near_p & supra_p) in
% the input arguments.
%
% It checks if the intensities are below the allowed intensity maximum
% (stim_max) while considering also the stimulus steps (stim_steps).
%
% If the input passes the checks it is displayed again for confirmation. If
% the default input was changed the detection probabilities are updated.
% 
%
% Input:
%   PF_params       - Psychometric function parameters matrix
%                     [threshold; slope; guess_rate; lapse_rate;]
%   near_p          - performance for near-threshold intensity (0-1)
%   supra_p         - performance for supra-threshold intensity (0-1)
%   stim_steps      - multiplicators for intensities (see nt.stim_steps)
%   stim_max        - allowed intensity maximum
%   dec_digits      - number of decimal digits values are rounded to
%
% Author:           Martin Grund
% Last update:      November 16, 2015


%% Compute intensities

PF = @(x) round(PAL_Quick(PF_params, x)*100);

near = round_dec(PAL_Quick(PF_params,near_p,'Inverse'),dec_digits);
supra = round_dec(PAL_Quick(PF_params,supra_p,'Inverse'),dec_digits);


%% Dialog & evaluation
while 1
    
    %%% Intensity input dialog %%%
    
    dlg.in = inputdlg({['Near-threshold (~' num2str(PF(near)) '%) in mA:'],...
                       ['Supra-threshold (~' num2str(PF(supra)) '%) in mA:']},...
                        'Stimulus intensities',...
                        1,...
                        {num2str(near),num2str(supra)});

    near = round_dec(str2num(dlg.in{1}),dec_digits);
    supra = round_dec(str2num(dlg.in{2}),dec_digits);

    
    %%% Evaluate input %%%

    if ~isempty(near) && ~isempty(supra)

        % Check if single numbers
        if numel(near) > 1 || numel(supra) > 1
            e = errordlg(['We only accept single numbers (e.g., 3.50) with ' num2str(dec_digits) ' decimal digits maximum.']);
            uiwait(e);
            
        % Check if below intensity maximum
        elseif round_dec(abs(near)*max(stim_steps),dec_digits) > stim_max || round_dec(abs(supra)*max(stim_steps),dec_digits) > stim_max            
            e = errordlg(sprintf(['At least one of the intensities ' ...
                                  'is above the maximum (' num2str(stim_max) ' mA).\n\n' ...
                                  'Near-threshold (~' num2str(PF(near)) '%%): ' num2str(near) ' mA\n\n' ...
                                  'Supra-threshold (~' num2str(PF(supra)) '%%): ' num2str(supra) ' mA']));                                  
            uiwait(e);
            
        % Display entered intensities
        else
            choice = questdlg(sprintf(['Continue with these intensities?\n\n' ...
                                       'Near-threshold (~' num2str(PF(near)) '%%): ' num2str(near) ' mA\n\n' ...
                                       'Supra-threshold (~' num2str(PF(supra)) '%%): ' num2str(supra) ' mA']),...
                                       'Intensity check');
                                   
            if strcmp(choice,'No') || strcmp(choice,'Cancel')
                e = errordlg('You did not like your chosen intensities.');
                uiwait(e);
            else
                WaitSecs(1);
                break
            end
        end

    else
        e = errordlg(['We only accept single numbers (e.g., 3.50) with ' num2str(dec_digits) ' decimal digits maximum.']);
        uiwait(e);
    end
end