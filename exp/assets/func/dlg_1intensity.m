function intensity = dlg_1intensity(PF_params,intensity_d,intensity_name,stim_steps,stim_max,dec_digits)
% intensity = dlg_1intensity(PF_params,intensity_d,intensity_name,stim_steps,stim_max,dec_digits)
% creates an input dialog with an intensity in mA. The default value is 
% based on the psychometric function parameters (PF_params) and the 
% detection rate in the input arguments.
%
% It checks if the intensity is below the allowed intensity maximum
% (stim_max) while considering also the stimulus steps (stim_steps).
%
% If the input passes the checks it is displayed again for confirmation. If
% the default input was changed the detection probability is updated.
% 
%
% Input:
%   PF_params       - Psychometric function parameters matrix
%                     [threshold; slope; guess_rate; lapse_rate;]
%   intensity_d     - detection rate for intensity (0-1)
%   intensity_name  - name for intensity in dialog (e.g. "Near-threshold")
%   stim_steps      - multiplicators for intensity (see nt.stim_steps)
%   stim_max        - allowed intensity maximum
%   dec_digits      - number of decimal digits values are rounded to
%
% Author:           Martin Grund
% Last update:      October 20, 2016

%% Compute intensities

PF = @(x) round(PAL_Quick(PF_params, x)*100);

intensity = round_dec(PAL_Quick(PF_params,intensity_d,'Inverse'),dec_digits);


%% Dialog & evaluation
while 1
    
    %%% Intensity input dialog %%%
    
    dlg.in = inputdlg([intensity_name ' (~' num2str(PF(intensity)) '%) in mA:'],...
                      'Stimulus intensity',...
                      1,...
                      {num2str(intensity)});

    intensity = round_dec(str2num(dlg.in{1}),dec_digits);

    
    %%% Evaluate input %%%

    if ~isempty(intensity)

        % Check if single numbers
        if numel(intensity) > 1
            e = errordlg(['We only accept single numbers (e.g., 3.50) with maximum ' num2str(dec_digits) ' decimal digits.']);
            uiwait(e);
            
        % Check if below intensity maximum
        elseif round_dec(abs(intensity)*max(stim_steps),dec_digits) > stim_max            
            e = errordlg(sprintf(['The intensity is above the maximum (' num2str(stim_max) ' mA).\n\n' ...
                                  intensity_name ' (~' num2str(PF(intensity)) '%%): ' num2str(intensity) ' mA']));                                  
            uiwait(e);
            
        % Display entered intensities
        else
            choice = questdlg(sprintf(['Continue with this intensity?\n\n' ...
                                       intensity_name ' (~' num2str(PF(intensity)) '%%): ' num2str(intensity) ' mA']),...
                                       'Intensity check');
                                   
            if strcmp(choice,'No') || strcmp(choice,'Cancel')
                e = errordlg('You did not like your chosen intensity.');
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