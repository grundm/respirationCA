function apply_pulse_1F(aio_s,psi)
% apply_pulse_1F(aio_s,psi) allows to familiarize participants with
% electrical finger nerve stimulation in advance of automatic threshold
% assessment (psi1AFC_run). Experimentator can control via commond window 
% input requests at which finger which intensity shall be applied. 
% Stimulation waveform is defined by psi1AFC_setup_*.
%
% Input:
%   aio_s           - daq acquisition session object
%   psi           - settings structure (doc psi1AFC_setup_*)
%
% Author:           Martin Grund
% Last update:      December 13, 2018

%% Prepare output vectors

% Generate default waveform vector
stim_wave = rectpulse2(psi.pulse_t,1,aio_s.Rate,psi.pre_pulse_t,psi.wave_t);

% Generate TTL pulse waveform vector
TTL_wave = rectpulse2(psi.TTL_t,psi.TTL_V,aio_s.Rate,psi.pre_pulse_t,psi.wave_t);


%% Stimulation loop

while 1
   
    intensity = input('Stimulus intensity (mA): ');
    
    if isempty(intensity)
        break
    end
    
    % Overscripe test intensity with maximum intensity when they exceed stimulus range
    if intensity > max(psi.stim_range)
        disp(['Input (' num2str(intensity) ') above maximum intensity (' num2str(max(psi.stim_range)) ' mA).'])
        intensity = max(psi.stim_range);
        disp(['-> Applied intensity: ' num2str(intensity) ' mA'])
    end
    
    % Zero negative intensities
    if intensity < 0    
        disp(['Input (' num2str(intensity) ') below zero.'])
        intensity = 0;
        disp(['-> Applied intensity: ' num2str(intensity) ' mA'])
    end
    
    stop(aio_s);
    
    queueOutputData(aio_s,[stim_wave*intensity TTL_wave]);       

    startForeground(aio_s);

    stop(aio_s);
    
    WaitSecs(0.5);
end