% calculate_bstar_imaging_parameters.m
% Written by Nam Gyun Lee
% Email: namgyunl@usc.edu, ggang56@gmail.com (preferred)
% Started: 08/23/2022, Last modified: 06/15/2023

%% Gradient mode
switch grad_mode
    case 'Fast'
        max_grad = 24;      % Max gradient strength [mT/m]
        max_slew = 180.18;  % Maximum slew rate [mT/m/ms]
    case 'Normal'
        max_grad = 22;      % Max gradient strength [mT/m]
        max_slew = 100;     % Maximum slew rate [mT/m/ms]
    case 'Whisper'
        max_grad = 22;      % Max gradient strength [mT/m]
        max_slew = 50;      % Maximum slew rate [mT/m/ms]
end

%% Set system limits
sys = mr.opts('MaxGrad'       , max_grad, 'GradUnit', 'mT/m' , ...
              'MaxSlew'       , max_slew, 'SlewUnit', 'T/m/s', ...
              'rfRingdownTime', 20e-6 , ...
              'rfDeadTime'    , 100e-6, ...
              'adcDeadTime'   , 10e-6 , ...
              'B0'            , B0);

%% Calculate the real dwell time [sec]
%--------------------------------------------------------------------------
% real dwell time [sec]
% IDEA p219: dRealDwellTime denotes the dwell time with oversampling
%--------------------------------------------------------------------------
% round-down dwell time to 100 ns (sys.adcRasterTime  = 100 ns)
real_dwell_time = round((1 / bandwidth) / (readout_os_factor * base_resolution) * 1e7) * 1e-7;

real_dwell_time1 = real_dwell_time;
real_dwell_time2 = real_dwell_time;

bandwidth1 = 1 / (real_dwell_time1 * (readout_os_factor * base_resolution));
bandwidth2 = 1 / (real_dwell_time2 * (readout_os_factor * base_resolution));

%% Calculate the duration of an ADC event [sec]
adc_samples = nr_echoes * base_resolution * readout_os_factor; % true number of samples
adc_duration = adc_samples * real_dwell_time;

%% Calculate the gradient amplitude based on FOV [mT/m]
amplitude = 1 / (fov_read * readout_os_factor) / (sys.gamma * real_dwell_time * 1e-3); % [mT/m]

%% Calculate the ramp time of a trapezoid (RampTime) [sec]
ramp_time = ceil((sys.maxGrad / sys.maxSlew) / sys.gradRasterTime) * sys.gradRasterTime; % [sec]

%% Calculate the slew rate [mT/m/ms]
slew_rate = amplitude / ramp_time * 1e-3; % [mT/m] / [sec] * [sec/1e3msec] => *1e-3 [mT/m/msec]

%% Calculate the duration of a trapezoid (TotalTime) [sec]
total_time = ceil(adc_duration / (sys.gradRasterTime * 2)) * (sys.gradRasterTime * 2) / 2; % [sec]

%% Calculate spatial resolution [m]
%--------------------------------------------------------------------------
% Calculate the area of a gradient lobe [mT/m*sec]
%--------------------------------------------------------------------------
gradient_area = (total_time - ramp_time) * amplitude; % [mT/m*sec]

%--------------------------------------------------------------------------
% Calculate the maximum k-space value [cycle/m]
% [Hz/T] * [mT/m*sec] * [T/1e3mT] => [cycle/m]
%--------------------------------------------------------------------------
krmax = sys.gamma * gradient_area * 1e-3;

%--------------------------------------------------------------------------
% Calculate the spatial resolution [m]
%--------------------------------------------------------------------------
resolution = 1 / (2 * krmax);

%% Calculate TE1 [sec]
TE1 = rf_length / 2 + sys.rfRingdownTime + sys.adcDeadTime;

%% Calculate TE2 [sec]
TE2 = TE1 + total_time * 2;

%% Calculate TR [sec]
TR = sys.rfDeadTime + rf_length / 2 + TE2 + sys.adcDeadTime;

%% Calculate "TR delay" [sec]
TR_delay = TR - (rf_length / 2 + TE2);

%% Calculate the slice thickness [m]
slice_thickness = fov_read / slices_per_slab;

%% Calculate the number of readouts per interleaf
nr_readouts = floor(nr_radial_views / nr_interleaves); % number of readouts per interleaf
nr_radial_views = nr_readouts * nr_interleaves;

%% Display parameters
fprintf('----------------------- Contrast (Common) ------------------------\n');
fprintf('TR         = %4.2f [ms]\n', TR * 1e3);
fprintf('TE1        = %4.2f [ms]\n', TE1 * 1e3);
fprintf('TE2        = %4.2f [ms]\n', TE2 * 1e3);
fprintf('Flip angle = %4d [deg]\n', flip_angle);
fprintf('------------------------------------------------------------------\n');

fprintf('----------------------- Resolution (Common) ----------------------\n');
fprintf('FoV read        = %5.0f [mm]\n', fov_read * 1e3);
fprintf('FoV phase       = %5.1f [%%]\n', 100.0);
fprintf('Slice thickness = %5.2f [mm]\n', slice_thickness * 1e3);
fprintf('Base resolution = %5d\n', base_resolution);
fprintf('Radial views    = %5d\n', nr_radial_views);
fprintf('------------------------------------------------------------------\n');

fprintf('----------------------- Sequence (Part 1) ------------------------\n');
fprintf('Contrasts       = %7d [mm]\n', nr_echoes);
fprintf('Bandwidth 1     = %7.0f [Hz/Px]\n', bandwidth1);
fprintf('Bandwidth 2     = %7.0f [Hz/Px]\n', bandwidth2);
fprintf('Readout mode    = %s\n', readout_mode);
fprintf('------------------------------------------------------------------\n');

fprintf('----------------------- Sequence (Special) -----------------------\n');
fprintf('Trajectory       = %11s      RF length  = %4.0f [us]\n', trajectory, rf_length * 1e6);
fprintf('Readout mode     = %11s      RF dummies = %4.0f\n', 'Combined', rf_dummies);
fprintf('F-mod            = %11.1f      RF phase   = %4d [deg]\n', f_mod, rf_phase);
fprintf('Tilt             = %11.1f\n', tilt);
fprintf('Interleaves      = %11d\n', nr_interleaves);
fprintf('Gradient factor  = %11.2f      TR delay   = %4.0f [us]\n', gradient_factor, (TR - (rf_length / 2 + TE2)) * 1e6);
fprintf('Amplitude factor = %11.2f\n', amplitude_factor);
fprintf('Spoiler factor   = %11.2f\n', spoiler_factor);
fprintf('Slew rate        = %11.2f [mT/m/ms]\n', slew_rate);
fprintf('Amplitude        = %11.2f [mT/m]\n', amplitude);
fprintf('------------------------------------------------------------------\n');
fprintf('RUTime/TotalTime/ADC[0]/Resolution: %3.0f/%3.0f/%6.1f/%4.2f\n', ramp_time * 1e6, total_time * 1e6, adc_duration * 1e6, resolution * 1e3);
fprintf('dwell time (echo 1) = %4.0f [nsec]\n', real_dwell_time1 * 1e9);
fprintf('dwell time (echo 2) = %4.0f [nsec]\n', real_dwell_time2 * 1e9);
