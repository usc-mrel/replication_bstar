% calculate_bstar_imaging_parameters.m
% Written by Nam Gyun Lee
% Email: namgyunl@usc.edu, ggang56@gmail.com (preferred)
% Started: 08/23/2022, Last modified: 04/18/2024

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
    case 'FreeMax'
        max_grad = 26;      % Max gradient strength [mT/m]
        max_slew = 45;      % Maximum slew rate [mT/m/ms]
end

%% Set system limits
max_grad_derated = max_grad / sqrt(3);
max_slew_derated = max_slew / sqrt(3);

sys_derated = mr.opts('MaxGrad'       , max_grad_derated, 'GradUnit', 'mT/m' , ...
                      'MaxSlew'       , max_slew_derated, 'SlewUnit', 'T/m/s', ...
                      'rfRingdownTime', 20e-6 , ...
                      'rfDeadTime'    , 100e-6, ...
                      'adcDeadTime'   , 10e-6 , ...
                      'B0'            , B0);

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
real_dwell_time1 = real_dwell_time;
real_dwell_time2 = real_dwell_time;

%% Calculate the gradient amplitude based on FOV [mT/m]
amplitude = 1 / (fov_read * readout_os_factor) / (sys.gamma * real_dwell_time * 1e-3); % [mT/m]

%% Calculate the ramp time of a trapezoid (RampTime) [sec]
ramp_time = ceil((sys.maxGrad / sys.maxSlew) / sys.gradRasterTime) * sys.gradRasterTime; % [sec]

%% Calculate the slew rate [mT/m/ms]
slew_rate = amplitude / ramp_time * 1e-3; % [mT/m] / [sec] * [sec/1e3msec] => *1e-3 [mT/m/msec]

%% Calculate the maximum k-space value [rad/m]
if true_resolution_flag
    grid_scale_factor = 0.8;
    resolution_scale_factor = sqrt(pi / 4);
else
    grid_scale_factor = 1;
    resolution_scale_factor = 1;
end

krmax = (2 * pi) / (2 * resolution * resolution_scale_factor); % [rad/m]

%% Calculate the area of a gradient lobe [mT/m*sec]
%--------------------------------------------------------------------------
% [rad/m] * [cycle/2pi rad] / [Hz/T] * [1e3mT/T] => [mT/m*sec]
%--------------------------------------------------------------------------
gradient_area = krmax / (2 * pi * sys.gamma) * 1e3;

%% Calculate the total time of a gradient lobe [sec]
%--------------------------------------------------------------------------
%                              _______________
%                             /|             |\
%                            / |             | \
%                           /  |             |  \
%                     _____/   |             |   \______
%                          |<->|             |<->| ramp_time
%                          |<------------------->| total_time
%
% gradient_area = (total_time - ramp_time) * amplitude
%--------------------------------------------------------------------------
total_time = round(gradient_area / amplitude / (sys.gradRasterTime * 2)) * (sys.gradRasterTime * 2) + ramp_time;

%% Calculate the number of ADC samples per echo
adc_samples_per_echo = round(total_time / real_dwell_time);

%% Calculate the base_resolution
base_resolution = adc_samples_per_echo / readout_os_factor;

%% Calculate the duration of an ADC event [sec]
adc_samples = nr_echoes * base_resolution * readout_os_factor; % true number of samples
adc_duration = adc_samples * real_dwell_time;

%% Calculate "delay_pre" and "delay_post" [sec]
adc_duration_pre  = real_dwell_time * discard_pre;
adc_duration_post = real_dwell_time * discard_post;

delay_pre  = ceil((sys.adcDeadTime + adc_duration_pre)  / sys.gradRasterTime) * sys.gradRasterTime;
delay_post = ceil((sys.adcDeadTime + adc_duration_post) / sys.gradRasterTime) * sys.gradRasterTime;

%% Calculate TE1 [sec]
TE1 = rf_length / 2 + sys.rfRingdownTime + delay_pre;

%% Calculate TE2 [sec]
TE2 = TE1 + total_time * 2;

%% Calculate TR [sec]
TR = sys.rfDeadTime + rf_length / 2 + TE2 + delay_post;

%% Calculate "TR_delay" [sec]
TR_delay = TR - (rf_length / 2 + TE2);

%% Set "slices_per_slab"
slices_per_slab = base_resolution; % Slices per slab

%% Calculate the number of readouts per interleaf
nr_readouts = floor(nr_radial_views / nr_interleaves); % number of readouts per interleaf
nr_radial_views = nr_readouts * nr_interleaves;

%% Calculate the bandwidth [Hz/Px]
%--------------------------------------------------------------------------
% real dwell time [sec]
% IDEA p219: dRealDwellTime denotes the dwell time with oversampling
%--------------------------------------------------------------------------
bandwidth1 = 1 / (real_dwell_time1 * (readout_os_factor * base_resolution));
bandwidth2 = 1 / (real_dwell_time2 * (readout_os_factor * base_resolution));

%% Define the filename of a .txt file
filename_base = sprintf('bstar_%3.0fmm_%4.2fmm_true%d_TR%3.2fms_rf%3.0f_i%d_%1.0fk_FA%d_%s_self%d_fid%d_noise%d', ...
    fov_read * 1e3, resolution * 1e3, true_resolution_flag, TR * 1e3, rf_length * 1e6, nr_interleaves, nr_radial_views * 1e-3, flip_angle, trajectory, self_navigation_flag, fid_navigator_flag, noise_scan_flag);

%% Save input parameters as a .txt file
txt_file = fullfile(output_path, sprintf('%s.txt', filename_base));
[fid,message] = fopen(txt_file, 'w');

fprintf(fid, '----------------------- Contrast (Common) ------------------------\n');
fprintf(fid, 'TR         = %4.2f [ms]\n', TR * 1e3);
fprintf(fid, 'TE1        = %4.2f [ms]\n', TE1 * 1e3);
fprintf(fid, 'TE2        = %4.2f [ms]\n', TE2 * 1e3);
fprintf(fid, 'Flip angle = %4d [deg]\n', flip_angle);
fprintf(fid, '------------------------------------------------------------------\n');

fprintf(fid, '----------------------- Resolution (Common) ----------------------\n');
fprintf(fid, 'FoV read        = %5.0f [mm]\n', fov_read * 1e3);
fprintf(fid, 'FoV phase       = %5.1f [%%]\n', 100.0);
fprintf(fid, 'Base resolution = %5d\n', base_resolution);
fprintf(fid, 'Radial views    = %5d\n', nr_radial_views);
fprintf(fid, '------------------------------------------------------------------\n');

fprintf(fid, '----------------------- Sequence (Part 1) ------------------------\n');
fprintf(fid, 'Contrasts       = %7d\n', nr_echoes);
fprintf(fid, 'Bandwidth 1     = %7.0f [Hz/Px]\n', bandwidth1);
fprintf(fid, 'Bandwidth 2     = %7.0f [Hz/Px]\n', bandwidth2);
fprintf(fid, 'Readout mode    = %s\n', readout_mode);
fprintf(fid, '------------------------------------------------------------------\n');

fprintf(fid, '----------------------- Sequence (Special) -----------------------\n');
fprintf(fid, 'Trajectory       = %11s      RF length  = %4.0f [us]\n', trajectory, rf_length * 1e6);
fprintf(fid, 'Readout mode     = %11s      RF dummies = %4.0f\n', 'Combined', rf_dummies);
fprintf(fid, 'F-mod            = %11.1f      RF phase   = %4d [deg]\n', [], rf_phase);
fprintf(fid, 'Tilt             = %11.1f\n', []);
fprintf(fid, 'Interleaves      = %11d\n', nr_interleaves);
fprintf(fid, 'Gradient factor  = %11.2f      TR delay   = %4.0f [us]\n', [], (TR - (rf_length / 2 + TE2)) * 1e6);
fprintf(fid, 'Amplitude factor = %11.2f\n', []);
fprintf(fid, 'Spoiler factor   = %11.2f\n', []);
fprintf(fid, 'Slew rate        = %11.2f [mT/m/ms]\n', slew_rate);
fprintf(fid, 'Amplitude        = %11.2f [mT/m]\n', amplitude);
fprintf(fid, 'True resolution  = %11d \n', true_resolution_flag);
fprintf(fid, '------------------------------------------------------------------\n');
fprintf(fid, 'RUTime/TotalTime/ADC[0]/Resolution: %3.0f/%3.0f/%6.1f/%4.2f\n', ramp_time * 1e6, total_time * 1e6, adc_duration * 1e6, resolution * 1e3);
fprintf(fid, 'Real dwell time (echo 1) = %4.0f [nsec]\n', real_dwell_time1 * 1e9);
fprintf(fid, 'Real dwell time (echo 2) = %4.0f [nsec]\n', real_dwell_time2 * 1e9);

fclose(fid);

%% Display parameters
[fid,message] = fopen(txt_file, 'r');
txt = fread(fid, inf, '*char').';
fprintf('%s', txt);
fclose(fid);

%% Check the validity of imaging parameters
assert(slew_rate <= sys.maxSlew);
assert(ramp_time * 2 <= total_time);
