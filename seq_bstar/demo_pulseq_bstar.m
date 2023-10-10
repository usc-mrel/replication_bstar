% demo_pulseq_bstar.m
% Written by Nam Gyun Lee
% Email: namgyunl@usc.edu, ggang56@gmail.com (preferred)
% Started: 07/10/2022, Last modified: 10/09/2023

%--------------------------------------------------------------------------
% TWIX MUST BE ACQUIRED WITH
% Sequence => Part 1 => Dimension 2D
%--------------------------------------------------------------------------

%% Clean slate
close all; clear all; clc;

%% Set source directories
package_directory = 'D:\replication_bstar';
pulseq_directory = 'D:\pulseq\pulseq';

%% Add source directories to search path
addpath(genpath(package_directory));
addpath(genpath(pulseq_directory));

%% Load an input file
input_pulseq_bstar_protocol_20220727_1_38ms_1_60mm; % breath-hold

%% Calculate bSTAR imaging parameters
calculate_bstar_imaging_parameters;
pause;

%% Make an output directory
mkdir(output_directory);

%% Create a sequence object
seq = mr.Sequence(sys);
start_time = tic;

%% bSTAR pulse sequence
%--------------------------------------------------------------------------
%
%                      TR = 50 + 80 + 1100 + 110 = 1340 usec
%         |<---------------------------------------------------->|
%         |                                                      |
%         |                 TE2 = 1180 usec                      |
%         |   |<-------------------------------------->|         |
%         |   | TE1 = 80 usec      1100 usec           |110 usec |
%         |   |<------>|<----------------------------->|<------->|
%         |   |    rfRingdownTime = 20 usec             TR delay
%         |   |   |<>| |                               |adcDeadTime = 10 usec
%         |___|___|->| |<- adcDeadTime = 10 usec       | |        _______
%         |   |   |  | |   _________                   | |       |       |
%         |   |   |  | |  /|        \                  | |       |       |
%         |   |   |  | | / |         \                 | |       |       |
%         |   |   |  | |/  |          \                | |       |       |
% o-------+---+---+--+-/---|-----------\---------------/---------+-------+-> t
% |<----->|<----->|  | |<->|           |\             /| |<----->|<----->|
% deadTime| RF length| |RampTime       | \           / | deadTime| RF length
%  = 100  |=100 usec | |               |  \_________/  | | = 100    = 100 usec
% |<---------------->| |<------------->|               | |
%        delayTE     | | | TotalTime   |             | | |
%      = 220 usec    | | |                           | | |
%                    | | |xxxxxxxxxxxxxxxxxxxxxxxxxxx| | |
%                    | | |            ADC              | |
%                    |>| |<- shift_adc                 | |
%                    |<--------------------------------->|
%                    | |            delayTR            | |
% |<---------------->|<--------------------------------->|
%       block 1      |              block 2            | |
% o-------+-------+--+-+-------------------------------+-+-------+---------> t
% 0      100    200 220 230                        1330 1340   1440
%
%%--------------------------------------------------------------------------

%% Create an alpha-degree hard pulse event
rf = mr.makeBlockPulse(flip_angle * pi / 180, sys, 'Duration', rf_length);

%% Create an alpha/2-degree hard pulse event
rf_half = mr.makeBlockPulse(flip_angle / 2 * pi / 180, sys, 'Duration', rf_length);

%% Calculate timing (need to decide on the block structure already)
delayTE = round((rf.deadTime + rf_length / 2 + TE1 - sys.adcDeadTime) / sys.gradRasterTime) * sys.gradRasterTime;
delayTR = round((TR - delayTE) / sys.gradRasterTime) * sys.gradRasterTime;

%% Create a base bipolar gradient waveform in the RO direction
%--------------------------------------------------------------------------
% Create a positive trapezoid gradient event
%--------------------------------------------------------------------------
% gamma * amplitude: [Hz/T] * [mT/m] * [T/1e3 mT] => *1e-3 [Hz/m]
g_positive = mr.makeTrapezoid('x', 'riseTime', ramp_time, 'flatTime', total_time - 2 * ramp_time, 'fallTime', ramp_time, 'amplitude', sys.gamma * amplitude * 1e-3);

%--------------------------------------------------------------------------
% Create a negative trapezoid gradient event
%--------------------------------------------------------------------------
g_negative = mr.scaleGrad(g_positive,-1);
g_negative.delay = mr.calcDuration(g_positive);

%--------------------------------------------------------------------------
% Combine a positive gradient event and a negative gradient event
%--------------------------------------------------------------------------
g_bipolar = mr.addGradients({g_positive, g_negative}, sys);
g_bipolar.delay = sys.adcDeadTime;

%--------------------------------------------------------------------------
% Calculate a bipolar gradient waveform (PE-RO-SL order)
%--------------------------------------------------------------------------
grad_samples = length(g_bipolar.waveform);
bipolar_waveform = zeros(3, grad_samples, 'double');
bipolar_waveform(2,:) = g_bipolar.waveform; % PE-RO-SL order

%% Create a bipolar gradient event in the RO direction ([PE,RO,SL] = [y,x,z] in Pulseq)
gx_bipolar = g_bipolar;
gx_bipolar.channel = 'x';

%% Create a bipolar gradient event in the PE direction ([PE,RO,SL] = [y,x,z] in Pulseq)
gy_bipolar = g_bipolar;
gy_bipolar.channel = 'y';

%% Create a bipolar gradient event in the SL direction ([PE,RO,SL] = [y,x,z] in Pulseq)
gz_bipolar = g_bipolar;
gz_bipolar.channel = 'z';

%% Create an ADC event
% NOT WORKING?? A BUG?? Here, adc_delay must be a multiple of 1 us (=1000 ns) instead of 100 ns. 
%shift_adc = sys.adcDeadTime + round(((TE2 - TE1) - adc_duration) / 2 / sys.adcRasterTime) * sys.adcRasterTime; % [sec] 
shift_adc = round((2 * total_time - adc_duration) / 2 / (sys.adcRasterTime * 10)) * (sys.adcRasterTime * 10); % [sec]
adc_delay = sys.adcDeadTime + shift_adc;
adc = mr.makeAdc(adc_samples, 'Dwell', real_dwell_time, 'delay', adc_delay, 'system', sys);

%% Calculate 3D trajectories
%--------------------------------------------------------------------------
% Calculate a spiral phyllotaxis pattern (half-spokes, Delacoste 2018 MRM)
%--------------------------------------------------------------------------
if strcmp(trajectory, 'Phyllotaxis')
    [phi, theta] = calculate_spiral_phyllotaxis_pattern(nr_readouts, nr_interleaves, flag_self);
    if flag_self
        nr_projections_per_interleaf = nr_readouts + 1;
    else
        nr_projections_per_interleaf = nr_readouts;
    end

%--------------------------------------------------------------------------
% Calculate a wobbling Archimedean spiral pole (WASP) pattern
%--------------------------------------------------------------------------
elseif strcmp(trajectory, 'WASP')
    rpa_range = [0 10]; % [degree]
    [phi, theta] = calculate_wasp_pattern(nr_readouts, nr_interleaves, rpa_range);
    nr_projections_per_interleaf = nr_readouts;
end

%% Calculate rotated bipolar gradient waveforms in the GCS [Hz/m] [PE,RO,SL]
% Note that [y,x,z] in Pulseq corresponds to [PE,RO,SL]!
% [PE,RO,SL] = [y,x,z] in Pulseq
%--------------------------------------------------------------------------
%         ^ RO(2)                   ^ RO(2)                  ^ RO(2)
%         |                  theta  |                        |
%         |                       \ |                        | /|
%         |                        \|                        |/ |
%         +---------> PE(1) =>      +---------> PE(1)        +--|------> PE(1)
%        /                         /                        / \ |
%       /                         /                        /phi\|
%      v SL(3)                   v SL(3)                  v SL(3)
%
% All rotation matrices are right-handed.
% 1. Apply rotation about the PE direction by theta (polar)
% 2. Apply rotation about the RO direction by phi (azimuthal)
%
% R_PE(theta) = [1     0           0     ], R_RO(phi) = [ cos(phi) 0 sin(phi)]  
%               [0 cos(theta) -sin(theta)]              [    0     1     0   ]
%               [0 sin(theta)  cos(theta)]              [-sin(phi) 0 cos(phi)]
%--------------------------------------------------------------------------
g_gcs = zeros(3, grad_samples, nr_projections_per_interleaf * nr_interleaves, 'double');

for i = 1:(nr_projections_per_interleaf * nr_interleaves)
    cos_theta = cos(theta(i));
    sin_theta = sin(theta(i));
    cos_phi   = cos(phi(i));
    sin_phi   = sin(phi(i));

    %----------------------------------------------------------------------
    % Calculate the rotation matrix about the PE direction by theta
    %----------------------------------------------------------------------
    R_PE = [1        0             0     ;
            0    cos_theta    -sin_theta ;
            0    sin_theta     cos_theta];

    %----------------------------------------------------------------------
    % Calculate the rotation matrix about the RO direction by phi
    %----------------------------------------------------------------------
    R_RO = [ cos_phi    0    sin_phi ;
                0       1       0    ;
            -sin_phi    0    cos_phi];

    %----------------------------------------------------------------------
    % Calculate rotated bipolar waveforms
    %----------------------------------------------------------------------
    g_gcs(:,:,i) = R_RO * R_PE * bipolar_waveform;
end

%% Flip the sign of gradients in the PE and SL directions [PE,RO,SL]
%--------------------------------------------------------------------------
% The "old/compat” option maps Pulseq logical X, Y, and Z axes to the three 
% axes of Siemens% logical coordinate system (PE-RO-SL) and uses Siemens native 
% transformation from the logical gradient waveforms to the physical gradient
% waveforms. This also involves scaling all gradient axes by -1 followed by 
% additional scaling on the readout direction by -1. To counter these scaling 
% operations performed in the Pulseq interpreter, our code prepares a Pulseq 
% file by intentionally scaling gradient waveforms along the PE and SL 
% directions by -1
%--------------------------------------------------------------------------
g_gcs(1,:,:) = -g_gcs(1,:,:); % PE
g_gcs(3,:,:) = -g_gcs(3,:,:); % SL

%% Pre-register objects that do not change while looping
[~, rf.shapeIDs] = seq.registerRfEvent(rf); % the phase of the RF object will change, therefore we only per-register the shapes

%% Define sequence blocks for an alpha/2-TR/2 sequence
tstart = tic; fprintf('Defining blocks for an alpha/2-TR/2 sequence... ');
rf_phase_rad = rf_phase * pi / 180;
count = 1;

%--------------------------------------------------------------------------
% Set the phase of an RF event and an ADC event
%--------------------------------------------------------------------------
rf_half.phaseOffset = rf_phase_rad * mod(count,2);
adc.phaseOffset = rf_phase_rad * mod(count,2);

%--------------------------------------------------------------------------
% Add a new block to the sequence
%--------------------------------------------------------------------------
seq.addBlock(rf_half, mr.makeDelay(TR/2));
count = count + 1;
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%% Define sequence blocks for dummy pulses
%--------------------------------------------------------------------------
% Create a dummy gradient event in the PE direction (PRS)
% [PE,RO,SL] = [y,x,z] in Pulseq
%--------------------------------------------------------------------------
gy_dummy = gy_bipolar;
gy_dummy.waveform = g_gcs(1,:,1); % PE => 'y'

%--------------------------------------------------------------------------
% Create a dummy gradient event in the RO direction (PRS)
% [PE,RO,SL] = [y,x,z] in Pulseq
%--------------------------------------------------------------------------
gx_dummy = g_bipolar;
gx_dummy.waveform = g_gcs(2,:,1); % RO => 'x'

%--------------------------------------------------------------------------
% Create a dummy gradient event in the SL direction (PRS)
% [PE,RO,SL] = [y,x,z] in Pulseq
%--------------------------------------------------------------------------
gz_dummy = gz_bipolar;
gz_dummy.waveform = g_gcs(3,:,1); % SL => 'z'

for i = 1:rf_dummies
    tstart = tic; fprintf('Defining blocks for dummy pulses (%3d/%3d)... ', i, rf_dummies);

    %----------------------------------------------------------------------
    % Set the phase of an RF event and an ADC event
    %----------------------------------------------------------------------
    rf.phaseOffset = rf_phase_rad * mod(count,2);
    adc.phaseOffset = rf_phase_rad * mod(count,2);

    %----------------------------------------------------------------------
    % Add a new block to the sequence (block 1)
    %----------------------------------------------------------------------
    seq.addBlock(rf, mr.makeDelay(delayTE));

    %----------------------------------------------------------------------
    % Add a new block to the sequence (block 2)
    %----------------------------------------------------------------------
    seq.addBlock(gx_dummy, gy_dummy, gz_dummy, mr.makeDelay(delayTR));
    count = count + 1;
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
end

%% Define sequence blocks for MR data acquisition
for i = 1:(nr_projections_per_interleaf * nr_interleaves)
    tstart = tic; fprintf('Defining blocks for MR data acquisition (%5d/%5d)... ', i, (nr_projections_per_interleaf * nr_interleaves));

    %----------------------------------------------------------------------
    % Set the phase of an RF event and an ADC event
    %----------------------------------------------------------------------
    rf.phaseOffset = rf_phase_rad * mod(count,2);
    adc.phaseOffset = rf_phase_rad * mod(count,2);

    %----------------------------------------------------------------------
    % Create a bipolar gradient event in the PE direction (PRS)
    % [PE,RO,SL] = [y,x,z] in Pulseq
    %----------------------------------------------------------------------
    gy_bipolar.waveform = g_gcs(1,:,i); % PE => 'y'

    %----------------------------------------------------------------------
    % Create a bipolar gradient event in the RO direction (PRS)
    % [PE,RO,SL] = [y,x,z] in Pulseq
    %----------------------------------------------------------------------
    g_bipolar.waveform = g_gcs(2,:,i); % RO => 'x'

    %----------------------------------------------------------------------
    % Create a bipolar gradient event in the SL direction (PRS)
    % [PE,RO,SL] = [y,x,z] in Pulseq
    %----------------------------------------------------------------------
    gz_bipolar.waveform = g_gcs(3,:,i); % SL => 'z'

    %----------------------------------------------------------------------
    % Add a new block to the sequence (block 1)
    %----------------------------------------------------------------------
    seq.addBlock(rf, mr.makeDelay(delayTE));

    %----------------------------------------------------------------------
    % Add a new block to the sequence (block 2)
    %----------------------------------------------------------------------
    seq.addBlock(g_bipolar, gy_bipolar, gz_bipolar, adc, mr.makeDelay(delayTR));
    count = count + 1;
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
end

%% check whether the timing of the sequence is correct
[ok, error_report] = seq.checkTiming;

if (ok)
    fprintf('Timing check passed successfully\n');
else
    fprintf('Timing check failed! Error listing follows:\n');
    fprintf([error_report{:}]);
    fprintf('\n');
end

%% Prepare sequence export
seq.setDefinition('Amplitude', amplitude);
seq.setDefinition('AmplitudeFactor', amplitude_factor);
seq.setDefinition('Bandwidth1', bandwidth1);
seq.setDefinition('Bandwidth2', bandwidth2);
seq.setDefinition('BaseResolution', base_resolution);
seq.setDefinition('Contrasts', nr_echoes);
seq.setDefinition('F-mod', f_mod);
seq.setDefinition('FOV', [fov_read fov_read fov_read]);
seq.setDefinition('FlipAngle', flip_angle);
seq.setDefinition('GradientFactor', gradient_factor);
seq.setDefinition('Interleaves', nr_interleaves);
seq.setDefinition('Name', 'bSTAR');
seq.setDefinition('RFDummies', rf_dummies);
seq.setDefinition('RFLength', rf_length);
seq.setDefinition('RFPhase', rf_phase);
seq.setDefinition('RadialViews', nr_radial_views);
seq.setDefinition('RampTime', ramp_time);
seq.setDefinition('ReadoutMode', readout_mode);
seq.setDefinition('ReadoutOSFactor', readout_os_factor);
seq.setDefinition('Readouts', nr_readouts);
seq.setDefinition('RealDwellTime1', real_dwell_time1);
seq.setDefinition('RealDwellTime2', real_dwell_time2);
seq.setDefinition('Resolution', resolution);
seq.setDefinition('SelfNavigation', flag_self);
seq.setDefinition('SlewRate', slew_rate);
seq.setDefinition('SpoilerFactor', spoiler_factor);
seq.setDefinition('TE1', TE1);
seq.setDefinition('TE2', TE2);
seq.setDefinition('TR', TR);
seq.setDefinition('Tilt', tilt);
seq.setDefinition('TotalTime', total_time);
seq.setDefinition('Trajectory', trajectory);

seq_filename = sprintf('bstar_ecg%d_TR%3.2fms_%3.2fmm_b%3.0f_rf%3.0f_i%d_%1.0fk_FA%d_self%d_%s.seq', 0, TR * 1e3, resolution * 1e3, bandwidth1, rf_length * 1e6, nr_interleaves, nr_radial_views * 1e-3, flip_angle, flag_self, trajectory);
seq_file = fullfile(output_directory, seq_filename);
seq.write(seq_file);   % Output sequence for scanner

%% Plot sequence and k-space diagrams
seq.plot('timeRange', [0 1200] * TR);

if 0
% k-space trajectory calculation
[ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing] = seq.calculateKspacePP();

% plot 3D k-space
figure; plot3(ktraj(1,:), ktraj(2,:), ktraj(3,:), 'b'); % a 3D k-space plot
axis('equal'); % enforce aspect ratio for the correct trajectory display
hold; plot3(ktraj_adc(1,:), ktraj_adc(2,:), ktraj_adc(3,:), 'r.'); % plot the sampling points
title('full k-space trajectory (k_x x k_y x k_z)');

%% very optional slow step, but useful for testing during development e.g. for the real TE, TR or for staying within slewrate limits  
rep = seq.testReport;
fprintf([rep{:}]);
end
