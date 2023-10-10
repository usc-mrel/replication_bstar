% demo_pulseq_Zhao_2020_MRM.m
% Written by Nam Gyun Lee
% Email: namgyunl@usc.edu, ggang56@gmail.com (preferred)
% Started: 09/12/2022, Last modified: 09/26/2022

%% Clean slate
close all; clear all; clc;

%% Set source directories
pulseq_directory = 'D:\pulseq\pulseq';

%% Add source directories to search path
addpath(genpath(pulseq_directory));

%% Define imaging parameters
%--------------------------------------------------------------------------
% Sequence (Special)
%--------------------------------------------------------------------------
TR              = 500;     % TR [msec]
TE              = 10;      % TE [msec]
nr_averages     = 1;       % Averages
flip_angle      = 90;      % Flip angle [deg]
nr_repetitions  = 2;       % gradient polarity: 1=positive only, 2=positive and negative

%--------------------------------------------------------------------------
% Resolution (Common)
%--------------------------------------------------------------------------
vector_size = 512; % Vector size (number of samples without oversampling)

%--------------------------------------------------------------------------
% Sequence (Common)
%--------------------------------------------------------------------------
nr_echoes           = 2;      % Contrasts
nr_dummies          = 0;      % Preparation scans
phase_cycling       = 'None'; % Phase cycling
bandwidth           = 100000; % Bandwidth [Hz]
remove_oversampling = 0;      % Remove oversampling: 1=yes, 0=no

%--------------------------------------------------------------------------
% Sequence (Special)
%--------------------------------------------------------------------------
ramp_time                 = 140;     % RampTime        [usec]
total_time                = 520;     % TotalTime       [usec]
adc_grad_delay            = 500;     % ADC-Grad Delay  [usec]
amplitude                 = 6;       % Amplitude       [mT/m]
max_slew                  = 137.06;  % Max SlewRate    [mT/m/ms]
rf_duration               = 8200;    % RF duration     [usec]
slice_thickness           = 5;       % Slice Thickness [mm]
fov                       = 220;     % Field of view   [mm]
nr_phase_encoding_steps_1 = 64;      % number of phase-encoding steps in the readout direction        (ISMRMRD: phase encoding line number)
nr_phase_encoding_steps_2 = 64;      % number of phase-encoding steps in the phase-encoding direction (ISMRMRD: partition encoding number)

%% Convert the units of input parameters
TR = TR * 1e-3; % [msec] * [sec/1e3 msec] => *1e-3 [sec]
TE = TE * 1e-3; % [msec] * [sec/1e3 msec] => *1e-3 [sec]

ramp_time       = ramp_time * 1e-6;       % [usec] * [sec/1e6 usec] => *1e-6 [sec]
total_time      = total_time * 1e-6;      % [usec] * [sec/1e6 usec] => *1e-6 [sec]
adc_grad_delay  = adc_grad_delay * 1e-6;  % [usec] * [sec/1e6 usec] => *1e-6 [sec]
rf_duration     = rf_duration * 1e-6;     % [usec] * [sec/1e6 usec] => *1e-6 [sec]
slice_thickness = slice_thickness * 1e-3; % [mm] * [m/1e3 mm] => *1e-3 [m]
fov             = fov * 1e-3;             % [mm] * [m/1e3 mm] => *1e-3 [m]

%% Calculate dependent parameters
%--------------------------------------------------------------------------
% readout oversampling factor
%--------------------------------------------------------------------------
if remove_oversampling
    readout_os_factor = 1;
else
    readout_os_factor = 2;
end

%--------------------------------------------------------------------------
% real dwell time [sec]
% IDEA p219: dRealDwellTime denotes the dwell time with oversampling
%--------------------------------------------------------------------------
% round-up dwell time to 100 ns (sys.adcRasterTime  = 100 ns)
real_dwell_time = round((1 / bandwidth) / readout_os_factor * 1e7) * 1e-7;

%--------------------------------------------------------------------------
% Acquisition duration [sec]
%--------------------------------------------------------------------------
adc_duration = real_dwell_time * vector_size * readout_os_factor;

%% Define an output directory
output_directory = pwd;

%% Display input parameters
fprintf('----------------------- Sequence (Common) -----------------------\n');
fprintf('Preparation scans = %3d      Phase cycling        = %s\n', nr_dummies, phase_cycling);
fprintf('                             [\bAcquisition delay    = %6.1f [ms]]\b\n', 0);
fprintf('                             Bandwidth            = %6.0f [Hz]\n', bandwidth);
fprintf('                             [\bAcquisition duration = %6.0f [ms]]\b\n', adc_duration * 1e3);
fprintf('------------------------------------------------------------------\n');

fprintf('----------------------- Sequence (Special) -----------------------\n');
fprintf('RampTime       = %4.0f [us]      RF Duration     = %4.0f [us]\n', ramp_time * 1e6, rf_duration * 1e6);
fprintf('TotalTime      = %4.0f [us]      Slice Thickness = %4.2f [mm]\n', total_time * 1e6, slice_thickness * 1e3);
fprintf('ADC-Grad Delay = %4.0f [us]      Slice Shift     = %4.0f [mm]\n', adc_grad_delay * 1e6, 0 * 1e3);
fprintf('Amplitude      = %4.1f [mT/m]\n', amplitude);
fprintf('[\bDwell Time     = %4.0f [ns]]\b\n', real_dwell_time * 1e9);
fprintf('Max SlewRate   = %4.0f [mT/m/ms]\n', max_slew);
fprintf('------------------------------------------------------------------\n');
pause;

%% Define MRI scanner specifications
B0 = 0.55; % main field strength [T]
grad_mode = 'whisper';

switch grad_mode
    case 'fast'
        max_grad = 24;      % Max gradient strength [mT/m]
        max_slew = 180.18;  % Maximum slew rate [mT/m/ms]
    case 'normal'
        max_grad = 22;      % Max gradient strength [mT/m]
        max_slew = 100;     % Maximum slew rate [mT/m/ms]
    case 'whisper'
        max_grad = 22;      % Max gradient strength [mT/m]
        max_slew = 50;      % Maximum slew rate [mT/m/ms]
end

%% Set system limits
sys = mr.opts('MaxGrad'       , max_grad, 'GradUnit', 'mT/m' , ...
              'MaxSlew'       , max_slew, 'SlewUnit', 'T/m/s', ...
              'rfRingdownTime', 20e-6 , ...
              'rfDeadtime'    , 100e-6, ...
              'adcDeadTime'   , 10e-6 , ...
              'B0'            , B0);

%% Pulse sequence
%-------------------------------------------------------------------------------
%                                     TR = 500,000 usec
%    |<----------------------------------------------------------->|
%    |                  TE = 6000 usec                             |
%    |         |<------------------------>|                        |
%    | |       |       |        |      |  |                        |
%    | |      _|_ rf   |        |      |  |                        |
%    | |     / | \     |        |      |  |                        |
% ___|_|    /  |  \    |________|______|__|________________________|
%    | |\__/   |   \__/|        |      |  |                        |
%    | |       |       |        |      |  |                        |
%    | |       |       |        |      |  |                        |
% Gz | |_______|_______| |      |      |  |                        |
%    | /       | gz_ss \ |      |      |  |                        |
% ___|/        |        \|      |______|__|________________________|
%    |         |         \     /|      |  |                        |
%    |         |         |\___/ |      |  |                        |
%    |         |         |gz_reph|     |  |ADC-Grad delay          |
%    |<----------------->|<---->|      |  |= 500 usec              |
%    |     4020 usec     |      |      |  |    |TotalTime          |
%    |                   |      |      |  |    |= 1100usec         |
%    |                   |  __  |      |  |<-->|<-------->|        |
% Gx |                   | /__\ |      |  |    | ___      |        |
% ___|___________________|/____\|______|__|____|/   \     |________|
%    |                   |\ __ /|      |  |    |     \___/|        |
%    |                   | \__/ |      |  |    |gx_bipolar|        |
%    |                   |      |      |  |    |          |        |
%    |                   |  __  |      |  |    |                   |
% Gy |                   | /__\ |      |  |    |                   |
% ___|___________________|/____\|______|__|____|___________________|
%    |                   |\ __ /|      |  |    |                   |
%    |                   | \__/ |      |  |    |                   |
%    |                   |      |      |  |    |                   |
%    |<----------------->|<----------->|<------------------------->|
%    |  mr.calcDuration  |  delayTE    |  |    |  delayTR          |
%    |      (gz_ss)      |             |   ______________________  |
%    |         |         |             |  |    |    ADC          | |
% ___|_________|_________|_____________|__|    |                 |_|
%    |         |         | adcDeadTime |<>|    |                   |
%    |         |         | = 10 usec   |  |    |                   |
%    |<----------------->|<----------->|<------------------------->|
%           block 1          block 2   |  |    |   block 3         |
%----o---------+---------+-------------+--+----+--------------------> t
%            2010      4020         8000 8010 8510
%-------------------------------------------------------------------------------

%% Create an alpha-degree slice selection pulse and gradient
[rf, gz_ss, gz_reph] = mr.makeSincPulse(flip_angle * pi / 180, 'Duration', rf_duration, ...
    'sliceThickness', slice_thickness, 'apodization', 0.5, 'timeBwProduct', 4, 'system', sys);

%% Define phase encoding gradients ([PE,RO,SL] = [y,x,z] in Pulseq)
deltak = 1 / fov; % [cycle/cm]
pe1_steps = (-floor(nr_phase_encoding_steps_1/2):ceil(nr_phase_encoding_steps_1/2)-1).' * deltak; % RO (= Gx)
pe2_steps = (-floor(nr_phase_encoding_steps_2/2):ceil(nr_phase_encoding_steps_2/2)-1).' * deltak; % PE (= Gy)

%% Flip the sign of the PE direction [PE,RO,SL]
% This step is necessary to make use of coordinate transformations in Siemens
% for datasets acquired with Pulseq.
pe2_steps = -pe2_steps; % PE

%% Create bipolar readout gradients ([PE,RO,SL] = [y,x,z] in Pulseq)
%--------------------------------------------------------------------------
% Create a trapezoid event
%--------------------------------------------------------------------------
% gamma * amplitude: [Hz/T] * [mT/m] * [T/1e3 mT] => *1e-3 [Hz/m]
gx_positive = mr.makeTrapezoid('x', 'riseTime', ramp_time, 'flatTime', total_time - 2 * ramp_time, ...
    'fallTime', ramp_time, 'amplitude', sys.gamma * amplitude * 1e-3);

%--------------------------------------------------------------------------
% Calculate a bipolar event
%--------------------------------------------------------------------------
gx_negative = mr.scaleGrad(gx_positive,-1);
gx_negative.delay = mr.calcDuration(gx_positive);
gx_bipolar = mr.addGradients({gx_positive, gx_negative}, sys);
gx_bipolar.delay = sys.adcDeadTime + adc_grad_delay;

%% Create an ADC readout event
adc_samples = vector_size * readout_os_factor;
adc_delay = sys.adcDeadTime;
adc = mr.makeAdc(adc_samples, 'Dwell', real_dwell_time, 'delay', adc_delay, 'system', sys);

%% Calculate timing (need to decide on the block structure already)
delayTE = round((TE - (gz_ss.flatTime / 2 + gz_ss.fallTime + sys.adcDeadTime)) / sys.gradRasterTime) * sys.gradRasterTime;
delayTR = round((TR - (mr.calcDuration(gz_ss) + delayTE)) / sys.gradRasterTime) * sys.gradRasterTime;

%% Create a sequence object
seq = mr.Sequence(sys);
start_time = tic;

%% Define sequence blocks
% all LABELS / counters and flags are automatically initialized to 0 in the beginning, no need to define initial 0's  
% so we will just increment LIN after the ADC event (e.g. during the spoiler)
%--------------------------------------------------------------------------
% ISMRMRD header
% uint16_t kspace_encode_step_1;    /**< e.g. phase encoding line number */
% uint16_t kspace_encode_step_2;    /**< e.g. partition encoding number */
% uint16_t average;                 /**< e.g. signal average number */
% uint16_t slice;                   /**< e.g. imaging slice number */
% uint16_t contrast;                /**< e.g. echo number in multi-echo */
% uint16_t phase;                   /**< e.g. cardiac phase number */
% uint16_t repetition;              /**< e.g. dynamic number for dynamic scanning */
% uint16_t set;                     /**< e.g. flow encoding set */
% uint16_t segment;                 /**< e.g. segment number for segmented acquisition */
% uint16_t user[ISMRMRD_USER_INTS]; /**< Free user parameters */
%--------------------------------------------------------------------------
lbl_inc_lin   = mr.makeLabel('INC', 'LIN', 1); % lin == line
lbl_inc_par   = mr.makeLabel('INC', 'PAR', 1); % par == partition
lbl_inc_avg   = mr.makeLabel('INC', 'AVG', 1); % avg == average
lbl_inc_rep   = mr.makeLabel('INC', 'REP', 1); % rep == repetition
lbl_reset_lin = mr.makeLabel('SET', 'LIN', 0);
lbl_reset_par = mr.makeLabel('SET', 'PAR', 0);
lbl_reset_avg = mr.makeLabel('SET', 'AVG', 0);
lbl_reset_rep = mr.makeLabel('SET', 'REP', 0);

%--------------------------------------------------------------------------
% Repetition (REP)
%--------------------------------------------------------------------------
for idx4 = 1:nr_repetitions

    %----------------------------------------------------------------------
    % Set the polarity of a target gradient
    %----------------------------------------------------------------------
    g_test = gx_bipolar;
    if mod(idx4-1,2) == 0     % positive gradient at off-center
        g_test = mr.scaleGrad(g_test,1);
    elseif mod(idx4-1,2) == 1 % negative gradient at off-center
        g_test = mr.scaleGrad(g_test,-1);
    end

    %----------------------------------------------------------------------
    % 2nd phase ecnoding (PAR)
    %----------------------------------------------------------------------
    for idx2 = 1:nr_phase_encoding_steps_2
        gy_pre = mr.makeTrapezoid('y', 'Area', pe2_steps(idx2), 'Duration', mr.calcDuration(gz_reph), 'system', sys);

        %------------------------------------------------------------------
        % 1st phase ecnoding (LIN)
        %------------------------------------------------------------------
        for idx1 = 1:nr_phase_encoding_steps_1
            gx_pre = mr.makeTrapezoid('x', 'Area', pe1_steps(idx1), 'Duration', mr.calcDuration(gz_reph), 'system', sys);

            %--------------------------------------------------------------
            % Average (AVG)
            %--------------------------------------------------------------
            for idx3 = 1:nr_averages
                tstart = tic; fprintf('(REP=%d/%d) Defining sequence blocks (LIN=%2d/%2d, PAR=%2d/%2d, AVG=%2d/%2d)... ', idx4, nr_repetitions, idx1, nr_phase_encoding_steps_1, idx2, nr_phase_encoding_steps_2, idx3, nr_averages);
                seq.addBlock(rf, gz_ss);
                seq.addBlock(gz_reph, gx_pre, gy_pre, mr.makeDelay(delayTE));
                seq.addBlock(g_test, adc, mr.makeDelay(delayTR));

                %----------------------------------------------------------
                % Update AVG counter
                %----------------------------------------------------------
                if idx3 ~= nr_averages
                    seq.addBlock(lbl_inc_avg);
                else
                    seq.addBlock(lbl_reset_avg);
                end
                fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
            end

            %--------------------------------------------------------------
            % Update LIN counter
            %--------------------------------------------------------------
            if idx1 ~= nr_phase_encoding_steps_1
                seq.addBlock(lbl_inc_lin);
            else
                seq.addBlock(lbl_reset_lin);
            end
        end

        %------------------------------------------------------------------
        % Update PAR counter
        %------------------------------------------------------------------
        if idx2 ~= nr_phase_encoding_steps_2
            seq.addBlock(lbl_inc_par);
        else
            seq.addBlock(lbl_reset_par);
        end
    end

    %----------------------------------------------------------------------
    % Update REP counter
    %----------------------------------------------------------------------
    if idx4 ~= nr_repetitions
        seq.addBlock(lbl_inc_rep);
    else
        seq.addBlock(lbl_reset_rep);
    end
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

%% prepare sequence export
seq.setDefinition('ADCGradDelay', adc_grad_delay);
seq.setDefinition('Amplitude', amplitude);
seq.setDefinition('Averages', nr_averages);
seq.setDefinition('Bandwidth', bandwidth);
seq.setDefinition('Contrasts', nr_echoes);
seq.setDefinition('FlipAngle', flip_angle);
seq.setDefinition('FOV', [fov fov slice_thickness]);
seq.setDefinition('MaxSlewRate', max_slew);
seq.setDefinition('Name', 'Zhao 2020 MRM v2');
seq.setDefinition('PhaseEncodingSteps1', nr_phase_encoding_steps_1);
seq.setDefinition('PhaseEncodingSteps2', nr_phase_encoding_steps_2);
seq.setDefinition('RFDuration', rf_duration);
seq.setDefinition('RampTime', ramp_time);
seq.setDefinition('ReadoutOSFactor', readout_os_factor);
seq.setDefinition('RealDwellTime', real_dwell_time);
seq.setDefinition('Repetitions', nr_repetitions);
seq.setDefinition('TE', TE);
seq.setDefinition('TR', TR);
seq.setDefinition('TotalTime', total_time);
seq.setDefinition('VectorSize', vector_size);

seq_filename = sprintf('Zhao_v2_%3.0fkHz_RT%3.0fus_TT%3.0fus_delay%3.0fus_amp%1.0f_Smax%2.0f_slice%4.2fmm_%dx%d.seq', ...
    bandwidth * 1e-3, ramp_time * 1e6, total_time * 1e6, adc_grad_delay * 1e6, amplitude, max_slew, slice_thickness * 1e3, nr_phase_encoding_steps_1, nr_phase_encoding_steps_2);
seq_path = fullfile(output_directory, seq_filename);
seq.write(seq_path);   % Output sequence for scanner

%% plot sequence and k-space diagrams
seq.plot('timeRange', [0 1] * TR, 'label', 'AVG');

if 0
% k-space trajectory calculation
[ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing] = seq.calculateKspacePP();

% plot k-spaces
figure; plot3(ktraj(1,:), ktraj(2,:), ktraj(3,:), 'b'); % a 3D k-space plot
axis('equal'); % enforce aspect ratio for the correct trajectory display
hold; plot3(ktraj_adc(1,:), ktraj_adc(2,:), ktraj_adc(3,:), 'r.'); % plot the sampling points
title('full k-space trajectory (k_x x k_y x k_z)');

%% very optional slow step, but useful for testing during development e.g. for the real TE, TR or for staying within slewrate limits  
rep = seq.testReport;
fprintf([rep{:}]);
end
