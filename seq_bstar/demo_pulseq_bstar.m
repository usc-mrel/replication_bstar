% demo_pulseq_bstar.m
% Written by Nam Gyun Lee
% Email: namgyunl@usc.edu, ggang56@gmail.com (preferred)
% Started: 04/18/2024, Last modified: 04/25/2024

%--------------------------------------------------------------------------
% TWIX MUST BE ACQUIRED WITH
% Sequence => Part 1 => Dimension 2D
%--------------------------------------------------------------------------

%% Clean slate
close all; clear all; clc;

%% Set source paths
package_path = 'D:\replication_bstar_v2';
pulseq_path = 'D:\pulseq';
ismrmrd_path = 'D:\ismrmrd';

%% Add source paths to search path
addpath(genpath(package_path));
addpath(genpath(pulseq_path))
addpath(genpath(ismrmrd_path));

%% Start a stopwatch timer
start_time = tic;

%% Define the full path of a .json file
%json_file = 'D:\replication_bstar_v2\seq_bstar\pulseq_bstar_sequence_parameters_0_90mm.json';
%json_files{1} = 'D:\replication_bstar_v2\seq_bstar\pulseq_bstar_sequence_parameters_1_60mm_17sec.json';
%json_files{2} = 'D:\replication_bstar_v2\seq_bstar\pulseq_bstar_sequence_parameters_1_60mm_23sec.json';
json_files{1} = 'D:\replication_bstar_v2\seq_bstar\pulseq_bstar_sequence_parameters_0_90mm_free_breathing.json';
nr_json_files = length(json_files);

for json_number = 1:nr_json_files
json_file = json_files{json_number};

%% Read a .json file
tstart = tic; fprintf('%s: Reading a .json file: %s... ', datetime, json_file);
fid = fopen(json_file);
json_txt = fread(fid, [1 inf], 'char=>char'); 
fclose(fid);
json = jsondecode(json_txt);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%% Define imaging parameters
%--------------------------------------------------------------------------
% Contrast (Common)
%--------------------------------------------------------------------------
flip_angle = json.flip_angle; % Flip angle [deg]

%--------------------------------------------------------------------------
% Resolution (Common)
%--------------------------------------------------------------------------
fov_read        = json.fov_read;        % FoV read [m]
resolution      = json.resolution;      % Spatial resolution [m]
nr_readouts     = json.nr_readouts;     % Number of half-spokes per interleaf
nr_interleaves  = json.nr_interleaves;  % Number of interleaves
nr_radial_views = json.nr_radial_views; % Radial views (i.e., total number of readouts)

%--------------------------------------------------------------------------
% Sequence (Part 1)
%--------------------------------------------------------------------------
nr_echoes = json.nr_echoes; % Contrasts

%--------------------------------------------------------------------------
% Sequence (Part 2)
%--------------------------------------------------------------------------
grad_mode    = json.grad_mode;    % Gradient mode: 'Fast', 'Normal' or 'Whisper'
readout_mode = json.readout_mode; % Readout mode

%--------------------------------------------------------------------------
% Sequence (Special)
%--------------------------------------------------------------------------
trajectory     = json.trajectory;     % Trajectory: Archimedean, WASP, CardioWASP, Phyllotaxis
rf_length      = json.rf_length;      % RF length [sec]
rf_dummies     = json.rf_dummies;     % RF dummies
rf_phase       = json.rf_phase;       % RF phase [deg]

%--------------------------------------------------------------------------
% Define a readout oversampling factor
%--------------------------------------------------------------------------
readout_os_factor = json.readout_os_factor;

%--------------------------------------------------------------------------
% Define parameters for an ADC window
%--------------------------------------------------------------------------
discard_pre     = json.discard_pre;     % number of samples to be discarded at the beginning of acquisition
discard_post    = json.discard_post;    % number of samples to be discarded at the end of acquisition
real_dwell_time = json.real_dwell_time; % Real dwell time [sec]

%--------------------------------------------------------------------------
% Define parameters for an FID navigator
%--------------------------------------------------------------------------
fid_navigator_flag = json.fid_navigator_flag; % FID navigator: 0=no, 1=yes
rf_length_nav      = json.rf_length_nav;      % RF duration [sec]
adc_samples_nav    = json.adc_samples_nav;    % number of ADC samples
T_delay            = json.T_delay;            % Delay between the end of an FID navigator and the start of a scan [sec]
TE_nav             = json.TE_nav;             % Navigator TE [sec]
flip_angle_nav     = json.flip_angle_nav;     % Flip angle [deg]
dwell_time_nav     = json.dwell_time_nav;     % Dwell time without readout oversampling [sec]

%--------------------------------------------------------------------------
% Self-navigation
%--------------------------------------------------------------------------
self_navigation_flag = json.self_navigation_flag; % self-navigation (SI projections): 0=no, 1=yes

%--------------------------------------------------------------------------
% True resolution
%--------------------------------------------------------------------------
true_resolution_flag = json.true_resolution_flag;

%--------------------------------------------------------------------------
% Noise scan
%--------------------------------------------------------------------------
noise_scan_flag = json.noise_scan_flag; % noise only scan: 0=no, 1=yes

%--------------------------------------------------------------------------
% Orientation mapping for a Siemens Pulseq interpreter
%--------------------------------------------------------------------------
orientation_mapping = json.orientation_mapping;

%--------------------------------------------------------------------------
% Define main field strength [T]
%--------------------------------------------------------------------------
B0 = json.B0; % [T]

%% Define an output path
[~,output_base,~] = fileparts(json_file);
[y,m,d] = ymd(datetime);
output_path = sprintf('%s_%d%02d%02d', output_base, y, m, d);

%% Make an output path
mkdir(output_path);

%% Calculate bSTAR imaging parameters
calculate_bstar_imaging_parameters;

return

%% Create a sequence object
seq = mr.Sequence(sys);

%% bSTAR pulse sequence
%--------------------------------------------------------------------------------------------------------------|
%                                             TR                                                               |
% |<------------------------------------------------------------------------------------------>|               |
% |           |                                    TE2                                         |               |   
% |           |<------------------------------------------------------------------------->|<---------->|       | 
% |           |    TE1                                                                    |  TR_delay  |       | 
% |           |<--------->|                                                               |    |       |       | 
% |           |    rfRingdownTime = 20 usec                                               |    |       |       |
% |       |   |   |<>|    |                                                               |    |       |   |   |
% |       |___|___|  |    |                                                               |    |       |___|___|
% |       |   |   |  |    |                                                               |    |       |   |   |
% |       |   |   |  |    |                                                               |    |       |   |   |
% | RF    |   |   |  |    |                                                               |    |       |   |   |
% |_______|___|___|__|____|_______________________________________________________________|____|_______|___|___|
% |<----->|<----->|  |<-->| delay_pre                     |                    delay_post |<-->|<----->|<----->|
% deadTime| RF length|    |   _________________________   |                               |    deadTime| RF length
% |                  | |  |  /|                        \  |                               |  | |               |
% |                  | |  | / |                         \ |                               |  | |               |
% |__________________|_|__|/  |                          \|                               |__|_|_______________|
% |                  | |  |<->|                           |\                             /|  | |               |
% |                  | |  |ramp_time                      | \                           / |  | |               |
% |                  | |  |<----------------------------->|  \_________________________/  |  | |               |
% |                  | |  |          total_time           |                               |  | |               |
% |                  | |  |                               |                               |  | |               |
% |                  | |  |                               |                               |  | |               |
% |__________________|_|__|_______________________________|_______________________________|__|_|_______________|
% |                  | |  |                                                               |  | |               |
% |                  | |  |                                                               |  | |               |
% |                  | |  |                                                               |  | |               |
% |    adcDeadTime ->| |<-|<------------------------------------------------------------->|->| |<- adcDeadTime |
% |     = 10 usec    | |  |                         grad_duration                         |  | |               |
% |                  | |xx|xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx|xx| |               |
% |      discard_pre ->|  |<-                            ADC                            ->|  |<- discard_post  |
% |<---------------->|<----------------------------------------------------------------------->|               |
% |     delayTE      |                                 delayTR                                 |               |
% |<---------------->|<----------------------------------------------------------------------->|               |
% |     block 1      |                                 block 2                                 |               |
%-o------------------+-------------------------------------------------------------------------+------------> t|
%--------------------------------------------------------------------------------------------------------------|

%% Create an alpha-degree hard pulse event
rf = mr.makeBlockPulse(flip_angle * pi / 180, sys, 'Duration', rf_length);

%% Create an alpha/2-degree hard pulse event
rf_half = mr.makeBlockPulse(flip_angle / 2 * pi / 180, sys, 'Duration', rf_length);

%% Create a base bipolar gradient waveform in the RO direction [Hz/m]
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
g_bipolar.delay = delay_pre;

%--------------------------------------------------------------------------
% Calculate a bipolar gradient waveform (PE-RO-SL order)
%--------------------------------------------------------------------------
grad_samples = length(g_bipolar.waveform);
g_base_seq = zeros(3, grad_samples, 'double');
g_base_seq(2,:) = g_bipolar.waveform; % PE-RO-SL order

%% Interpolate nominal gradient waveforms in the GCS [Hz/m] [PE,RO,SL]
t_seq = g_bipolar.tt;
t_adc = ((0:adc_samples-1).' + 0.5) * real_dwell_time; % art: ADC Raster Time

g_base_adc = zeros(3, adc_samples, 'double');
g_base_adc(1,:) = interp1(t_seq, g_base_seq(1,:), t_adc, 'linear', 'extrap');
g_base_adc(2,:) = interp1(t_seq, g_base_seq(2,:), t_adc, 'linear', 'extrap');
g_base_adc(3,:) = interp1(t_seq, g_base_seq(3,:), t_adc, 'linear', 'extrap');

%% Display interpolated gradient waveforms
FontSize = 12;

figure('Color', 'w', 'Position', [322 528 918 450]);
ColorOrder = get(gca, 'colororder');
hold on;
plot(t_seq * 1e6, g_base_seq(2,:) / (sys.gamma * 1e-3), '.-', 'MarkerSize', 14, 'Color', ColorOrder(1,:));
plot(t_adc * 1e6, g_base_adc(2,:) / (sys.gamma * 1e-3), '.-', 'MarkerSize', 10, 'Color', ColorOrder(2,:));
grid on; grid minor;
set(gca, 'Box', 'On', 'TickLabelInterpreter', 'latex', 'FontSize', FontSize);
xlabel('Time [usec]', 'Interpreter', 'latex', 'FontSize', FontSize);
ylabel('Amplitude [mT/m]', 'Interpreter', 'latex', 'FontSize', FontSize);
set(gca, 'MinorGridLineStyle', '-');
title({sprintf('RUTime/TotalTime/ADC[0]/Resolution: %3.0f/%3.0f/%6.1f/%4.2f', ramp_time * 1e6, total_time * 1e6, adc_duration * 1e6, resolution * 1e3), ...
       sprintf('amplitude: %3.1f mT/m, slew rate: %6.2f T/m/sec', amplitude, slew_rate)}, 'FontWeight', 'normal');
legend('Trapezoid', 'ADC raster centers', 'Interpreter', 'latex', 'FontSize', FontSize);

fig_filename = sprintf('bipolar_%3.0f_%3.0f_%6.1f_%4.2f_true%d_amp%3.1f', ramp_time * 1e6, total_time * 1e6, adc_duration * 1e6, resolution * 1e3, true_resolution_flag, amplitude);
export_fig(fullfile(output_path, fig_filename), '-r400', '-tif');
close gcf;

%% Calculate the (phi,theta) parameterization of 3D radial k-space trajectories
%--------------------------------------------------------------------------
% Calculate a spiral phyllotaxis pattern (half-spokes, Delacoste 2018 MRM)
%--------------------------------------------------------------------------
if strcmp(trajectory, 'Phyllotaxis')
    [phi, theta] = calculate_spiral_phyllotaxis_pattern(nr_readouts, nr_interleaves, 0);

%--------------------------------------------------------------------------
% Calculate a wobbling Archimedean spiral pole (WASP) pattern
%--------------------------------------------------------------------------
elseif strcmp(trajectory, 'WASP')
    rpa_range = [0 10]; % [degree]
    [phi, theta] = calculate_wasp_pattern(nr_readouts, nr_interleaves, rpa_range);
end

%% Calculate rotated bipolar gradient waveforms in the GCS [Hz/m] ([PE,RO,SL] = [y,x,z] in Pulseq)
%--------------------------------------------------------------------------
% theta: polar angle     (angle from RO to SL)
% phi  : azimuthal angle (angle from SL to PE)
%
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
g_gcs_seq = zeros(3, grad_samples, nr_readouts, nr_interleaves, 'double');
g_gcs_adc = zeros(3, adc_samples , nr_readouts, nr_interleaves, 'double');

for seg = 1:nr_interleaves % equivalent to "segment" in ISMRMRD format
    for lin = 1:nr_readouts

        %------------------------------------------------------------------
        % Calculate a linear index
        %------------------------------------------------------------------
        i = lin + (seg - 1) * nr_readouts;

        %------------------------------------------------------------------
        % Calculate cosine and sine values
        %------------------------------------------------------------------
        cos_theta = cos(theta(i));
        sin_theta = sin(theta(i));
        cos_phi   = cos(phi(i));
        sin_phi   = sin(phi(i));

        %------------------------------------------------------------------
        % Calculate the rotation matrix about the PE direction by theta
        %------------------------------------------------------------------
        R_PE = [1        0             0     ;
                0    cos_theta    -sin_theta ;
                0    sin_theta     cos_theta];

        %------------------------------------------------------------------
        % Calculate the rotation matrix about the RO direction by phi
        %------------------------------------------------------------------
        R_RO = [ cos_phi    0    sin_phi ;
                    0       1       0    ;
                -sin_phi    0    cos_phi];

        %------------------------------------------------------------------
        % Calculate rotated bipolar waveforms
        %------------------------------------------------------------------
        g_gcs_seq(:,:,lin,seg) = R_RO * R_PE * g_base_seq;
        g_gcs_adc(:,:,lin,seg) = R_RO * R_PE * g_base_adc;
    end
end

%% Flip the signs of gradients in the PE and SL directions [PE,RO,SL]
%--------------------------------------------------------------------------
% The "old/compat” option maps Pulseq logical X, Y, and Z axes to the three 
% axes of Siemens logical coordinate system (PE-RO-SL) and uses Siemens native 
% transformation from the logical gradient waveforms to the physical gradient
% waveforms. This also involves scaling all gradient axes by -1 followed by 
% additional scaling on the readout direction by -1. To counter these scaling 
% operations performed in the Pulseq interpreter, our code prepares a Pulseq 
% file by intentionally scaling gradient waveforms along the PE and SL 
% directions by -1.
%--------------------------------------------------------------------------
if strcmp(orientation_mapping, 'old/compat')
    PE_sign = -1;
    RO_sign = +1;
    SL_sign = -1;
else
    PE_sign = +1;
    RO_sign = +1;
    SL_sign = +1;
end

%% Calculate a nominal k-space trajectory of self-navigation in the GCS [rad/m] [PE,RO,SL]
if self_navigation_flag
    % [Hz/m] * [sec] * [2*pi rad/cycle] => * 2*pi [rad/m]
    k_base_adc = cumsum(g_base_adc * (2 * pi) * sys.gradRasterTime, 2); % 3 x adc_samples [rad/m]
end

%% Calculate nominal k-space trajectories in the GCS [rad/m] [PE,RO,SL]
% [Hz/m] * [sec] * [2*pi rad/cycle] => * 2*pi [rad/m]
k_gcs_adc  = cumsum(g_gcs_adc * (2 * pi) * sys.gradRasterTime, 2); % 3 x adc_samples x nr_radial_views [rad/m]

%% Set a slice normal vector
%--------------------------------------------------------------------------
% dNormalSag: Sagittal component of a slice normal vector (in PCS)
% dNormalCor: Coronal component of a slice normal vector (in PCS)
% dNormalTra: Transverse component of a slice normal vector (in PCS)
% dRotAngle: Slice rotation angle ("swap Fre/Pha")
%--------------------------------------------------------------------------
image_ori = 'sagittal';

if ~(exist('dNormalSag', 'var') && exist('dNormalCor', 'var') && exist('dNormalTra', 'var') && exist('dRotAngle', 'var'))
    switch image_ori
        case 'sagittal'
            dNormalSag = 1;
            dNormalCor = 0;
            dNormalTra = 0;
            dRotAngle  = 0; % left-hand rule about the slice normal vector [rad]
        case 'coronal'
            dNormalSag = 0;
            dNormalCor = 1;
            dNormalTra = 0;
            dRotAngle  = 0; % left-hand rule about the slice normal vector [rad]
        case 'axial'
            dNormalSag = 0;
            dNormalCor = 0;
            dNormalTra = 1;
            dRotAngle  = 0; % left-hand rule about the slice normal vector [rad]
        otherwise
    end
end

%% Calculate a rotation matrix from the GCS to the PCS
[R_gcs2pcs, phase_sign, read_sign, main_orientation] = siemens_calculate_transform_gcs_to_pcs(dNormalSag, dNormalCor, dNormalTra, dRotAngle);

%% Calculate a rotation matrix from the PCS to the DCS
R_pcs2dcs = siemens_calculate_transform_pcs_to_dcs('HFS');

%% Calculate a rotation matrix from the GCS to the DCS
R_gcs2dcs = R_pcs2dcs * R_gcs2pcs;

%% Calculate nominal gradient waveforms in the DCS [Hz/m] [x,y,z]
g_dcs_adc = zeros(3, adc_samples, nr_readouts, nr_interleaves, 'double');
for seg = 1:nr_interleaves
    for lin = 1:nr_readouts
        g_dcs_adc(:,:,lin,seg) = R_gcs2dcs * g_gcs_adc(:,:,lin,seg); % 3 x adc_samples
    end
end

%% Calculate the time courses of phase coefficients (adc_samples x Nl x nr_radial_views) [rad/m], [rad/m^2], [rad/m^3]
Nl = 19;
gamma = 4257.59 * (1e4 * 2 * pi); % gyromagnetic ratio for 1H [rad/sec/T]
tstart = tic; fprintf('%s: Calculating the time courses of phase coefficients... ', datetime);
% [Hz/m] / [Hz/T] * [1e3mT/T] => *1e3 [mT/m]
k = calculate_concomitant_field_coefficients(reshape(g_dcs_adc(1,:) / sys.gamma * 1e3, [adc_samples nr_radial_views]), reshape(g_dcs_adc(2,:) / sys.gamma * 1e3, [adc_samples nr_radial_views]), reshape(g_dcs_adc(3,:) / sys.gamma * 1e3, [adc_samples nr_radial_views]), Nl, B0, gamma, sys.gradRasterTime);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%% Calculate the average of k
k_bar = mean(k,3); % adc_samples x Nl

%% Display the phase coefficients
ell_list = [1:6, 8:9].';
figure('Color', 'w', 'Position', [6 178 1749 800]);
for i = 1:length(ell_list)
    ell = ell_list(i);
    ax = subplot(2,4,i);
    if ell == 1
        title_text = '$$k_{1,i}(t) = \int_{0}^{t} G_{x,i}(\tau) d \tau$$';
        ylabel_text = '[rad/m]';
        basis_text = '$$p_1(\mathbf{r}) = x$$';
    elseif ell == 2
        title_text = '$$k_{2,i}(t) = \int_{0}^{t} G_{y,i}(\tau) d \tau$$';
        ylabel_text = '[rad/m]';
        basis_text = '$$p_2(\mathbf{r}) = y$$';
    elseif ell == 3
        title_text = '$$k_{3,i}(t) = \int_{0}^{t} G_{z,i}(\tau) d \tau$$';
        ylabel_text = '[rad/m]';
        basis_text = '$$p_3(\mathbf{r}) = z$$';
    elseif ell == 4
        title_text = '$$k_{4,i}(t) = \int_{0}^{t} G_{z,i}^2(\tau) / (8 B_0) d \tau$$';
        ylabel_text = '[rad/m$^2$]';
        basis_text = '$$p_4(\mathbf{r}) = x^2$$';
    elseif ell == 5
        title_text = '$$k_{5,i}(t) = \int_{0}^{t} G_{z,i}^2(\tau) / (8 B_0) d \tau$$';
        ylabel_text = '[rad/m$^2$]';
        basis_text = '$$p_5(\mathbf{r}) = y^2$$';
    elseif ell == 6
        title_text = '$$k_{6,i}(t) = \int_{0}^{t} (G_{x,i}^2(\tau) + G_{y,i}^2(\tau)) / (2 B_0) d \tau$$';
        ylabel_text = '[rad/m$^2$]';
        basis_text = '$$p_6(\mathbf{r}) = z^2$$';
    elseif ell == 8
        title_text = '$$k_{8,i}(t) = \int_{0}^{t} -G_{y,i}(\tau) G_{z,i}(\tau) / (2 B_0) d \tau$$';
        ylabel_text = '[rad/m$^2$]';
        basis_text = '$$p_8(\mathbf{r}) = yz$$';
    elseif ell == 9
        title_text = '$$k_{9,i}(t) = \int_{0}^{t} -G_{x,i}(\tau) G_{z,i}(\tau) / (2 B_0) d \tau$$';
        ylabel_text = '[rad/m$^2$]';
        basis_text = '$$p_9(\mathbf{r}) = xz$$';
    end

    hold on;
    plot(t_adc * 1e3, reshape(k(:,ell,1:nr_readouts), [adc_samples nr_readouts]));
    if ell > 3
        plot(t_adc * 1e3, k_bar(:,ell), 'Color', 'r', 'LineWidth', 2);
    end
    YLim = ax.YLim;
    text(0.3, YLim(2) * 0.99, basis_text, 'Interpreter', 'latex', 'FontSize', 14, 'VerticalAlignment', 'top');
    xlabel('Time [msec]', 'Interpreter', 'latex', 'FontSize', 12);
    ylabel(ylabel_text, 'Interpreter', 'latex', 'FontSize', 12);
    title(title_text, 'Interpreter', 'latex', 'FontSize', 12);
    set(gca, 'Box', 'On');
    grid on; grid minor;
end
drawnow;

fig_filename = sprintf('concomitant_fields_%3.0f_%3.0f_%6.1f_%4.2f_true%d', ramp_time * 1e6, total_time * 1e6, adc_duration * 1e6, resolution * 1e3, true_resolution_flag);
export_fig(fullfile(output_path, fig_filename), '-r400', '-tif');
close gcf;

%% Display nominal k-space trajectories [rad/m]
if 0
krmax = max(k_gcs_adc(:));

FontSize = 14;
az = 35.31;
el = 32;

figure('Color', 'w', 'Position', [30 410 1210 568]);
ColorOrder = get(gca, 'ColorOrder');

readout_range = (1:nr_readouts).';
ax1 = subplot(1,2,1);
hold on;
grid on; grid minor;
view(az,el);
for i = 1:nr_readouts
    readout_number = readout_range(i);
    hPlot1 = plot3(ax1, k_gcs_adc(1,:,readout_number,1), k_gcs_adc(2,:,readout_number,1), k_gcs_adc(3,:,readout_number,1), '.', 'Color', ColorOrder(2,:));
    xlim([-krmax krmax]);
    ylim([-krmax krmax]);
    zlim([-krmax krmax]);
    axis square;
    xlabel('$k_{\mathrm{PE}}$ [rad/m]', 'Interpreter', 'latex', 'FontSize', FontSize);
    ylabel('$k_{\mathrm{RO}}$ [rad/m]', 'Interpreter', 'latex', 'FontSize', FontSize);
    zlabel('$k_{\mathrm{SL}}$ [rad/m]', 'Interpreter', 'latex', 'FontSize', FontSize);
    pause(0.1);
    set(hPlot1, 'Color', [0.5 0.5 0.5 0.1], 'MarkerSize', 2);
end
%%
ax2 = subplot(1,2,2);
hold on;
grid on; grid minor;
view(az,el);
for i = 1:nr_readouts
    readout_number = readout_range(i);
    hPlot1 = plot3(ax2, k_gcs_adc(1,:,readout_number,2), k_gcs_adc(2,:,readout_number,2), k_gcs_adc(3,:,readout_number,2), '.', 'Color', ColorOrder(2,:));
    xlim([-krmax krmax]);
    ylim([-krmax krmax]);
    zlim([-krmax krmax]);
    axis square;
    xlabel('$k_{\mathrm{PE}}$ [rad/m]', 'Interpreter', 'latex', 'FontSize', FontSize);
    ylabel('$k_{\mathrm{RO}}$ [rad/m]', 'Interpreter', 'latex', 'FontSize', FontSize);
    zlabel('$k_{\mathrm{SL}}$ [rad/m]', 'Interpreter', 'latex', 'FontSize', FontSize);
    pause(0.1);
    set(hPlot1, 'Color', [0.5 0.5 0.5 0.1], 'MarkerSize', 2);
end
end

%% Create an ADC event
adc_delay = delay_pre - adc_duration_pre;
adc = mr.makeAdc(discard_pre + adc_samples + discard_post, 'Dwell', real_dwell_time, 'delay', adc_delay, 'system', sys);

%% Calculate timing (need to decide on the block structure already)
delayTE = round((rf.deadTime + rf_length / 2 + TE1 - delay_pre) / sys.blockDurationRaster) * sys.blockDurationRaster;
delayTR = round((TR - delayTE) / sys.blockDurationRaster) * sys.blockDurationRaster;

%% Create a bipolar gradient event in the RO direction ([PE,RO,SL] = [y,x,z] in Pulseq)
gx_bipolar = g_bipolar;
gx_bipolar.channel = 'x';

%% Create a bipolar gradient event in the PE direction ([PE,RO,SL] = [y,x,z] in Pulseq)
gy_bipolar = g_bipolar;
gy_bipolar.channel = 'y';

%% Create a bipolar gradient event in the SL direction ([PE,RO,SL] = [y,x,z] in Pulseq)
gz_bipolar = g_bipolar;
gz_bipolar.channel = 'z';

%% Create dummy gradient events in the PE, RO, SL directions ([PE,RO,SL] = [y,x,z] in Pulseq)
gx_dummy = gx_bipolar;
gy_dummy = gy_bipolar;
gz_dummy = gz_bipolar;

if self_navigation_flag
    gy_dummy.waveform = PE_sign * g_base_seq(1,:,1); % PE => 'y'
    gx_dummy.waveform = RO_sign * g_base_seq(2,:,1); % RO => 'x'
    gz_dummy.waveform = SL_sign * g_base_seq(3,:,1); % SL => 'z'
else
    gy_dummy.waveform = PE_sign * g_gcs_seq(1,:,1); % PE => 'y'
    gx_dummy.waveform = RO_sign * g_gcs_seq(2,:,1); % RO => 'x'
    gz_dummy.waveform = SL_sign * g_gcs_seq(3,:,1); % SL => 'z'
end

%% Create an alpha-degree hard pulse event (FID navigator)
%--------------------------------------------------------------------------
% rf.delay includes rf.deadTime. That means rf.signal starts at rf.delay 
%--------------------------------------------------------------------------
rf_nav = mr.makeBlockPulse(flip_angle_nav * pi / 180, sys, 'Duration', rf_length_nav);

%% Create an ADC event (FID navigator)
adc_nav = mr.makeAdc(adc_samples_nav, 'Dwell', dwell_time_nav, 'delay', sys.adcDeadTime, 'system', sys);
adc_duration_nav = dwell_time_nav * adc_samples_nav;

%% Calculate timing (FID navigator)
%--------------------------------------------------------------------------
% delayTE0: delayTE0 = round((rf_nav.deadTime + rf_length_nav / 2 + TE_nav - sys.adcDeadTime) / sys.blockDurationRaster) * sys.blockDurationRaster;
% delayTR0: delayTR0 = round((nav_adc_duration + sys.adcDeadTime * 2) / sys.blockDurationRaster) * sys.blockDurationRaster;
%--------------------------------------------------------------------------
delayTE0 = round((rf_nav.deadTime + rf_length_nav / 2 + TE_nav - sys.adcDeadTime) / sys.blockDurationRaster) * sys.blockDurationRaster;
delayTR0 = round((adc_duration_nav + sys.adcDeadTime * 2) / sys.blockDurationRaster) * sys.blockDurationRaster;

%% Pre-register objects that do not change while looping
[~, rf.shapeIDs] = seq.registerRfEvent(rf); % the phase of the RF object will change, therefore we only per-register the shapes

%% Define sequence blocks for an FID navigator
if fid_navigator_flag
    tstart = tic; fprintf('%s: Defining blocks for an FID navigator... ', datetime);

    %----------------------------------------------------------------------
    % Add a new block to the sequence
    %----------------------------------------------------------------------
    seq.addBlock(rf_nav, mr.makeDelay(delayTE0));

    %----------------------------------------------------------------------
    % Add a new block to the sequence
    %----------------------------------------------------------------------
    seq.addBlock(adc_nav, mr.makeDelay(delayTR0));

    %----------------------------------------------------------------------
    % Add a new block to the sequence
    %----------------------------------------------------------------------
    seq.addBlock(mr.makeDelay(T_delay));
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
end

%% Define sequence blocks for an alpha/2-TR/2 sequence
tstart = tic; fprintf('%s: Defining blocks for an alpha/2-TR/2 sequence... ', datetime);

%--------------------------------------------------------------------------
% Set the RF phase increment [rad]
%--------------------------------------------------------------------------
rf_phase_rad = rf_phase * pi / 180;

%--------------------------------------------------------------------------
% Set the counter
%--------------------------------------------------------------------------
count = 1;

%--------------------------------------------------------------------------
% Set the phase of an RF event and an ADC event
%--------------------------------------------------------------------------
rf_half.phaseOffset = rf_phase_rad * mod(count,2);
adc.phaseOffset = rf_phase_rad * mod(count,2);

%--------------------------------------------------------------------------
% Add a new block to the sequence
%--------------------------------------------------------------------------
if ~noise_scan_flag
    seq.addBlock(rf_half, mr.makeDelay(TR / 2));
end

%--------------------------------------------------------------------------
% Increment the counter
%--------------------------------------------------------------------------
count = count + 1;
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%% Define sequence blocks for dummy pulses
for i = 1:rf_dummies
    tstart = tic; fprintf('%s: Defining blocks for dummy pulses (%3d/%3d)... ', datetime, i, rf_dummies);

    %----------------------------------------------------------------------
    % Set the phase of an RF event and an ADC event
    %----------------------------------------------------------------------
    rf.phaseOffset = rf_phase_rad * mod(count,2);
    adc.phaseOffset = rf_phase_rad * mod(count,2);

    %----------------------------------------------------------------------
    % Add a new block to the sequence (block 1)
    %----------------------------------------------------------------------
    if ~noise_scan_flag
        seq.addBlock(rf, mr.makeDelay(delayTE));
    end

    %----------------------------------------------------------------------
    % Add a new block to the sequence (block 2)
    %----------------------------------------------------------------------
    if ~noise_scan_flag
        seq.addBlock(gx_dummy, gy_dummy, gz_dummy, mr.makeDelay(delayTR));
    end

    %----------------------------------------------------------------------
    % Increment the counter
    %----------------------------------------------------------------------
    count = count + 1;
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
end

%% Define sequence blocks for radial dual-echo bSSSP acquisitions
for seg = 1:nr_interleaves

    %% Self-navigation
    if self_navigation_flag
        tstart = tic; fprintf('%s: Defining blocks for self-navigation (SEG=%3d/%3d)... ', datetime, seg, nr_interleaves);

        %------------------------------------------------------------------
        % Set the phase of an RF event and an ADC event
        %------------------------------------------------------------------
        rf.phaseOffset = rf_phase_rad * mod(count,2);
        adc.phaseOffset = rf_phase_rad * mod(count,2);

        %------------------------------------------------------------------
        % Create a bipolar gradient event in the PE direction (PRS)
        % [PE,RO,SL] = [y,x,z] in Pulseq
        %------------------------------------------------------------------
        gy_bipolar.waveform = PE_sign * g_base_seq(1,:); % PE => 'y'

        %------------------------------------------------------------------
        % Create a bipolar gradient event in the RO direction (PRS)
        % [PE,RO,SL] = [y,x,z] in Pulseq
        %------------------------------------------------------------------
        gx_bipolar.waveform = RO_sign * g_base_seq(2,:); % RO => 'x'

        %------------------------------------------------------------------
        % Create a bipolar gradient event in the SL direction (PRS)
        % [PE,RO,SL] = [y,x,z] in Pulseq
        %------------------------------------------------------------------
        gz_bipolar.waveform = SL_sign * g_base_seq(3,:); % SL => 'z'

        %------------------------------------------------------------------
        % Add a new block to the sequence (block 1)
        %------------------------------------------------------------------
        if ~noise_scan_flag
            seq.addBlock(rf, mr.makeDelay(delayTE));
        end

        %------------------------------------------------------------------
        % Add a new block to the sequence (block 2)
        %------------------------------------------------------------------
        if ~noise_scan_flag
            seq.addBlock(gx_bipolar, gy_bipolar, gz_bipolar, adc, mr.makeDelay(delayTR));
        else
            seq.addBlock(adc, mr.makeDelay(delayTR));
        end

        %------------------------------------------------------------------
        % Increment the counter
        %------------------------------------------------------------------
        count = count + 1;
        fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
    end

    %% bSTAR acquisitions (LIN)
    for lin = 1:nr_readouts
        tstart = tic; fprintf('%s: Defining blocks for bSTAR acquisitions (SEG=%3d/%3d)(LIN=%3d/%3d)... ', datetime, seg, nr_interleaves, lin, nr_readouts);

        %------------------------------------------------------------------
        % Set the phase of an RF event and an ADC event
        %------------------------------------------------------------------
        rf.phaseOffset = rf_phase_rad * mod(count,2);
        adc.phaseOffset = rf_phase_rad * mod(count,2);

        %------------------------------------------------------------------
        % Create a bipolar gradient event in the PE direction (PRS)
        % [PE,RO,SL] = [y,x,z] in Pulseq
        %------------------------------------------------------------------
        gy_bipolar.waveform = PE_sign * g_gcs_seq(1,:,lin,seg); % PE => 'y'

        %------------------------------------------------------------------
        % Create a bipolar gradient event in the RO direction (PRS)
        % [PE,RO,SL] = [y,x,z] in Pulseq
        %------------------------------------------------------------------
        gx_bipolar.waveform = RO_sign * g_gcs_seq(2,:,lin,seg); % RO => 'x'

        %------------------------------------------------------------------
        % Create a bipolar gradient event in the SL direction (PRS)
        % [PE,RO,SL] = [y,x,z] in Pulseq
        %------------------------------------------------------------------
        gz_bipolar.waveform = SL_sign * g_gcs_seq(3,:,lin,seg); % SL => 'z'

        %------------------------------------------------------------------
        % Add a new block to the sequence (block 1)
        %------------------------------------------------------------------
        if ~noise_scan_flag
            seq.addBlock(rf, mr.makeDelay(delayTE));
        end

        %------------------------------------------------------------------
        % Add a new block to the sequence (block 2)
        %------------------------------------------------------------------
        if ~noise_scan_flag
            seq.addBlock(gx_bipolar, gy_bipolar, gz_bipolar, adc, mr.makeDelay(delayTR));
        else
            seq.addBlock(adc, mr.makeDelay(delayTR));
        end

        %------------------------------------------------------------------
        % Increment the counter
        %------------------------------------------------------------------
        count = count + 1;
        fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
    end
end

%% Define sequence blocks for an FID navigator
if fid_navigator_flag
    tstart = tic; fprintf('%s: Defining blocks for an FID navigator... ', datetime);

    %----------------------------------------------------------------------
    % Add a new block to the sequence
    %----------------------------------------------------------------------
    seq.addBlock(mr.makeDelay(T_delay));

    %----------------------------------------------------------------------
    % Add a new block to the sequence
    %----------------------------------------------------------------------
    seq.addBlock(rf_nav, mr.makeDelay(delayTE0));

    %----------------------------------------------------------------------
    % Add a new block to the sequence
    %----------------------------------------------------------------------
    seq.addBlock(adc_nav, mr.makeDelay(delayTR0));
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
seq.setDefinition('Bandwidth1', bandwidth1);
seq.setDefinition('Bandwidth2', bandwidth2);
seq.setDefinition('BaseResolution', base_resolution);
seq.setDefinition('Contrasts', nr_echoes);
seq.setDefinition('FOV', [fov_read fov_read fov_read]);
seq.setDefinition('FlipAngle', flip_angle);
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
seq.setDefinition('SelfNavigationFlag', self_navigation_flag);
seq.setDefinition('SlewRate', slew_rate);
seq.setDefinition('SlicePositions', 0);
seq.setDefinition('SliceThickness', fov_read);
seq.setDefinition('SliceGap', 0);
seq.setDefinition('TE1', TE1);
seq.setDefinition('TE2', TE2);
seq.setDefinition('TR', TR);
seq.setDefinition('TotalTime', total_time);
seq.setDefinition('Trajectory', trajectory);
seq.setDefinition('TrueResolutionFlag', true_resolution_flag);

seq.setDefinition('FidNavigatorFlag', fid_navigator_flag);
seq.setDefinition('FlipAngleNav', flip_angle_nav);
seq.setDefinition('RFLengthNav', rf_length_nav);
seq.setDefinition('AdcSamplesNav', adc_samples_nav);
seq.setDefinition('T_delay', T_delay);
seq.setDefinition('TE_nav', TE_nav);

seq_file = fullfile(output_path, sprintf('%s.seq', filename_base));

%% Write a .seq file
tstart = tic; fprintf('%s: Writing a .seq file: %s... ', datetime, seq_file);
seq.write(seq_file); % Write to a pulseq file
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%% Plot sequence and k-space diagrams
seq.plot('timeRange', [0 5]);
%seq.plot('timeRange', [0 5], 'showBlocks', 1);

set(gcf, 'Position', [61 425 1858 553]);
export_fig(fullfile(output_path, filename_base), '-r300', '-tif');
close all;

%% k-space trajectory calculation
if 0
[ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing] = seq.calculateKspacePP();

% plot 3D k-space
figure; plot3(ktraj(1,:), ktraj(2,:), ktraj(3,:), 'b'); % a 3D k-space plot
axis('equal'); % enforce aspect ratio for the correct trajectory display
hold; plot3(ktraj_adc(1,:), ktraj_adc(2,:), ktraj_adc(3,:), 'r.'); % plot the sampling points
title('full k-space trajectory (k_x x k_y x k_z)');
end

%% very optional slow step, but useful for testing during development e.g. for the real TE, TR or for staying within slewrate limits  
if 0
rep = seq.testReport;
fprintf([rep{:}]);
end

%% Write ISMRMRD format for nominal k-space trajectories
traj_filename = sprintf('traj_%s.h5', filename_base);
ismrmrd_traj_file = fullfile(output_path, traj_filename);
dset = ismrmrd.Dataset(ismrmrd_traj_file);

%% Add a block of acquisitions at a time
%--------------------------------------------------------------------------
% ISMRMRD header from ismrmrd.h (Header for each MR acquisition)
%--------------------------------------------------------------------------
% uint16_t version;                                    /**< First unsigned int indicates the version */
% uint64_t flags;                                      /**< bit field with flags */
% uint32_t measurement_uid;                            /**< Unique ID for the measurement */
% uint32_t scan_counter;                               /**< Current acquisition number in the measurement */
% uint32_t acquisition_time_stamp;                     /**< Acquisition clock */
% uint32_t physiology_time_stamp[ISMRMRD_PHYS_STAMPS]; /**< Physiology time stamps, e.g. ecg, breating, etc. */
% uint16_t number_of_samples;                          /**< Number of samples acquired */
% uint16_t available_channels;                         /**< Available coils */
% uint16_t active_channels;                            /**< Active coils on current acquisiton */
% uint64_t channel_mask[ISMRMRD_CHANNEL_MASKS];        /**< Mask to indicate which channels are active. Support for 1024 channels */
% uint16_t discard_pre;                                /**< Samples to be discarded at the beginning of  acquisition */
% uint16_t discard_post;                               /**< Samples to be discarded at the end of acquisition */
% uint16_t center_sample;                              /**< Sample at the center of k-space */
% uint16_t encoding_space_ref;                         /**< Reference to an encoding space, typically only one per acquisition */
% uint16_t trajectory_dimensions;                      /**< Indicates the dimensionality of the trajectory vector (0 means no trajectory) */
% float sample_time_us;                                /**< Time between samples in micro seconds, sampling BW */
% float position[3];                                   /**< Three-dimensional spatial offsets from isocenter */
% float read_dir[3];                                   /**< Directional cosines of the readout/frequency encoding */
% float phase_dir[3];                                  /**< Directional cosines of the phase */
% float slice_dir[3];                                  /**< Directional cosines of the slice direction */
% float patient_table_position[3];                     /**< Patient table off-center */
% ISMRMRD_EncodingCounters idx;                        /**< Encoding loop counters, see above */
% int32_t user_int[ISMRMRD_USER_INTS];                 /**< Free user parameters */
% float user_float[ISMRMRD_USER_FLOATS];               /**< Free user parameters */
%--------------------------------------------------------------------------
% Where EncodingCounters are defined as:
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
nr_acquisitions = nr_radial_views + (nr_interleaves * self_navigation_flag) + (2 * fid_navigator_flag);
acqblock = ismrmrd.Acquisition(nr_acquisitions);

%--------------------------------------------------------------------------
% Set the header elements that don't change
%--------------------------------------------------------------------------
acqblock.head.version(:)               = 1;
acqblock.head.number_of_samples(:)     = adc_samples;
acqblock.head.available_channels(:)    = 1;
acqblock.head.active_channels(:)       = 1;
acqblock.head.discard_pre(:)           = discard_pre;
acqblock.head.discard_post(:)          = discard_post;
acqblock.head.trajectory_dimensions(:) = 3;
acqblock.head.sample_time_us(:)        = real_dwell_time * 1e6; % [usec]

%% Loop over the acquisitions, set the header, set the data and append
scan_counter = 0;

%--------------------------------------------------------------------------
% FID navigator
%--------------------------------------------------------------------------
if fid_navigator_flag

    %----------------------------------------------------------------------
    % scan_counter: Current acquisition number in the measurement
    %----------------------------------------------------------------------
    scan_counter = scan_counter + 1;
    acqblock.head.scan_counter(scan_counter) = scan_counter;

    %----------------------------------------------------------------------
    % flags: bit field with flags
    %----------------------------------------------------------------------
    acqblock.head.flagClearAll(scan_counter);
    acqblock.head.flagSet('ACQ_IS_NAVIGATION_DATA', scan_counter);

    %----------------------------------------------------------------------
    % number_of_samples: Number of samples acquired
    %----------------------------------------------------------------------
    acqblock.head.number_of_samples(scan_counter) = adc_samples_nav;

    %----------------------------------------------------------------------
    % discard_pre: Samples to be discarded at the beginning of acquisition
    %----------------------------------------------------------------------
    acqblock.head.discard_pre(scan_counter) = 0;

    %----------------------------------------------------------------------
    % discard_pre: Samples to be discarded at the end of acquisition
    %----------------------------------------------------------------------
    acqblock.head.discard_post(scan_counter) = 0;

    %----------------------------------------------------------------------
    % trajectory_dimensions: Indicates the dimensionality of the
    % trajectory vector (0 means no trajectory)
    %----------------------------------------------------------------------
    acqblock.head.trajectory_dimensions(scan_counter) = 0;

    %----------------------------------------------------------------------
    % discard_pre: Samples to be discarded at the end of acquisition
    %----------------------------------------------------------------------
    acqblock.head.sample_time_us(scan_counter) = dwell_time_nav * 1e6;

    %----------------------------------------------------------------------
    % average: signal average number
    %----------------------------------------------------------------------
    acqblock.head.idx.average(scan_counter) = 0;

    %----------------------------------------------------------------------
    % contrast: echo number in multi-echo
    %----------------------------------------------------------------------
    acqblock.head.idx.contrast(scan_counter) = 0;

    %----------------------------------------------------------------------
    % segment: segment number for segmented acquisition
    %----------------------------------------------------------------------
    acqblock.head.idx.segment(scan_counter) = 0;

    %----------------------------------------------------------------------
    % traj
    %----------------------------------------------------------------------
    acqblock.traj{scan_counter} = [];

    %----------------------------------------------------------------------
    % data
    %----------------------------------------------------------------------
    acqblock.data{scan_counter} = zeros(adc_samples_nav, 1, 'single');
end

%--------------------------------------------------------------------------
% Interleaves (SEG)
%--------------------------------------------------------------------------
for seg = 1:nr_interleaves

    %----------------------------------------------------------------------
    % Self-navigation
    %----------------------------------------------------------------------
    if self_navigation_flag

        %------------------------------------------------------------------
        % scan_counter: Current acquisition number in the measurement
        %------------------------------------------------------------------
        scan_counter = scan_counter + 1;
        acqblock.head.scan_counter(scan_counter) = scan_counter;

        %------------------------------------------------------------------
        % flags: bit field with flags
        %------------------------------------------------------------------
        acqblock.head.flagClearAll(scan_counter);
        acqblock.head.flagSet('ACQ_USER1', scan_counter);

        %------------------------------------------------------------------
        % kspace_encode_step_1: phase encoding line number
        %------------------------------------------------------------------
        acqblock.head.idx.kspace_encode_step_1(scan_counter) = 0; % starts from 0

        %------------------------------------------------------------------
        % kspace_encode_step_2: partition encoding number
        %------------------------------------------------------------------
        acqblock.head.idx.kspace_encode_step_2(scan_counter) = 0; % starts from 0

        %------------------------------------------------------------------
        % average: signal average number
        %------------------------------------------------------------------
        acqblock.head.idx.average(scan_counter) = 0;

        %------------------------------------------------------------------
        % contrast: echo number in multi-echo
        %------------------------------------------------------------------
        acqblock.head.idx.contrast(scan_counter) = 0;

        %------------------------------------------------------------------
        % segment: segment number for segmented acquisition
        %------------------------------------------------------------------
        acqblock.head.idx.segment(scan_counter) = seg - 1;

        %------------------------------------------------------------------
        % traj
        %------------------------------------------------------------------
        acqblock.traj{scan_counter} = single(k_base_adc); % [rad/m]

        %------------------------------------------------------------------
        % data
        %------------------------------------------------------------------
        acqblock.data{scan_counter} = zeros(adc_samples, 1, 'single');
        fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
    end

    %----------------------------------------------------------------------
    % bSTAR acquisitions (LIN)
    %----------------------------------------------------------------------
    for lin = 1:nr_readouts
        tstart = tic; fprintf('%s: Filling an acquisition block (SEG=%3d/%3d)(LIN=%3d/%3d)... ', datetime, seg, nr_interleaves, lin, nr_readouts);

        %------------------------------------------------------------------
        % scan_counter: Current acquisition number in the measurement
        %------------------------------------------------------------------
        scan_counter = scan_counter + 1;
        acqblock.head.scan_counter(scan_counter) = scan_counter;

        %------------------------------------------------------------------
        % kspace_encode_step_1: phase encoding line number
        %------------------------------------------------------------------
        acqblock.head.idx.kspace_encode_step_1(scan_counter) = lin - 1; % starts from 0

        %------------------------------------------------------------------
        % kspace_encode_step_2: partition encoding number
        %------------------------------------------------------------------
        acqblock.head.idx.kspace_encode_step_2(scan_counter) = 0; % starts from 0

        %------------------------------------------------------------------
        % average: signal average number
        %------------------------------------------------------------------
        acqblock.head.idx.average(scan_counter) = 0;

        %------------------------------------------------------------------
        % contrast: echo number in multi-echo
        %------------------------------------------------------------------
        acqblock.head.idx.contrast(scan_counter) = 0;

        %------------------------------------------------------------------
        % segment: segment number for segmented acquisition
        %------------------------------------------------------------------
        acqblock.head.idx.segment(scan_counter) = seg - 1;

        %------------------------------------------------------------------
        % traj
        %------------------------------------------------------------------
        acqblock.traj{scan_counter} = single(k_gcs_adc(:,:,lin,seg)); % [rad/m]

        %------------------------------------------------------------------
        % data
        %------------------------------------------------------------------
        i = lin + (seg - 1) * nr_readouts;
        data = zeros(adc_samples, 1, 'single');
        data(1) = phi(i);
        data(2) = theta(i);
        acqblock.data{scan_counter} = data;
        fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
    end
end

%--------------------------------------------------------------------------
% FID navigator
%--------------------------------------------------------------------------
if fid_navigator_flag

    %----------------------------------------------------------------------
    % scan_counter: Current acquisition number in the measurement
    %----------------------------------------------------------------------
    scan_counter = scan_counter + 1;
    acqblock.head.scan_counter(scan_counter) = scan_counter;

    %----------------------------------------------------------------------
    % flags: bit field with flags
    %----------------------------------------------------------------------
    acqblock.head.flagClearAll(scan_counter);
    acqblock.head.flagSet('ACQ_IS_NAVIGATION_DATA', scan_counter);

    %----------------------------------------------------------------------
    % number_of_samples: Number of samples acquired
    %----------------------------------------------------------------------
    acqblock.head.number_of_samples(scan_counter) = adc_samples_nav;

    %----------------------------------------------------------------------
    % discard_pre: Samples to be discarded at the beginning of acquisition
    %----------------------------------------------------------------------
    acqblock.head.discard_pre(scan_counter) = 0;

    %----------------------------------------------------------------------
    % discard_pre: Samples to be discarded at the end of acquisition
    %----------------------------------------------------------------------
    acqblock.head.discard_post(scan_counter) = 0;

    %----------------------------------------------------------------------
    % trajectory_dimensions: Indicates the dimensionality of the
    % trajectory vector (0 means no trajectory)
    %----------------------------------------------------------------------
    acqblock.head.trajectory_dimensions(scan_counter) = 0;

    %----------------------------------------------------------------------
    % discard_pre: Samples to be discarded at the end of acquisition
    %----------------------------------------------------------------------
    acqblock.head.sample_time_us(scan_counter) = dwell_time_nav * 1e6;

    %----------------------------------------------------------------------
    % average: signal average number
    %----------------------------------------------------------------------
    acqblock.head.idx.average(scan_counter) = 0;

    %----------------------------------------------------------------------
    % contrast: echo number in multi-echo
    %----------------------------------------------------------------------
    acqblock.head.idx.contrast(scan_counter) = 0;

    %----------------------------------------------------------------------
    % segment: segment number for segmented acquisition
    %----------------------------------------------------------------------
    acqblock.head.idx.segment(scan_counter) = 0;

    %----------------------------------------------------------------------
    % traj
    %----------------------------------------------------------------------
    acqblock.traj{scan_counter} = [];

    %----------------------------------------------------------------------
    % data
    %----------------------------------------------------------------------
    acqblock.data{scan_counter} = zeros(adc_samples_nav, 1, 'single');
end

%--------------------------------------------------------------------------
% Append the acquisition block
%--------------------------------------------------------------------------
dset.appendAcquisition(acqblock);

%% Create an XML header
%--------------------------------------------------------------------------
% We create a matlab struct and then serialize it to xml.
% Look at the xml schema to see what the field names should be
%--------------------------------------------------------------------------
header = [];

%--------------------------------------------------------------------------
% Experimental Conditions (Required)
%--------------------------------------------------------------------------
header.experimentalConditions.H1resonanceFrequency_Hz = B0 * sys.gamma; % [T] * [Hz/T] = > [T]

%--------------------------------------------------------------------------
% Acquisition System Information (Optional)
%--------------------------------------------------------------------------
header.acquisitionSystemInformation.systemVendor = 'Pulseq';
header.acquisitionSystemInformation.systemModel = 'v1.4.1';

%--------------------------------------------------------------------------
% Trajectory
%--------------------------------------------------------------------------
header.encoding.trajectory = 'radial';

%--------------------------------------------------------------------------
% Encoding Space (Nkx, Nky, Nkz)
%--------------------------------------------------------------------------
header.encoding.encodedSpace.fieldOfView_mm.x = fov_read * readout_os_factor * 1e3; % [m] * [1e3mm/m] => *1e3 [mm]
header.encoding.encodedSpace.fieldOfView_mm.y = fov_read * readout_os_factor * 1e3; % [m] * [1e3mm/m] => *1e3 [mm]
header.encoding.encodedSpace.fieldOfView_mm.z = fov_read * readout_os_factor * 1e3; % [m] * [1e3mm/m] => *1e3 [mm]

encoded_matrix_size = ceil(fov_read * readout_os_factor / resolution);

header.encoding.encodedSpace.matrixSize.x = encoded_matrix_size;
header.encoding.encodedSpace.matrixSize.y = encoded_matrix_size;
header.encoding.encodedSpace.matrixSize.z = encoded_matrix_size;

%--------------------------------------------------------------------------
% Recon Space (Nx, Ny, Nz)
%--------------------------------------------------------------------------
header.encoding.reconSpace.fieldOfView_mm.x = fov_read * 1e3; % [m] * [1e3mm/m] => *1e3 [mm]
header.encoding.reconSpace.fieldOfView_mm.y = fov_read * 1e3; % [m] * [1e3mm/m] => *1e3 [mm]
header.encoding.reconSpace.fieldOfView_mm.z = fov_read * 1e3; % [m] * [1e3mm/m] => *1e3 [mm]

recon_matrix_size = ceil(fov_read / resolution / grid_scale_factor);

header.encoding.reconSpace.matrixSize.x = recon_matrix_size;
header.encoding.reconSpace.matrixSize.y = recon_matrix_size;
header.encoding.reconSpace.matrixSize.z = recon_matrix_size;

%--------------------------------------------------------------------------
% Encoding Limits
%--------------------------------------------------------------------------
header.encoding.encodingLimits.kspace_encoding_step_1.minimum = 0;
header.encoding.encodingLimits.kspace_encoding_step_1.maximum = nr_readouts - 1;
header.encoding.encodingLimits.kspace_encoding_step_1.center = 0;

header.encoding.encodingLimits.average.minimum = 0;
header.encoding.encodingLimits.average.maximum = 0;
header.encoding.encodingLimits.average.center = 0;

header.encoding.encodingLimits.contrast.minimum = 0;
header.encoding.encodingLimits.contrast.maximum = 0;
header.encoding.encodingLimits.contrast.center = 0;

header.encoding.encodingLimits.segment.minimum = 0;
header.encoding.encodingLimits.segment.maximum = nr_interleaves - 1;
header.encoding.encodingLimits.segment.center = 0;

%--------------------------------------------------------------------------
% Trajectory Description (Optional)
%--------------------------------------------------------------------------
header.encoding.trajectoryDescription.identifier = '';

% Long
idx = 1;
header.encoding.trajectoryDescription.userParameterLong(idx).value = true_resolution_flag;
header.encoding.trajectoryDescription.userParameterLong(idx).name = 'true_resolution_flag';

idx = idx + 1;
header.encoding.trajectoryDescription.userParameterLong(idx).value = fid_navigator_flag;
header.encoding.trajectoryDescription.userParameterLong(idx).name = 'fid_navigator_flag';

idx = idx + 1;
header.encoding.trajectoryDescription.userParameterLong(idx).value = noise_scan_flag;
header.encoding.trajectoryDescription.userParameterLong(idx).name = 'noise_scan_flag';

idx = idx + 1;
header.encoding.trajectoryDescription.userParameterLong(idx).value = self_navigation_flag;
header.encoding.trajectoryDescription.userParameterLong(idx).name = 'self_navigation_flag';

idx = idx + 1;
header.encoding.trajectoryDescription.userParameterLong(idx).value = grad_samples;
header.encoding.trajectoryDescription.userParameterLong(idx).name = 'grad_samples';

idx = idx + 1;
header.encoding.trajectoryDescription.userParameterLong(idx).value = adc_samples;
header.encoding.trajectoryDescription.userParameterLong(idx).name = 'adc_samples';

idx = idx + 1;
header.encoding.trajectoryDescription.userParameterLong(idx).value = ramp_time;
header.encoding.trajectoryDescription.userParameterLong(idx).name = 'ramp_time';

idx = idx + 1;
header.encoding.trajectoryDescription.userParameterLong(idx).value = total_time;
header.encoding.trajectoryDescription.userParameterLong(idx).name = 'total_time';

% Double
idx = 1;
header.encoding.trajectoryDescription.userParameterDouble(idx).value = grid_scale_factor;
header.encoding.trajectoryDescription.userParameterDouble(idx).name = 'grid_scale_factor';

idx = idx + 1;
header.encoding.trajectoryDescription.userParameterDouble(idx).value = sys.gradRasterTime;
header.encoding.trajectoryDescription.userParameterDouble(idx).name = 'grad_raster_time';

idx = idx + 1;
header.encoding.trajectoryDescription.userParameterDouble(idx).value = real_dwell_time;
header.encoding.trajectoryDescription.userParameterDouble(idx).name = 'real_dwell_time';

idx = idx + 1;
header.encoding.trajectoryDescription.userParameterDouble(idx).value = TE_nav;
header.encoding.trajectoryDescription.userParameterDouble(idx).name = 'TE_nav';

idx = idx + 1;
header.encoding.trajectoryDescription.userParameterDouble(idx).value = dwell_time_nav;
header.encoding.trajectoryDescription.userParameterDouble(idx).name = 'dwell_time_nav';

%% Serialize and write to the data set
xmlstring = ismrmrd.xml.serialize(header);
dset.writexml(xmlstring);

%% Write the dataset
dset.close();

%% Read an ISMRMRD file
tstart = tic; fprintf('%s: Reading an ISMRMRD file: %s... ', datetime, ismrmrd_traj_file);
if exist(ismrmrd_traj_file, 'file')
    traj_dset = ismrmrd.Dataset(ismrmrd_traj_file, 'dataset');
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
else
    error('File %s does not exist.  Please generate it.' , ismrmrd_traj_file);
end

%% Get imaging parameters from an XML header
traj_header = ismrmrd.xml.deserialize(traj_dset.readxml);

%--------------------------------------------------------------------------
% Encoding Space (Nkx, Nky, Nkz)
%--------------------------------------------------------------------------
encoded_fov(1) = traj_header.encoding.encodedSpace.fieldOfView_mm.x * 1e-3; % [m] RO
encoded_fov(2) = traj_header.encoding.encodedSpace.fieldOfView_mm.y * 1e-3; % [m] PE
encoded_fov(3) = traj_header.encoding.encodedSpace.fieldOfView_mm.z * 1e-3; % [m] SL

Nkx = traj_header.encoding.encodedSpace.matrixSize.x; % number of readout samples in k-space
Nky = traj_header.encoding.encodedSpace.matrixSize.y; % number of phase-encoding steps in k-space
Nkz = traj_header.encoding.encodedSpace.matrixSize.z; % number of slice-encoding steps in k-space

encoded_resolution = encoded_fov ./ [Nkx Nky Nkz]; % [m]

%--------------------------------------------------------------------------
% Recon Space (Nx, Ny, Nz)
%--------------------------------------------------------------------------
recon_fov(1) = traj_header.encoding.reconSpace.fieldOfView_mm.x * 1e-3; % [m] RO
recon_fov(2) = traj_header.encoding.reconSpace.fieldOfView_mm.y * 1e-3; % [m] PE
recon_fov(3) = traj_header.encoding.reconSpace.fieldOfView_mm.z * 1e-3; % [m] SL

Nx = traj_header.encoding.reconSpace.matrixSize.x; % number of samples in image space (RO)
Ny = traj_header.encoding.reconSpace.matrixSize.y; % number of samples in image space (PE)
Nz = traj_header.encoding.reconSpace.matrixSize.z; % number of samples in image space (SL)

recon_resolution = recon_fov ./ [Nx Ny Nz]; % [m]

%% Read acquisitions
tstart = tic; fprintf('%s: Reading acquisitions... ', datetime);
raw_traj = traj_dset.readAcquisition(); % read all the acquisitions
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

acq_is_noise_measurement = raw_traj.head.flagIsSet('ACQ_IS_NOISE_MEASUREMENT');
acq_is_navigation_data   = raw_traj.head.flagIsSet('ACQ_IS_NAVIGATION_DATA'); % FID navigator
acq_user1                = raw_traj.head.flagIsSet('ACQ_USER1'); % self-navigation

nav_traj  = raw_traj.select(find(acq_is_navigation_data));
self_traj = raw_traj.select(find(acq_user1));
img_traj = raw_traj.select(find(~(acq_is_navigation_data | acq_user1)));

figure('Color', 'w', 'Position',  [3 557 1237 421]);
subplot(1,2,1);
stem(acq_is_navigation_data);
title('ACQ\_IS\_NAVIGATION\_DATA');
grid on; grid minor;

subplot(1,2,2);
stem(acq_user1);
title('ACQ\_USER1');
grid on; grid minor;

%% Parse an ISMRMRD header
nr_phase_encoding_steps = double(max(img_traj.head.idx.kspace_encode_step_1)) + 1; % nr_interleaves for spiral imaging
nr_slice_encoding_steps = double(max(img_traj.head.idx.kspace_encode_step_2)) + 1;
nr_averages             = double(max(img_traj.head.idx.average)) + 1;
nr_slices               = double(max(img_traj.head.idx.slice)) + 1;
nr_contrasts            = double(max(img_traj.head.idx.contrast)) + 1;
nr_phases               = double(max(img_traj.head.idx.phase)) + 1;
nr_repetitions          = double(max(img_traj.head.idx.repetition)) + 1;
nr_sets                 = double(max(img_traj.head.idx.set)) + 1;
nr_segments             = double(max(img_traj.head.idx.segment)) + 1;
center_sample           = double(max(img_traj.head.center_sample));
nr_channels             = double(max(img_traj.head.active_channels));

%--------------------------------------------------------------------------
% Get the dimensionality of the trajectory vector (0 means no trajectory)
%--------------------------------------------------------------------------
trajectory_dimensions = double(max(img_traj.head.trajectory_dimensions));

%--------------------------------------------------------------------------
% grad_samples 
%--------------------------------------------------------------------------
grad_samples = traj_header.encoding.trajectoryDescription.userParameterLong(5).value;

%--------------------------------------------------------------------------
% adc_samples
%--------------------------------------------------------------------------
adc_samples = traj_header.encoding.trajectoryDescription.userParameterLong(6).value;

%--------------------------------------------------------------------------
% Get the gradient raster time in [sec]
%--------------------------------------------------------------------------
grad_raster_time = traj_header.encoding.trajectoryDescription.userParameterDouble(2).value; % [sec]

%--------------------------------------------------------------------------
% Get the real dwell time in [sec]
%--------------------------------------------------------------------------
real_dwell_time = traj_header.encoding.trajectoryDescription.userParameterDouble(3).value; % [sec]

%--------------------------------------------------------------------------
% Calculate the readout time, tau [sec]
%--------------------------------------------------------------------------
tau = adc_samples * real_dwell_time; % readout time [sec]

%% Display an ISMRMRD header
fprintf('========================= ISMRMRD header ========================\n');
fprintf('encoded_fov [mm]   = %8.4f %8.4f %8.4f\n', encoded_fov(1) * 1e3, encoded_fov(2) * 1e3, encoded_fov(3) * 1e3);
fprintf('Nkx Nky Nkz        = %d      %d      %d\n', Nkx, Nky, Nkz);
fprintf('encoded_resolution = %8.4f %8.4f %8.4f\n', encoded_resolution(1) * 1e3, encoded_resolution(2) * 1e3, encoded_resolution(3) * 1e3);
fprintf('-----------------------------------------------------------------\n');
fprintf('recon_fov [mm]     = %8.4f %8.4f %8.4f\n', recon_fov(1) * 1e3, recon_fov(2) * 1e3, recon_fov(3) * 1e3);
fprintf('Nx Ny Nz           = %d      %d        %d\n', Nx, Ny, Nz);
fprintf('recon_resolution   = %8.4f %8.4f %8.4f\n', recon_resolution(1) * 1e3, recon_resolution(2) * 1e3, recon_resolution(3) * 1e3);
fprintf('-----------------------------------------------------------------\n');
fprintf('trajectory              = %s\n', traj_header.encoding.trajectory);
fprintf('grad_samples            = %d\n', grad_samples);
fprintf('adc_samples             = %d\n', adc_samples);
fprintf('discard_pre             = %d\n', discard_pre);
fprintf('discard_post            = %d\n', discard_post);
fprintf('center_sample           = %d\n', center_sample);
fprintf('nr_channels             = %d\n', nr_channels);
fprintf('nr_phase_encoding_steps = %d\n', nr_phase_encoding_steps);
fprintf('nr_slice_encoding_steps = %d\n', nr_slice_encoding_steps);
fprintf('nr_averages             = %d\n', nr_averages);
fprintf('nr_slices               = %d\n', nr_slices);
fprintf('nr_contrasts            = %d\n', nr_contrasts);
fprintf('nr_phases               = %d\n', nr_phases);
fprintf('nr_repetitions          = %d\n', nr_repetitions);
fprintf('nr_sets                 = %d\n', nr_sets);
fprintf('nr_segments             = %d\n', nr_segments);
fprintf('grad_raster_time        = %5.2f [usec]\n', grad_raster_time * 1e6);
fprintf('real_dwell_time         = %5.2f [usec]\n', real_dwell_time * 1e6);
fprintf('readout_time (tau)      = %5.2f [msec]\n', tau * 1e3);
fprintf('=================================================================\n');

end

%%
return
phi2 = zeros(nr_radial_views, 1, 'single');
theta2 = zeros(nr_radial_views, 1, 'single');

for i = 1:nr_radial_views
    phi2(i) = img_traj.data{i}(1);
    theta2(i) = img_traj.data{i}(2);
end

[phi phi2 phi - phi2]
