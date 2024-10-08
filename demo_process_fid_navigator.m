% demo_process_fid_navigator.m
% Written by Nam Gyun Lee
% Email: namgyunl@usc.edu, ggang56@gmail.com (preferred)
% Started: 02/27/2024, Last modified: 09/17/2024

%% Clean slate
close all; clearvars -except json_files json_number json_file nr_jsob_files; clc;

%% Start a stopwatch timer
start_time = tic;

%% Read a .json file
tstart = tic; fprintf('%s: Reading a .json file: %s... ', datetime, json_file);
fid = fopen(json_file); 
json_txt = fread(fid, [1 inf], 'char=>char'); 
fclose(fid);
json = jsondecode(json_txt);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%--------------------------------------------------------------------------
% Define the full path of a filename
%--------------------------------------------------------------------------
if ispc
    siemens_twix_file  = strrep(json.siemens_twix_file, '/', '\');
    ismrmrd_data_file  = strrep(json.ismrmrd_data_file, '/', '\');
    ismrmrd_noise_file = strrep(json.ismrmrd_noise_file, '/', '\');
    ismrmrd_traj_file  = strrep(json.ismrmrd_traj_file, '/', '\');
    output_path        = strrep(json.output_path, '/', '\');
else
    siemens_twix_file  = json.siemens_twix_file;
    ismrmrd_data_file  = json.ismrmrd_data_file;
    ismrmrd_noise_file = json.ismrmrd_noise_file;
    ismrmrd_traj_file  = json.ismrmrd_traj_file;
    output_path        = json.output_path;
end

%--------------------------------------------------------------------------
% Define the BART path
%--------------------------------------------------------------------------
bart_path = json.bart_path;

%--------------------------------------------------------------------------
% Reconstruction parameters
%--------------------------------------------------------------------------
recon_matrix_size = json.recon_conf.recon_matrix_size.'; % reconstruction matrix size

%% Make an output path
mkdir(output_path);

%% Set up BART commands
%--------------------------------------------------------------------------
% Define a BART command
%--------------------------------------------------------------------------
if ispc
    command_prefix = 'wsl';
else
    command_prefix = '';
end
bart_command = sprintf('%s %s/bart', command_prefix, bart_path);

%--------------------------------------------------------------------------
% Translate from a Windows path to a WSL path 
%--------------------------------------------------------------------------
if ispc
    bart_output_path = strrep(output_path, '\', '/');
    bart_output_path = sprintf('/mnt/%s/%s/', lower(bart_output_path(1)), bart_output_path(4:end));
else
    bart_output_path = sprintf('%s/', output_path);
end

%% Read an ISMRMRD file
%--------------------------------------------------------------------------
% k-space trajectories
%--------------------------------------------------------------------------
tstart = tic; fprintf('%s: Reading an ISMRMRD file: %s... ', datetime, ismrmrd_traj_file);
if exist(ismrmrd_traj_file, 'file')
    traj_dset = ismrmrd.Dataset(ismrmrd_traj_file, 'dataset');
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
else
    error('File %s does not exist.  Please generate it.' , ismrmrd_traj_file);
end

%--------------------------------------------------------------------------
% k-space data
%--------------------------------------------------------------------------
tstart = tic; fprintf('%s: Reading an ISMRMRD file: %s... ', datetime, ismrmrd_data_file);
if exist(ismrmrd_data_file, 'file')
    data_dset = ismrmrd.Dataset(ismrmrd_data_file, 'dataset');
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
else
    error('File %s does not exist.  Please generate it.' , ismrmrd_data_file);
end

%% Get imaging parameters from an XML header
traj_header = ismrmrd.xml.deserialize(traj_dset.readxml);
data_header = ismrmrd.xml.deserialize(data_dset.readxml);

%--------------------------------------------------------------------------
% measurement information
%--------------------------------------------------------------------------
patient_position        = data_header.measurementInformation.patientPosition;
relative_table_position = data_header.measurementInformation.relativeTablePosition;

%--------------------------------------------------------------------------
% Sequence parameters
%--------------------------------------------------------------------------
TR         = data_header.sequenceParameters.TR * 1e-3;     % [msec] * [sec/1e3msec] => [sec]
TE         = data_header.sequenceParameters.TE * 1e-3;     % [msec] * [sec/1e3msec] => [sec]
flip_angle = data_header.sequenceParameters.flipAngle_deg; % [degrees]

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

%--------------------------------------------------------------------------
% Number of receive channels
%--------------------------------------------------------------------------
nr_channels = data_header.acquisitionSystemInformation.receiverChannels;

%% Read acquisitions (traj)
tstart = tic; fprintf('%s: Reading acquisitions (traj)... ', datetime);
raw_traj = traj_dset.readAcquisition(); % read all the acquisitions
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%--------------------------------------------------------------------------
% Get data type
%--------------------------------------------------------------------------
acq_is_navigation_data = raw_traj.head.flagIsSet('ACQ_IS_NAVIGATION_DATA'); % FID navigator

figure('Color', 'w');
stem(acq_is_navigation_data);
title('ACQ\_IS\_NAVIGATION\_DATA (FID navigator)', 'Interpreter', 'latex');
grid on;

%% Read acquisitions (FID navigator)
nav_list = find(acq_is_navigation_data);
nr_navigators = length(nav_list);

nav_data = cell(nr_navigators,1);
for idx = 1:nr_navigators
    tstart = tic; fprintf('%s: Reading acquisitions (data)... ', datetime);
    nav_data{idx} = data_dset.readAcquisition(nav_list(idx)); % read all the acquisitions
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
end

%% Get k-space trajectories for each data type
nav_traj = raw_traj.select(find(acq_is_navigation_data));
img_traj = raw_traj.select(find(~acq_is_navigation_data));

%% Parse an ISMRMRD header
%--------------------------------------------------------------------------
% ISMRMRD header
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
% uint16_t discard_pre;                                /**< Samples to be discarded at the beginning of acquisition */
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
nr_phase_encoding_steps = double(max(img_traj.head.idx.kspace_encode_step_1)) + 1;
nr_slice_encoding_steps = double(max(img_traj.head.idx.kspace_encode_step_2)) + 1;
nr_averages             = double(max(img_traj.head.idx.average)) + 1;
nr_slices               = double(max(img_traj.head.idx.slice)) + 1;
nr_contrasts            = double(max(img_traj.head.idx.contrast)) + 1;
nr_phases               = double(max(img_traj.head.idx.phase)) + 1;
nr_repetitions          = double(max(img_traj.head.idx.repetition)) + 1;
nr_sets                 = double(max(img_traj.head.idx.set)) + 1;
nr_segments             = double(max(img_traj.head.idx.segment)) + 1;
center_sample           = double(max(img_traj.head.center_sample));

%--------------------------------------------------------------------------
% Get the dimensionality of the trajectory vector (0 means no trajectory)
%--------------------------------------------------------------------------
trajectory_dimensions = double(max(img_traj.head.trajectory_dimensions));

%--------------------------------------------------------------------------
% discard_pre
%--------------------------------------------------------------------------
discard_pre = double(max(img_traj.head.discard_pre));

%--------------------------------------------------------------------------
% discard_post
%--------------------------------------------------------------------------
discard_post = double(max(img_traj.head.discard_post));

%--------------------------------------------------------------------------
% true_resolution_flag
%--------------------------------------------------------------------------
true_resolution_flag = traj_header.encoding.trajectoryDescription.userParameterLong(1).value;

%--------------------------------------------------------------------------
% fid_navigator_flag
%--------------------------------------------------------------------------
fid_navigator_flag = traj_header.encoding.trajectoryDescription.userParameterLong(2).value;

%--------------------------------------------------------------------------
% noise_scan_flag
%--------------------------------------------------------------------------
noise_scan_flag = traj_header.encoding.trajectoryDescription.userParameterLong(3).value;

%--------------------------------------------------------------------------
% self_navigation_flag
%--------------------------------------------------------------------------
self_navigation_flag = traj_header.encoding.trajectoryDescription.userParameterLong(4).value;

%--------------------------------------------------------------------------
% grad_samples
%--------------------------------------------------------------------------
grad_samples = traj_header.encoding.trajectoryDescription.userParameterLong(5).value;

%--------------------------------------------------------------------------
% adc_samples
%--------------------------------------------------------------------------
adc_samples = traj_header.encoding.trajectoryDescription.userParameterLong(6).value;

%--------------------------------------------------------------------------
% ramp_time
%--------------------------------------------------------------------------
ramp_time = traj_header.encoding.trajectoryDescription.userParameterLong(7).value;
% wrong! Long can't be used in here

%--------------------------------------------------------------------------
% total_time
%--------------------------------------------------------------------------
total_time = traj_header.encoding.trajectoryDescription.userParameterLong(8).value;
% wrong! Long can't be used in here

%--------------------------------------------------------------------------
% grid_scale_factor
%--------------------------------------------------------------------------
grid_scale_factor = traj_header.encoding.trajectoryDescription.userParameterDouble(1).value;

%--------------------------------------------------------------------------
% Get the gradient raster time in [sec]
%--------------------------------------------------------------------------
grad_raster_time = traj_header.encoding.trajectoryDescription.userParameterDouble(2).value; % [sec]

%--------------------------------------------------------------------------
% Get the real dwell time in [sec]
%--------------------------------------------------------------------------
real_dwell_time = traj_header.encoding.trajectoryDescription.userParameterDouble(3).value; % [sec]

%--------------------------------------------------------------------------
% TE_nav in [sec]
%--------------------------------------------------------------------------
TE_nav = traj_header.encoding.trajectoryDescription.userParameterDouble(4).value; % [sec]

%--------------------------------------------------------------------------
% dwell_time_nav in [sec]
%--------------------------------------------------------------------------
dwell_time_nav = traj_header.encoding.trajectoryDescription.userParameterDouble(5).value; % [sec]

%--------------------------------------------------------------------------
% Calculate the gradient duration [sec]
%--------------------------------------------------------------------------
grad_duration = grad_samples * grad_raster_time;

%--------------------------------------------------------------------------
% Calculate the readout time (tau) [sec]
%--------------------------------------------------------------------------
tau = adc_samples * real_dwell_time;

%% Display an ISMRMRD header
fprintf('========================= ISMRMRD header ========================\n');
fprintf('encoded_fov [mm]   = %8.4f %8.4f %8.4f\n', encoded_fov(1) * 1e3, encoded_fov(2) * 1e3, encoded_fov(3) * 1e3);
fprintf('Nkx Nky Nkz        = %d      %d        %d\n', Nkx, Nky, Nkz);
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
fprintf('grad raster time        = %5.2f [usec]\n', grad_raster_time * 1e6);
fprintf('gradient duration       = %5.2f [msec]\n', grad_duration * 1e3);
fprintf('real dwell time         = %5.2f [usec]\n', real_dwell_time * 1e6);
fprintf('readout time            = %5.2f [msec]\n', tau * 1e3);
fprintf('=================================================================\n');

%% Calculate the receiver noise matrix (Nc x Nc)
[Psi,inv_L] = calculate_receiver_noise_matrix(ismrmrd_noise_file);

%% Define parameters for non-Cartesian reconstruction
%--------------------------------------------------------------------------
% Set the number of samples in the row (r) direction of the RCS
%--------------------------------------------------------------------------
N1 = recon_matrix_size(1);

%--------------------------------------------------------------------------
% Set the number of samples in the column (c) direction of the RCS
%--------------------------------------------------------------------------
N2 = recon_matrix_size(2);

%--------------------------------------------------------------------------
% Set the number of samples in the slice (S) direction of the RCS
%--------------------------------------------------------------------------
N3 = recon_matrix_size(3);

%--------------------------------------------------------------------------
% Set the total number of voxels in image space
%--------------------------------------------------------------------------
N  = N1 * N2 * N3;

%% Get FID navigator data (1 x nr_nav_samples x nr_navigators x nr_channels)
%--------------------------------------------------------------------------
% ( 1)not used   x ( 2)readout dimension x ( 3)number of TRs  x ( 4)COIL_DIM
% ( 5)MAPS_DIM   x ( 6)TE_DIM            x ( 7)COEFF_DIM      x ( 8)COEFF2_DIM
% ( 9)ITER_DIM   x (10)CSHIFT_DIM        x (11)TIME_DIM       x (12)TIME2_DIM
% (13)LEVEL_DIM  x (14)SLICE_DIM         x (15)AVG_DIM        x (16)BATCH_DIM
%--------------------------------------------------------------------------
if fid_navigator_flag
    %----------------------------------------------------------------------
    % Get the number of ADC samples for an FID navigator
    %----------------------------------------------------------------------
    nav_adc_samples = double(max(nav_traj.head.number_of_samples));

    %----------------------------------------------------------------------
    % Calculate the index range of readout samples
    %----------------------------------------------------------------------
    nav_discard_pre = 10;
    samples_index_range = ((nav_discard_pre + 1):nav_adc_samples).';
    nr_nav_samples = length(samples_index_range);

    %----------------------------------------------------------------------
    % Calculate the number of FID navigators
    %----------------------------------------------------------------------
    nr_navigators = length(nav_data);

    %----------------------------------------------------------------------
    % Preallocate FID navigator data
    %----------------------------------------------------------------------
    fid = complex(zeros([1 nr_nav_samples nr_navigators nr_channels], 'double'));

    %----------------------------------------------------------------------
    % Sort one acquisition at a time
    %----------------------------------------------------------------------
    for nav = 1:nr_navigators
        tstart = tic; fprintf('%s: Reading FID navigator data (%4d/%4d)... ', datetime, nav, nr_navigators);

        %------------------------------------------------------------------
        % Prewhiten k-space data
        %------------------------------------------------------------------
        readout = nav_data{nav}.data{:}; % number_of_samples x nr_channels
        readout = (inv_L * readout.').';

        %------------------------------------------------------------------
        % Collect FID navigator data
        %------------------------------------------------------------------
        fid(1, :, nav, :) = readout(samples_index_range,:); % nav_adc_samples x nr_channels => nr_nav_samples x nr_channels
        fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
    end
end

%% Perform phase unwrapping of FIDs [rad]
fid_phase_unwrapped = unwrap(angle(fid), [], 2);

%% Calculate a times axis for an FID navigator [sec]
t_nav = TE_nav + (samples_index_range - 1 + 0.5) * dwell_time_nav;

%% Perform a least-squares linear fit
%--------------------------------------------------------------------------
% [phi(t_1)]   [2 * pi * t_1 1][    f0    ]
% [phi(t_2)] = [2 * pi * t_2 1][phi_system]
% [phi(t_N)]   [2 * pi * t_N 1]
%  => y = A * x
%  x = inv(A) * y or x = A \ y
%--------------------------------------------------------------------------
A = [2 * pi * t_nav ones(nr_nav_samples, 1)];
f0 = zeros([1 1 nr_navigators nr_channels], 'double'); % [Hz]
phi_system = zeros([1 1 nr_navigators nr_channels], 'double'); % [rad]

for nav = 1:nr_navigators
    for chan = 1:nr_channels
        tstart = tic; fprintf('%s: Performing a least-squares linear fit (NAV=%d/%d)(CHAN=%2d/%2d)... ', datetime, nav, nr_navigators, chan, nr_channels);
        y = reshape(fid_phase_unwrapped(1, :, nav, chan), [nr_nav_samples 1]);
        x_fit = A \ y;
        f0(1, 1, nav, chan) = x_fit(1);
        phi_system(1, 1, nav, chan) = x_fit(2);
        y_fit = A * x_fit;
        fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
    end
end

%% Display results
FontSize = 12;

for nav = 1:nr_navigators
    figure('Color', 'w', 'Position', [4 7 1236 971]);
    for chan = 1:nr_channels
    y = reshape(fid_phase_unwrapped(1, :, nav, chan), [nr_nav_samples 1]);
    y_fit = 2 * pi * f0(1, 1, nav, chan) * t_nav + phi_system(1, 1, nav, chan);
    subplot(4,4,chan);
    hold on;
    plot(t_nav * 1e3, y, '.-', 'MarkerSize', 12);
    plot(t_nav * 1e3, y_fit, '.-', 'MarkerSize', 12);
    grid on; grid minor;
    set(gca, 'Box', 'On', 'TickLabelInterpreter', 'latex', 'FontSize', FontSize);
    title({sprintf('Ch %d, %5.2f Hz/%4.2f rad', chan, f0(1, 1, nav, chan), phi_system(1, 1, nav, chan))}, 'Interpreter', 'latex', 'FontSize', FontSize);
    xlabel('Time [msec]', 'Interpreter', 'latex', 'FontSize', FontSize);
    ylabel('Phase [rad]', 'Interpreter', 'latex', 'FontSize', FontSize);
    end
    sgtitle({sprintf('$f_0$ navigator %d, least-squares linear fit: $y = 2\\pi f_0 t + \\phi_{\\mathrm{system}}$', nav)}, 'Interpreter', 'latex', 'FontSize', FontSize + 4);
    export_fig(fullfile(output_path, sprintf('f0_navigator%d_least_squares_fit', nav)), '-r300', '-tif', '-c[120, 300, 200, 330]'); % [top,right,bottom,left]
    close gcf;
end

%% Calculate the correlation matrix
f0_ = reshape(f0(1, 1, :, :) , [nr_navigators nr_channels]);
cor = zeros(nr_channels, nr_channels, 'double');
for idx2 = 1:nr_channels
    for idx1 = 1:nr_channels
        cor(idx1,idx2) = f0_(:,idx1).' * f0_(:,idx2) / norm(f0_(:,idx1)) / norm(f0_(:,idx2));
    end
end

FontSize = 12;
figure('Color', 'w');
imagesc(cor);
axis image;
title('Correlation matrix', 'FontSize', FontSize, 'FontWeight', 'normal');
colormap(jet(256));
colorbar;
export_fig(fullfile(output_path, sprintf('correlation_matrix')), '-r300', '-tif');
close gcf;

%% Select receive channels that are highly correlated
coil_selection = true(nr_channels, 1);
threshold = median(cor(:)); % 0.9975
count = 1;
for chan = 1:nr_channels
    if length(find(cor(chan,:) < threshold)) > nr_channels / 2
        coil_selection(chan) = 0;
        count = count + 1;
    end
end

%% Write a .cfl file
%--------------------------------------------------------------------------
% fid (1 x nr_nav_samples x nr_navigators x nr_channels)
%--------------------------------------------------------------------------
cfl_file = fullfile(output_path, 'fid');
tstart = tic; fprintf('%s: Writing a .cfl file: %s... ', datetime, cfl_file);
writecfl(cfl_file, fid(1, :, :, coil_selection));
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%% Perform coil compression
%--------------------------------------------------------------------------
% Usage: cc [-p d] [-M] [-r d:d:d] [-A] [-S] [-G] [-E] <kspace> <coeff|proj_kspace>
% 
% Performs coil compression.
% 
% -p N    perform compression to N virtual channels
% -M      output compression matrix
% -r S    size of calibration region
% -A      use all data to compute coefficients
% -S      type: SVD
% -G      type: Geometric
% -E      type: ESPIRiT
% -h      help
%--------------------------------------------------------------------------
fid_file    = strcat(bart_output_path, 'fid');
fid_cc_file = strcat(bart_output_path, 'fid_cc');
command = sprintf('%s cc -p 1 -S %s %s', bart_command, fid_file, fid_cc_file);
tstart = tic; fprintf('%s:[BART] Performing coil compression:\n%s\n', datetime, command);
[status_cc,result_cc] = system(command);
fprintf('%s: done! (%6.4f/%6.4f sec)\n', datetime, toc(tstart), toc(start_time));

%% Read a .cfl file
cfl_file = fullfile(output_path, 'fid_cc');
tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
fid_cc = readcfl(cfl_file);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%% Perform phase unwrapping of FIDs (1 x nr_nav_samples x nr_navigators)
fid_cc_phase_unwrapped = unwrap(angle(fid_cc), [], 2);

%% Perform a least-squares linear fit
%--------------------------------------------------------------------------
% [phi(t_1)]   [2 * pi * t_1 1][    f0    ]
% [phi(t_2)] = [2 * pi * t_2 1][phi_system]
% [phi(t_N)]   [2 * pi * t_N 1]
%  => y = A * x
%  x = inv(A) * y or x = A \ y
%--------------------------------------------------------------------------
A = [2 * pi * t_nav ones(nr_nav_samples, 1)];
f0_cc = zeros([1 1 nr_navigators], 'double'); % [Hz]
phi_system_cc = zeros([1 1 nr_navigators], 'double'); % [rad]

for nav = 1:nr_navigators
    tstart = tic; fprintf('%s: Performing a least-squares linear fit (NAV=%d/%d)... ', datetime, nav, nr_navigators);
    y = reshape(fid_cc_phase_unwrapped(1, :, nav), [nr_nav_samples 1]);
    x_fit = A \ y;
    f0_cc(1, 1, nav) = x_fit(1);
    phi_system_cc(1, 1, nav) = x_fit(2);
    y_fit = A * x_fit;
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
end

%% Calculate a time axis [sec]
t_nav_data = zeros(nr_navigators, 1, 'double');
for idx = 1:nr_navigators
    t_nav_data(idx) = nav_data{idx}.head.acquisition_time_stamp * 2.5e-3; % [sec]
end
t_nav_data = t_nav_data - t_nav_data(1);

%% Perform a least-squares linear fit
%--------------------------------------------------------------------------
% [f0_cc(t_1)]   [t_1 1][f0_slope ]
% [f0_cc(t_2)] = [t_2 1][f0_offset]
% [f0_cc(t_N)]   [t_N 1]
%  => y = A * x
%  x = inv(A) * y or x = A \ y
%--------------------------------------------------------------------------
A = [t_nav_data ones(nr_navigators, 1)];
f0_nav_data = reshape(f0_cc, [nr_navigators 1]); % [Hz]
x_fit = A \ f0_nav_data;
f0_slope = x_fit(1); % [Hz]
f0_offset = x_fit(2); % [Hz]

%% Write a .cfl file
%--------------------------------------------------------------------------
% t_nav_data (nr_navigators x 1)
%--------------------------------------------------------------------------
cfl_file = fullfile(output_path, 't_nav_data');
tstart = tic; fprintf('%s: Writing a .cfl file: %s... ', datetime, cfl_file);
writecfl(cfl_file, t_nav_data);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%--------------------------------------------------------------------------
% f0_nav_data (nr_navigators x 1)
%--------------------------------------------------------------------------
cfl_file = fullfile(output_path, 'f0_nav_data');
tstart = tic; fprintf('%s: Writing a .cfl file: %s... ', datetime, cfl_file);
writecfl(cfl_file, f0_nav_data);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%% Display results
MarkerSize = 12;
figure('Color', 'w');
ColorOrder = get(gca, 'ColorOrder');
hold on;
%hPlot2 = plot(t_img_data, f0_img_data, '.', 'Color', ColorOrder(2,:), 'MarkerSize', MarkerSize);
hPlot1 = plot(t_nav_data, f0_nav_data, '.', 'Color', ColorOrder(2,:), 'MarkerSize', MarkerSize + 10);
hPlot3 = plot(t_nav_data, A * x_fit, '.-', 'Color', ColorOrder(1,:), 'MarkerSize', MarkerSize);
grid on;
xlabel('Time [sec]', 'Interpreter', 'latex', 'FontSize', FontSize);
ylabel('$f_0$ drift [Hz]', 'Interpreter', 'latex', 'FontSize', FontSize);
set(gca, 'Box', 'On');
legend([hPlot1 hPlot3], 'FID navigator', 'FID navigator fit', 'Spiral arm', 'Interpreter', 'latex', 'FontSize', FontSize, 'Location', 'Northwest');
title({'Linear interpolation of center frequency drifts', sprintf('measured before and after a scan')}, 'FontSize', FontSize, 'FontWeight', 'normal');

export_fig(fullfile(output_path, sprintf('linear_interpolation_center_frequency_drift')), '-r300', '-tif', '-c[0, 100, 0, 20]'); % [top,right,bottom,left]
close gcf;
