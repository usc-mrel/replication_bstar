% demo_step1_bh_prepare_ksp.m
% Written by Nam Gyun Lee
% Email: namgyunl@usc.edu, ggang56@gmail.com (preferred)
% Started: 04/25/2024, Last modified: 04/28/2024

%% Clean slate
close all; clearvars -except json_number nr_json_files json_files json_file; clc;

%% Start a stopwatch timer
start_time = tic;

%% Read a .json file
tstart = tic; fprintf('%s: Reading a .json file: %s... ', datetime, json_file);
fid = fopen(json_file);
json_txt = fread(fid, [1 inf], 'char=>char');
fclose(fid);
json = jsondecode(json_txt);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%% Get parameters from a .json file
%--------------------------------------------------------------------------
% Define the full path of a filename
%--------------------------------------------------------------------------
if ispc
    ismrmrd_data_file  = strrep(json.ismrmrd_data_file, '/', '\');
    ismrmrd_noise_file = strrep(json.ismrmrd_noise_file, '/', '\');
    ismrmrd_traj_file  = strrep(json.ismrmrd_traj_file, '/', '\');
    output_path        = strrep(json.output_path, '/', '\');
else
    ismrmrd_data_file  = json.ismrmrd_data_file;
    ismrmrd_noise_file = json.ismrmrd_noise_file;
    ismrmrd_traj_file  = json.ismrmrd_traj_file;
    output_path        = json.output_path;
end

%--------------------------------------------------------------------------
% coil selection
%--------------------------------------------------------------------------
if isfield(json, 'coil_selection')
    coil_selection = json.coil_selection;
end

%--------------------------------------------------------------------------
% Define the number of echoes
%--------------------------------------------------------------------------
nr_echoes = 2;

%% Make an output path
mkdir(output_path);

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

%% Read acquisitions
tstart = tic; fprintf('%s: Reading acquisitions... ', datetime);
raw_traj = traj_dset.readAcquisition(); % read all the acquisitions
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

tstart = tic; fprintf('%s: Reading acquisitions... ', datetime);
raw_data = data_dset.readAcquisition(); % read all the acquisitions
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%--------------------------------------------------------------------------
% Get data type
%--------------------------------------------------------------------------
acq_is_navigation_data = raw_traj.head.flagIsSet('ACQ_IS_NAVIGATION_DATA'); % FID navigator
acq_user1              = raw_traj.head.flagIsSet('ACQ_USER1'); % self-navigation

%% Get navigator data and imaging data
nav_data = raw_data.select(find(acq_is_navigation_data));
img_data = raw_data.select(find(~(acq_is_navigation_data | acq_user1)));

%% Get k-space trajectories for each data type
nav_traj  = raw_traj.select(find(acq_is_navigation_data));
self_traj = raw_traj.select(find(acq_user1));
img_traj  = raw_traj.select(find(~(acq_is_navigation_data | acq_user1)));

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
adc_samples  = double(max(img_traj.head.number_of_samples));
discard_pre  = double(max(img_traj.head.discard_pre));
discard_post = double(max(img_traj.head.discard_post));

%--------------------------------------------------------------------------
% Number of half-spokes per interleaf
%--------------------------------------------------------------------------
nr_readouts = traj_header.encoding.encodingLimits.kspace_encoding_step_1.maximum + 1;

%--------------------------------------------------------------------------
% Number of interleaves
%--------------------------------------------------------------------------
nr_interleaves = traj_header.encoding.encodingLimits.segment.maximum + 1;

%--------------------------------------------------------------------------
% Number of receive channels
%--------------------------------------------------------------------------
nr_channels = data_header.acquisitionSystemInformation.receiverChannels;

%--------------------------------------------------------------------------
% Define the total number of half-spokes
%--------------------------------------------------------------------------
nr_radial_views = nr_readouts * nr_interleaves;

%% Calculate the receiver noise matrix
[Psi,inv_L] = calculate_receiver_noise_matrix(ismrmrd_noise_file);

%% Get k-space data (1 x adc_samples x nr_radial_views x nr_channels)
%--------------------------------------------------------------------------
% ( 1)not used   x ( 2)readout dimension x ( 3)number of TRs  x ( 4)COIL_DIM
% ( 5)MAPS_DIM   x ( 6)TE_DIM            x ( 7)COEFF_DIM      x ( 8)COEFF2_DIM
% ( 9)ITER_DIM   x (10)CSHIFT_DIM        x (11)TIME_DIM       x (12)TIME2_DIM
% (13)LEVEL_DIM  x (14)SLICE_DIM         x (15)AVG_DIM        x (16)BATCH_DIM
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Calculate the index range of readout samples for reconstruction
%--------------------------------------------------------------------------
samples_index_range = ((discard_pre + 1):(adc_samples + discard_pre)).';

%--------------------------------------------------------------------------
% Preallocate k-space data (adc_samples x nr_channels x nr_radial_views)
%--------------------------------------------------------------------------
ksp = complex(zeros([adc_samples nr_channels nr_radial_views], 'single'));

%--------------------------------------------------------------------------
% Sort one acquisition at a time
%--------------------------------------------------------------------------
for idx = 1:nr_radial_views
    tstart = tic; fprintf('%s: Reading k-space data (%5d/%5d)... ', datetime, idx, nr_radial_views);

    %----------------------------------------------------------------------
    % Get a phase encoding line number [1, ..., nr_readouts]
    %----------------------------------------------------------------------
    phase_encoding_line_number = img_traj.head.idx.kspace_encode_step_1(idx) + 1;

    %----------------------------------------------------------------------
    % Get a segment number [1, ..., nr_interleaves]
    %----------------------------------------------------------------------
    segment_number = img_traj.head.idx.segment(idx) + 1;

    %----------------------------------------------------------------------
    % Calculate a view number [1, ..., nr_radial_views]
    %----------------------------------------------------------------------
    view_number = phase_encoding_line_number + (segment_number - 1) * nr_readouts;

    %----------------------------------------------------------------------
    % Prewhiten k-space data
    %----------------------------------------------------------------------
    readout = img_data.data{idx}(samples_index_range,:); % number_of_samples x nr_channels => adc_samples x nr_channels
    readout = (inv_L * readout.').';

    %----------------------------------------------------------------------
    % Collect k-space data
    %----------------------------------------------------------------------
    ksp(:,:,view_number) = readout; % adc_samples x nr_channels
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
end

%% Permute "ksp"
ksp = reshape(permute(ksp, [1 3 2]), [1 adc_samples nr_radial_views nr_channels]);

%% Perform coil selection
if exist('coil_selection', 'var')
    ksp = ksp(:,:,:,coil_selection);
    nr_channels = length(coil_selection);
end

%% Write a .cfl file
adc_samples_per_echo = adc_samples / nr_echoes;

for eco = 1:nr_echoes
    %----------------------------------------------------------------------
    % Calculate the index range of samples per echo
    %----------------------------------------------------------------------
    echo_range = (1:adc_samples_per_echo) + (eco - 1) * adc_samples_per_echo;

    %----------------------------------------------------------------------
    % ksp (1 x adc_samples_per_echo x nr_radial_views x nr_channels)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path, sprintf('ksp_eco%d', eco));
    tstart = tic; fprintf('%s:(ECO=%d/%d) Writing a .cfl file: %s... ', datetime, eco, nr_echoes, cfl_file);
    writecfl(cfl_file, ksp(:,echo_range,:,:));
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
end
