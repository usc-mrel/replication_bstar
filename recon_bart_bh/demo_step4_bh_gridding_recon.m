% demo_step4_bh_gridding_recon.m
% Written by Nam Gyun Lee
% Email: namgyunl@usc.edu, ggang56@gmail.com (preferred)
% Started: 04/26/2024, Last modified: 04/28/2024

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
    siemens_twix_file = strrep(json.siemens_twix_file, '/', '\');
    ismrmrd_traj_file = strrep(json.ismrmrd_traj_file, '/', '\');
    output_path       = strrep(json.output_path, '/', '\');
else
    siemens_twix_file = json.siemens_twix_file;
    ismrmrd_traj_file = json.ismrmrd_traj_file;
    output_path       = json.output_path;
end

%--------------------------------------------------------------------------
% Define the BART path
%--------------------------------------------------------------------------
bart_path = json.bart_path;

%--------------------------------------------------------------------------
% Reconstruction parameters
%--------------------------------------------------------------------------
recon_matrix_size   = json.recon_conf.recon_matrix_size.';   % reconstruction matrix size
dicom_matrix_size   = json.recon_conf.dicom_matrix_size.';   % DICOM matrix size
recon_interp_factor = json.recon_conf.recon_interp_factor.'; % reconstruction interpolation factor: recon_resolution = encoded_resolution / recon_interp_factor;

%--------------------------------------------------------------------------
% Define the number of echoes
%--------------------------------------------------------------------------
nr_echoes = 2;

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

recon_resolution = encoded_resolution;

%% Read a Siemens .dat file
fprintf('%s: Reading a Siemens .dat file: %s\n', datetime, siemens_twix_file);
twix = mapVBVD(siemens_twix_file);
if length(twix) > 1
    twix = twix{end};
end

%% Define reconstruction parameters
traj_scale_factor = round(recon_matrix_size ./ recon_interp_factor);
N1 = recon_matrix_size(1);
N2 = recon_matrix_size(2);
N3 = recon_matrix_size(3);

%% Perform NUFFT reconstruction per echo
for eco = 1:nr_echoes

    %% Define filenames
    traj_unscaled_filename = sprintf('traj_unscaled_eco%d', eco);
    traj_filename          = sprintf('traj_eco%d_%dx%dx%d', eco, traj_scale_factor(1), traj_scale_factor(2), traj_scale_factor(3));
    ksp_dcf_filename       = sprintf('ksp_dcf_eco%d', eco);
    cimg_filename          = sprintf('cimg_eco%d_%dx%dx%d', eco, N1, N2, N3);
    img_rss_filename       = sprintf('img_rss_eco%d_%dx%dx%d', eco, N1, N2, N3);

    traj_file    = strcat(bart_output_path, traj_filename);
    dcf_file     = strcat(bart_output_path, sprintf('dcf_eco%d_effMtx%d', eco, N1));
    ksp_file     = strcat(bart_output_path, sprintf('ksp_eco%d', eco));
    ksp_dcf_file = strcat(bart_output_path, ksp_dcf_filename);
    cimg_file    = strcat(bart_output_path, cimg_filename);
    img_rss_file = strcat(bart_output_path, img_rss_filename);

    %% Read "unscaled" k-space trajectories [-0.5,0.5]
    cfl_file = fullfile(output_path, traj_unscaled_filename);
    tstart = tic; fprintf('%s:(ECO=%d/%d) Reading a .cfl file: %s... ', datetime, eco, nr_echoes, cfl_file);
    traj_unscaled = readcfl(cfl_file);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %% Scale k-space trajectories for BART: [-0.5,0.5] => [-0.5,0.5] * N
    traj = zeros(size(traj_unscaled), 'single');
    traj(1,:,:) = traj_unscaled(1,:,:) * traj_scale_factor(1);
    traj(2,:,:) = traj_unscaled(2,:,:) * traj_scale_factor(2);
    traj(3,:,:) = traj_unscaled(3,:,:) * traj_scale_factor(3);

    %% Write as a .cfl file
    cfl_file = fullfile(output_path, traj_filename);
    tstart = tic; fprintf('%s:(ECO=%d/%d) Writing a .cfl file: %s... ', datetime, eco, nr_echoes, cfl_file);
    writecfl(cfl_file, traj);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %% Multiply k-space with dcf
    %----------------------------------------------------------------------
    % Usage: fmac [-A] [-C] [-s d] <input1> [<input2>] <output>
    %
    % Multiply <input1> and <input2> and accumulate in <output>.
    % If <input2> is not specified, assume all-ones.
    %
    % -A      add to existing output (instead of overwriting)
    % -C      conjugate input2
    % -s b    squash dimensions selected by bitmask b
    % -h      help
    %----------------------------------------------------------------------
    command = sprintf('%s fmac %s %s %s', bart_command, ksp_file, dcf_file, ksp_dcf_file);
    tstart = tic; fprintf('%s:(ECO=%d/%d)[BART] Multiplying ksp with dcf:\n%s\n', datetime, eco, nr_echoes, command);
    [status_fmac,result_fmac] = system(command);
    fprintf('%s: done! (%6.4f/%6.4f sec)\n', datetime, toc(tstart), toc(start_time));

    %% Perform NUFFT reconstruction
    %----------------------------------------------------------------------
    % Usage: nufft [-a] [-i] [-x d:d:d] [-t] [-r] [-c] [-l f] [-m d] [-P] [-s]
    %              [-g] [-1] [--lowmem] [--no-precomp] [-B <file>] [-p <file>]
    %              <traj> <input> <output>
    %
    % Perform non-uniform Fast Fourier Transform.
    %
    % -a              adjoint
    % -i              inverse
    % -x x:y:z        dimensions
    % -t              Toeplitz embedding for inverse NUFFT
    % -r              turn-off Toeplitz embedding for inverse NUFFT
    % -c              Preconditioning for inverse NUFFT
    % -l lambda       l2 regularization
    % -m iter         max. number of iterations (inverse only)
    % -P              periodic k-space
    % -s              DFT
    % -g              GPU (only inverse)
    % -1              use/return oversampled grid
    % --lowmem        Use low-mem mode of the nuFFT
    % --no-precomp    Use low-low-mem mode of the nuFFT
    % -B file         temporal (or other) basis
    % -p file         weighting of nufft
    % -h              help
    %----------------------------------------------------------------------
    command = sprintf('%s nufft -a -x %d:%d:%d -t %s %s %s', bart_command, N1, N2, N3, traj_file, ksp_dcf_file, cimg_file);
    tstart = tic; fprintf('%s:(ECO=%d/%d)[BART] Performing (adjoint) NUFFT reconstruction:\n%s\n', datetime, eco, nr_echoes, command);
    [status_nufft,result_nufft] = system(command);
    fprintf('%s: done! (%6.4f/%6.4f sec)\n', datetime, toc(tstart), toc(start_time));

    %% Perform the root-sum-of-squares of coil images
    command = sprintf('%s rss 8 %s %s', bart_command, cimg_file, img_rss_file);
    tstart = tic; fprintf('%s:(ECO=%d/%d)[BART] Performing the root-sum-of-squares of coil images:\n%s\n', datetime, eco, nr_echoes, command);
    [status_rss,result_rss] = system(command);
    fprintf('%s: done! (%6.4f/%6.4f sec)\n', datetime, toc(tstart), toc(start_time));

    %% Delete a .cfl file
    delete(fullfile(output_path, sprintf('%s*', ksp_dcf_filename)));
    delete(fullfile(output_path, sprintf('%s*', cimg_filename)));

    %% Save as a .dicom file
    recon_type = sprintf('rss_eco%d', eco);
    cfl_file = fullfile(output_path, img_rss_filename);
    dicom_directory = fullfile(output_path, sprintf('%s_dicom%d', img_rss_filename, dicom_matrix_size(1)));
    save_images_as_dicom_files(twix, dicom_matrix_size, recon_resolution, recon_type, dicom_directory, cfl_file);
end
