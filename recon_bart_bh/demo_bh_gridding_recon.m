% demo_bh_gridding_recon.m
% Written by Nam Gyun Lee
% Email: namgyunl@usc.edu, ggang56@gmail.com (preferred)
% Started: 07/27/2022, Last modified: 08/02/2023

%% Clean slate
close all; clearvars -except json_nr json_files json_file; clc;

%% Read a .json file
start_time = tic;
tstart = tic; fprintf('%s: Reading a .json file: %s... \n', datetime, json_file);
fid = fopen(json_file); 
json_txt = fread(fid, [1 inf], 'char=>char'); 
fclose(fid);
json = jsondecode(json_txt);
fprintf('%s: done! (%6.4f/%6.4f sec)\n', datetime, toc(tstart), toc(start_time));

%--------------------------------------------------------------------------
% Define the full path of a filename
%--------------------------------------------------------------------------
if ispc
    siemens_twix_file  = strrep(json.siemens_twix_file, '/', '\');
    ismrmrd_data_file  = strrep(json.ismrmrd_data_file, '/', '\');
    ismrmrd_noise_file = strrep(json.ismrmrd_noise_file, '/', '\');
    trj_file           = strrep(json.trj_file, '/', '\');
    output_path        = strrep(json.output_path, '/', '\');
else
    siemens_twix_file  = json.siemens_twix_file;
    ismrmrd_data_file  = json.ismrmrd_data_file;
    ismrmrd_noise_file = json.ismrmrd_noise_file;
    trj_file           = json.trj_file;
    output_path        = json.output_path;
end

%--------------------------------------------------------------------------
% Define the BART directory
%--------------------------------------------------------------------------
bart_path = json.bart_path;

%--------------------------------------------------------------------------
% Reconstruction parameters
%--------------------------------------------------------------------------
recon_matrix_size   = json.recon_parameters.recon_matrix_size.'; % reconstruction matrix size
dicom_matrix_size   = json.recon_parameters.dicom_matrix_size.'; % DICOM matrix size
recon_interp_factor = json.recon_parameters.recon_interp_factor; % reconstruction interpolation factor: recon_resolution = encoded_resolution / recon_interp_factor;
max_iter            = json.recon_parameters.max_iter;            % max. number of iterations for pics in BART

%--------------------------------------------------------------------------
% Sequence parameters
%--------------------------------------------------------------------------
nr_contrasts = json.seq.nr_contrasts;
resolution   = json.seq.resolution;

%% Set up a BART command
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

%% Read a Siemens .dat file
fprintf('%s: Reading a Siemens .dat file: %s\n', datetime, siemens_twix_file);
twix = mapVBVD(siemens_twix_file);
if length(twix) > 1
    twix = twix{end};
end

%% Calculate reconstruction parameters
traj_scale_factor = ceil(recon_matrix_size / recon_interp_factor);
N1 = recon_matrix_size(1);
N2 = recon_matrix_size(2);
N3 = recon_matrix_size(3);

%% Calculate imaging parameters from a .seq file
encoded_resolution = ones(1,3) * resolution;

recon_resolution = encoded_resolution ./ recon_interp_factor;

%% NUFFT reconstruction
for echo = 1:nr_contrasts
    %% Define filenames
    traj_unscaled_filename = sprintf('traj_echo%d_unscaled', echo);
    traj_filename          = sprintf('traj_echo%d_%dx%dx%d', echo, traj_scale_factor(1), traj_scale_factor(2), traj_scale_factor(3));
    img_rss_filename       = sprintf('img_rss_echo%d_%dx%dx%d_res%5.3f', echo, N1, N2, N3, recon_resolution(1) * 1e3);

    traj_file    = strcat(bart_output_path, traj_filename);
    dcf_file     = strcat(bart_output_path, sprintf('dcf_echo%d', echo));
    ksp_file     = strcat(bart_output_path, sprintf('ksp_echo%d', echo));
    ksp_dcf_file = strcat(bart_output_path, sprintf('ksp_dcf_echo%d', echo));
    cimg_file    = strcat(bart_output_path, sprintf('cimg_echo%d_%dx%dx%d_res%5.3f', echo, N1, N2, N3, recon_resolution(1) * 1e3));
    img_rss_file = strcat(bart_output_path, img_rss_filename);

    %% Scale k-space trajectories
    %----------------------------------------------------------------------
    % Read "unscaled" k-space trajectories [-0.5,0.5] (3 x Nk x NiNs)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path, traj_unscaled_filename);
    tstart = tic; fprintf('%s:(e=%d/%d) Reading a .cfl file: %s... \n', datetime, echo, nr_contrasts, cfl_file);
    traj_unscaled = readcfl(cfl_file);
    fprintf('%s: done! (%6.4f/%6.4f sec)\n', datetime, toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % Scale k-space trajectories for BART (3 x Nk x NiNs)
    % [-0.5,0.5] => [-0.5,0.5] * N
    %----------------------------------------------------------------------
    traj = zeros(size(traj_unscaled), 'single');
    traj(1,:,:) = traj_unscaled(1,:,:) * traj_scale_factor(1);
    traj(2,:,:) = traj_unscaled(2,:,:) * traj_scale_factor(2); % AP direction
    traj(3,:,:) = traj_unscaled(3,:,:) * traj_scale_factor(3);

    %----------------------------------------------------------------------
    % Save scaled k-space trajectories as a .cfl file
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path, traj_filename);
    tstart = tic; fprintf('%s:(e=%d/%d) Writing a .cfl file: %s... \n', datetime, echo, nr_contrasts, cfl_file);
    writecfl(cfl_file, traj);
    fprintf('%s: done! (%6.4f/%6.4f sec)\n', datetime, toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % Multiply k-space with dcf
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
    %------------------------------------------------------------------
    command = sprintf('%s fmac %s %s %s', bart_command, ksp_file, dcf_file, ksp_dcf_file);
    tstart = tic; fprintf('%s:(e=%d/%d)[BART] Multiplying ksp with dcf... \n%s\n', datetime, echo, nr_contrasts, command);
    [status_fmac,result_fmac] = system(command);
    fprintf('%s: done! (%6.4f/%6.4f sec)\n', datetime, toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % Perform NUFFT reconstruction
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
    %------------------------------------------------------------------
    command = sprintf('%s nufft -a -x %d:%d:%d -t %s %s %s', bart_command, N1, N2, N3, traj_file, ksp_dcf_file, cimg_file);
    tstart = tic; fprintf('%s:(e=%d/%d)[BART] Performing (adjoint) NUFFT reconstruction (high-resolution)... \n%s\n', datetime, echo, nr_contrasts, command);
    [status_nufft,result_nufft] = system(command);
    fprintf('%s: done! (%6.4f/%6.4f sec)\n', datetime, toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % Perform the root-sum-of-squares of coil images
    %----------------------------------------------------------------------
    command = sprintf('%s rss 8 %s %s', bart_command, cimg_file, img_rss_file);
    tstart = tic; fprintf('%s:(e=%d/%d)[BART] Performing the root-sum-of-squares of coil images... \n%s\n', datetime, echo, nr_contrasts, command);
    [status_rss,result_rss] = system(command);
    fprintf('%s: done! (%6.4f/%6.4f sec)\n', datetime, toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % Save as .dicom files
    %----------------------------------------------------------------------
    recon_type = sprintf('rss_echo%d', echo);
    cfl_file = fullfile(output_path, img_rss_filename);
    dicom_directory = fullfile(output_path, sprintf('%s_dicom%d', img_rss_filename, dicom_matrix_size(1)));
    save_images_as_dicom_files(twix, dicom_matrix_size, recon_resolution, recon_type, dicom_directory, cfl_file);
end

%% Combine both echoes
img_rss_combined = complex(zeros(N1, N2, N3, 'single'));
for echo = 1:nr_contrasts
    img_rss_filename = sprintf('img_rss_echo%d_%dx%dx%d_res%5.3f', echo, N1, N2, N3, recon_resolution(1) * 1e3);
    img_rss_combined = img_rss_combined + readcfl(fullfile(output_path, img_rss_filename));
end

%--------------------------------------------------------------------------
% Define the full path of a .cfl file without file extension
%--------------------------------------------------------------------------
img_rss_combined_filename = sprintf('img_rss_combined_%dx%dx%d_res%5.3f', N1, N2, N3, recon_resolution(1) * 1e3);
cfl_file = fullfile(output_path, img_rss_combined_filename);

%--------------------------------------------------------------------------
% Save as a .cfl file
%--------------------------------------------------------------------------
tstart = tic; fprintf('%s: Writing a .cfl file: %s... \n', datetime, cfl_file);
writecfl(cfl_file, img_rss_combined);
fprintf('%s: done! (%6.4f/%6.4f sec)\n', datetime, toc(tstart), toc(start_time));

%--------------------------------------------------------------------------
% Save as .dicom files
%--------------------------------------------------------------------------
recon_type = 'rss_combined';
dicom_directory = fullfile(output_path, sprintf('%s_dicom%d', img_rss_combined_filename, dicom_matrix_size(1)));
save_images_as_dicom_files(twix, dicom_matrix_size, recon_resolution, recon_type, dicom_directory, cfl_file);

