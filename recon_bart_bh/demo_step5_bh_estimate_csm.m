% demo_step5_bh_estimate_csm.m
% Written by Nam Gyun Lee
% Email: namgyunl@usc.edu, ggang56@gmail.com (preferred)
% Started: 04/27/2024, Last modified: 04/28/2024

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
    siemens_twix_file  = strrep(json.siemens_twix_file, '/', '\');
    ismrmrd_data_file  = strrep(json.ismrmrd_data_file, '/', '\');
    ismrmrd_noise_file = strrep(json.ismrmrd_noise_file, '/', '\');
    output_path        = strrep(json.output_path, '/', '\');
else
    siemens_twix_file  = json.siemens_twix_file;
    ismrmrd_data_file  = json.ismrmrd_data_file;
    ismrmrd_noise_file = json.ismrmrd_noise_file;
    output_path        = json.output_path;
end

%--------------------------------------------------------------------------
% Define the BART path
%--------------------------------------------------------------------------
bart_path = json.bart_path;

%--------------------------------------------------------------------------
% Reconstruction parameters
%--------------------------------------------------------------------------
recon_matrix_size   = json.recon_conf.recon_matrix_size;   % reconstruction matrix size
dicom_matrix_size   = json.recon_conf.dicom_matrix_size;   % DICOM matrix size
recon_interp_factor = json.recon_conf.recon_interp_factor; % reconstruction interpolation factor: recon_resolution = encoded_resolution / recon_interp_factor;
cal_size            = json.recon_conf.cal_size;            % size of the calibration region
max_iter            = json.recon_conf.max_iter;            % max. number of iterations for pics in BART

%--------------------------------------------------------------------------
% NLINV parameters
%--------------------------------------------------------------------------
noir_conf = json.noir_conf;

%--------------------------------------------------------------------------
% Sequence parameters
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

%% Define reconstruction parameters
traj_scale_factor = round(recon_matrix_size ./ recon_interp_factor);
N1 = recon_matrix_size(1);
N2 = recon_matrix_size(2);
N3 = recon_matrix_size(3);

%% Estimate coil sensitivity maps per echo
for eco = 1:nr_echoes

    %% Define filenames
    traj_unscaled_filename = sprintf('traj_unscaled_eco%d', eco);
    traj_filename          = sprintf('traj_eco%d_%dx%dx%d', eco, traj_scale_factor(1), traj_scale_factor(2), traj_scale_factor(3));
    cimg_filename          = sprintf('cimg_eco%d_%dx%dx%d', eco, N1, N2, N3);
    sens_fiename           = sprintf('sens_eco%d_%dx%dx%d', eco, N1, N2, N3);
    img_calib_filename     = sprintf('img_eco%d_calib', eco);
    sens_calib_filename    = sprintf('sens_eco%d_calib', eco);

    traj_file       = strcat(bart_output_path, traj_filename);
    ksp_file        = strcat(bart_output_path, sprintf('ksp_eco%d', eco));
    img_calib_file  = strcat(bart_output_path, img_calib_filename);
    sens_calib_file = strcat(bart_output_path, sens_calib_filename);

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

    %% Write a .cfl file
    cfl_file = fullfile(output_path, traj_filename);
    tstart = tic; fprintf('%s:(ECO=%d/%d) Writing a .cfl file: %s... ', datetime, eco, nr_echoes, cfl_file);
    writecfl(cfl_file, traj);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %% Estimate coil sensitivity maps using NLINV
    %----------------------------------------------------------------------
    % Usage: nlinv [-i d] [-d d] [-c] [-N] [-m d] [-U] [-f f] [-p <file>]
    %              [-t <file>] [-I <file>] [-g] [-S] [--lowmem] [-x d:d:d]
    %              <kspace> <output> [<sensitivities>]
    %
    % Jointly estimate image and sensitivities with nonlinear
    % inversion using {iter} iteration steps. Optionally outputs
    % the sensitivities.
    %
    % -i iter     Number of Newton steps
    % -d level    Debug level
    % -c          Real-value constraint
    % -N          Do not normalize image with coil sensitivities
    % -m nmaps    Number of ENLIVE maps to use in reconstruction
    % -U          Do not combine ENLIVE maps in output
    % -f FOV      restrict FOV
    % -p file     pattern / transfer function
    % -t file     kspace trajectory
    % -I file     File for initialization
    % -g          use gpu
    % -S          Re-scale image after reconstruction
    % --lowmem    Use low-mem mode of the nuFFT
    % -x x:y:z    Explicitly specify image dimensions
    % -h          help
    %----------------------------------------------------------------------
    command = sprintf('%s nlinv -a %d -b %d -S -d5 -i%d -x %d:%d:%d -t %s %s %s %s', bart_command, noir_conf.a, noir_conf.b, noir_conf.max_iter, cal_size(1), cal_size(2), cal_size(3), traj_file, ksp_file, img_calib_file, sens_calib_file);
    tstart = tic; fprintf('%s:(ECO=%d/%d)[BART] Estimating coil sensitivities using NLINV:\n%s\n', datetime, eco, nr_echoes, command);
    [status_nlinv,result_nlinv] = system(command);
    fprintf('%s: done! (%6.4f/%6.4f sec)\n', datetime, toc(tstart), toc(start_time));

    %% Read a .cfl file
    cfl_file = fullfile(output_path, sens_calib_filename);
    tstart = tic; fprintf('%s:(ECO=%d/%d) Reading a .cfl file: %s... ', datetime, eco, nr_echoes, cfl_file);
    sens_calib = readcfl(cfl_file);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %% Tranform from image space to k-space (image space <=> k-space)
    %----------------------------------------------------------------------
    % Usage: fft [-u] [-i] [-n] bitmask <input> <output>
    %
    % Performs a fast Fourier transform (FFT) along selected dimensions.
    %
    % -u    unitary
    % -i    inverse
    % -n    un-centered
    % -h    help
    %----------------------------------------------------------------------
    tstart = tic; fprintf('%s:(ECO=%d/%d) Transforming from image space to k-space... ', datetime, eco, nr_echoes);
    ksp_sens_calib = sens_calib;
    clear sens_calib;
    for dim = 1:3
        ksp_sens_calib = 1 / sqrt(cal_size(dim) * 2) * fftshift(fft(ifftshift(ksp_sens_calib, dim), [], dim), dim);
    end
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %% Zeropad to a 2X grid in k-space
    tstart = tic; fprintf('%s:(ECO=%d/%d) Zeropadding to a 2X grid in k-space... ', datetime, eco, nr_echoes);
    idx1_range = (-floor(cal_size(1)*2/2):ceil(cal_size(1)*2/2)-1).' + floor(recon_matrix_size(1)*2/2) + 1;
    idx2_range = (-floor(cal_size(2)*2/2):ceil(cal_size(2)*2/2)-1).' + floor(recon_matrix_size(2)*2/2) + 1;
    idx3_range = (-floor(cal_size(3)*2/2):ceil(cal_size(3)*2/2)-1).' + floor(recon_matrix_size(3)*2/2) + 1;
    ksp_sens_calib_zpad = complex(zeros([recon_matrix_size * 2 size(ksp_sens_calib,4)], 'single'));
    ksp_sens_calib_zpad(idx1_range, idx2_range, idx3_range, :) = ksp_sens_calib;
    clear ksp_sens_calib;
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %% Tranform from k-space to image space (image space <=> k-space)
    tstart = tic; fprintf('%s:(ECO=%d/%d) Transforming from k-space to image space... ', datetime, eco, nr_echoes);
    sens_2Xgrid = ksp_sens_calib_zpad;
    clear ksp_sens_calib_zpad;
    for dim = 1:3
        sens_2Xgrid = sqrt(recon_matrix_size(dim) * 2) * fftshift(ifft(ifftshift(sens_2Xgrid, dim), [], dim), dim);
    end
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %% Crop a 2X grid to a 1X grid in image space
    tstart = tic; fprintf('%s:(ECO=%d/%d) Cropping a 2X grid to a 1X grid in image space... ', datetime, eco, nr_echoes);
    idx1_range = (-floor(recon_matrix_size(1)/2):ceil(recon_matrix_size(1)/2)-1).' + floor(recon_matrix_size(1)*2/2) + 1;
    idx2_range = (-floor(recon_matrix_size(2)/2):ceil(recon_matrix_size(2)/2)-1).' + floor(recon_matrix_size(2)*2/2) + 1;
    idx3_range = (-floor(recon_matrix_size(3)/2):ceil(recon_matrix_size(3)/2)-1).' + floor(recon_matrix_size(3)*2/2) + 1;
    sens_1Xgrid = sens_2Xgrid(idx1_range, idx2_range, idx3_range, :);
    clear sens_2Xgrid;
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %% Normalize coil sensitivity maps estimated by NLINV
    tstart = tic; fprintf('%s:(ECO=%d/%d) Normalizing coil sensitivity maps estimated by NLINV... ', datetime, eco, nr_echoes);
    sens = bsxfun(@rdivide, sens_1Xgrid, sqrt(sum(abs(sens_1Xgrid).^2, 4)));
    clear sens_1Xgrid;
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %% Write a .cfl file
    cfl_file = fullfile(output_path, sens_fiename);
    tstart = tic; fprintf('%s:(ECO=%d/%d) Writing a .cfl file: %s... ', datetime, eco, nr_echoes, cfl_file);
    writecfl(cfl_file, sens);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %% Delete a .cfl file
    delete(fullfile(output_path, sprintf('%s*', img_calib_filename)));
    delete(fullfile(output_path, sprintf('%s*', sens_calib_filename)));
end
