% demo_bh_estimate_csm.m
% Written by Nam Gyun Lee
% Email: namgyunl@usc.edu, ggang56@gmail.com (preferred)
% Started: 07/27/2022, Last modified: 08/04/2023

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
cal_size            = json.recon_parameters.cal_size.';          % size of the calibration region
max_iter            = json.recon_parameters.max_iter;            % max. number of iterations for pics in BART
kappa               = json.recon_parameters.kappa;               % kappa: kappa=0 => W=I

%--------------------------------------------------------------------------
% NLINV parameters
%--------------------------------------------------------------------------
noir_conf = json.noir_conf;

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

%% Reconstruct each binned k-space data separately
for echo = 1:nr_contrasts
    %% Define filenames
    traj_unscaled_filename    = sprintf('traj_echo%d_unscaled', echo);
    traj_filename             = sprintf('traj_echo%d_%dx%dx%d', echo, traj_scale_factor(1), traj_scale_factor(2), traj_scale_factor(3));
    cimg_filename             = sprintf('cimg_echo%d_%dx%dx%d_res%5.3f', echo, N1, N2, N3, recon_resolution(1) * 1e3);
    sens_ecalib_fiename       = sprintf('sens_echo%d_ecalib_%dx%dx%d_res%5.3f', echo, N1, N2, N3, recon_resolution(1) * 1e3);
    sens_nlinv_fiename        = sprintf('sens_echo%d_nlinv_%dx%dx%d_res%5.3f', echo, N1, N2, N3, recon_resolution(1) * 1e3);
    img_nufft_ecalib_filename = sprintf('img_nufft_echo%d_%dx%dx%d_res%5.3f_ecalib', echo, N1, N2, N3, recon_resolution(1) * 1e3);
    img_nufft_nlinv_filename  = sprintf('img_nufft_echo%d_%dx%dx%d_res%5.3f_nlinv_a%d_b%d_i%d', echo, N1, N2, N3, recon_resolution(1) * 1e3, noir_conf.a, noir_conf.b, noir_conf.max_iter);

    traj_file        = strcat(bart_output_path, traj_filename);
    ksp_file         = strcat(bart_output_path, sprintf('ksp_echo%d', echo));
    cimg_calib_file  = strcat(bart_output_path, sprintf('cimg_echo%d_calib', echo));
    kgrid_file       = strcat(bart_output_path, sprintf('kgrid_echo%d', echo));
    kgrid_zpad_file  = strcat(bart_output_path, sprintf('kgrid_echo%d_zpad', echo));
    sens_ecalib_file = strcat(bart_output_path, sens_ecalib_fiename);
    ev_maps_file     = strcat(bart_output_path, sprintf('ev_maps_echo%d', echo));

    img_nlinv_calib_file  = strcat(bart_output_path, sprintf('img_echo%d_nlinv_calib', echo));
    sens_nlinv_calib_file = strcat(bart_output_path, sprintf('sens_echo%d_nlinv_calib', echo));
    sens_nlinv_file       = strcat(bart_output_path, sens_nlinv_fiename);

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

    %% Estimate coil sensitivity maps using ESPIRiT calibration
    if 0
        %------------------------------------------------------------------
        % Reconstruct a low-resolution image
        %------------------------------------------------------------------
        % Usage: nufft [-a] [-i] [-x d:d:d] [-t] [-r] [-c] [-l f] [-m d] [-P] [-s] [-g]
        %              [-1] [--lowmem] [--no-precomp] [-B <file>] [-p <file>]
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
        command = sprintf('%s nufft -i -r -x %d:%d:%d %s %s %s', bart_command, cal_size(1), cal_size(2), cal_size(3), traj_file, ksp_file, cimg_calib_file);
        tstart = tic; fprintf('%s:(e=%d/%d)[BART] Performing (iterative) NUFFT reconstruction (low-resolution)... \n%s\n', datetime, echo, nr_contrasts, command);
        [status_nufft,result_nufft] = system(command);
        fprintf('%s: done! (%6.4f/%6.4f sec)\n', datetime, toc(tstart), toc(start_time));

        %------------------------------------------------------------------
        % Tranform from image space to k-space (image space <=> k-space)
        % Prepare gridded k-space using the FT definition in ESPIRiT
        %------------------------------------------------------------------
        command = sprintf('%s fft -u 7 %s %s', bart_command, cimg_calib_file, kgrid_file);
        tstart = tic; fprintf('%s:(e=%d/%d)[BART] Transforming from image space to k-space... \n%s\n', datetime, echo, nr_contrasts, command);
        [status_fft,result_fft] = system(command);
        fprintf('%s: done! (%6.4f/%6.4f sec)\n', datetime, toc(tstart), toc(start_time));

        %------------------------------------------------------------------
        % Zeropad to full size
        %------------------------------------------------------------------
        command = sprintf('%s resize -c 0 %d 1 %d 2 %d %s %s', bart_command, N1, N2, N3, kgrid_file, kgrid_zpad_file);
        tstart = tic; fprintf('%s:(e=%d/%d)[BART] Zeropadding to full size... \n%s\n', datetime, echo, nr_contrasts, command);
        [status_resize,result_resize] = system(command);
        fprintf('%s: done! (%6.4f/%6.4f sec)\n', datetime, toc(tstart), toc(start_time));

        %------------------------------------------------------------------
        % ESPIRiT calibration
        %------------------------------------------------------------------
        % Usage: ecalib [-t f] [-c f] [-k d:d:d] [-r d:d:d] [-m d] [-S] [-W] [-I] [-1]
        %               [-P] [-v f] [-a] [-d d] <kspace> <sensitivities> [<ev-maps>]
        %
        % Estimate coil sensitivities using ESPIRiT calibration.
        % Optionally outputs the eigenvalue maps.
        %
        % -t threshold     This determines the size of the null-space
        % -c crop_value    Crop the sensitivities if the eigenvalue is smaller than {crop_value}
        % -k ksize         kernel size
        % -r cal_size      Limits the size of the calibration region
        % -m maps          Number of maps to compute
        % -S               create maps with smooth transitions (Soft-SENSE)
        % -W               soft-weighting of the singular vectors
        % -I               intensity correction
        % -1               perform only first part of the calibration
        % -P               Do not rotate the phase with respect to the first principal component
        % -v variance      Variance of noise in data
        % -a               Automatically pick thresholds
        % -d level         Debug level
        % -h               help
        %------------------------------------------------------------------
        command = sprintf('%s ecalib -t 0.001 -c 0 -k6:6:6 -r%d:%d:%d -m 1 -d5 %s %s %s', bart_command, cal_size(1), cal_size(2), cal_size(3), kgrid_zpad_file, sens_ecalib_file, ev_maps_file);
        tstart = tic; fprintf('%s:(e=%d/%d)[BART] Estimating coil sensitivities using ESPIRiT calibration... \n%s\n', datetime, echo, nr_contrasts, command);
        [status_ecalib,result_ecalib] = system(command);
        fprintf('%s: done! (%6.4f/%6.4f sec)\n', datetime, toc(tstart), toc(start_time));
    end

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
    command = sprintf('%s nlinv -a %d -b %d -S -d5 -i%d -x %d:%d:%d -t %s %s %s %s', bart_command, noir_conf.a, noir_conf.b, noir_conf.max_iter, cal_size(1), cal_size(2), cal_size(3), traj_file, ksp_file, img_nlinv_calib_file, sens_nlinv_calib_file);
    tstart = tic; fprintf('%s:(e=%d/%d)[BART] Estimating coil sensitivities using NLINV... \n%s\n', datetime, echo, nr_contrasts, command);
    [status_nlinv,result_nlinv] = system(command);
    fprintf('%s: done! (%6.4f/%6.4f sec)\n', datetime, toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % Tranform from image space to k-space (image space <=> k-space)
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
    tstart = tic; fprintf('%s:(e=%d/%d) Transforming from image space to k-space... \n', datetime, echo, nr_contrasts);
    sens_nlinv_calib = readcfl(fullfile(output_path, sprintf('sens_echo%d_nlinv_calib', echo)));
    ksp_sens_nlinv_calib = sens_nlinv_calib;
    clear sens_nlinv_calib;
    for dim = 1:3
        ksp_sens_nlinv_calib = 1 / sqrt(cal_size(dim) * 2) * fftshift(fft(ifftshift(ksp_sens_nlinv_calib, dim), [], dim), dim);
    end
    fprintf('%s: done! (%6.4f/%6.4f sec)\n', datetime, toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % Zeropad to a 2X grid in k-space
    %----------------------------------------------------------------------
    tstart = tic; fprintf('%s:(e=%d/%d) Zeropadding to a 2X grid in k-space... \n', datetime, echo, nr_contrasts);
    idx1_range = (-floor(cal_size(1)*2/2):ceil(cal_size(1)*2/2)-1).' + floor(recon_matrix_size(1)*2/2) + 1;
    idx2_range = (-floor(cal_size(2)*2/2):ceil(cal_size(2)*2/2)-1).' + floor(recon_matrix_size(2)*2/2) + 1;
    idx3_range = (-floor(cal_size(3)*2/2):ceil(cal_size(3)*2/2)-1).' + floor(recon_matrix_size(3)*2/2) + 1;
    ksp_sens_nlinv_calib_zpad = complex(zeros([recon_matrix_size * 2 size(ksp_sens_nlinv_calib,4)], 'single'));
    ksp_sens_nlinv_calib_zpad(idx1_range, idx2_range, idx3_range, :) = ksp_sens_nlinv_calib;
    clear ksp_sens_nlinv_calib;
    fprintf('%s: done! (%6.4f/%6.4f sec)\n', datetime, toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % Tranform from k-space to image space (image space <=> k-space)
    %----------------------------------------------------------------------
    tstart = tic; fprintf('%s:(e=%d/%d) Transforming from k-space to image space... \n', datetime, echo, nr_contrasts);
    sens_nlinv_2Xgrid = ksp_sens_nlinv_calib_zpad;
    clear ksp_sens_nlinv_calib_zpad;
    for dim = 1:3
        sens_nlinv_2Xgrid = sqrt(recon_matrix_size(dim) * 2) * fftshift(ifft(ifftshift(sens_nlinv_2Xgrid, dim), [], dim), dim);
    end
    fprintf('%s: done! (%6.4f/%6.4f sec)\n', datetime, toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % Crop a 2X grid to a 1X grid in image space
    %----------------------------------------------------------------------
    tstart = tic; fprintf('%s:(e=%d/%d) Cropping from a 2X grid to a 1X grid in image space... \n%s\n', datetime, echo, nr_contrasts);
    idx1_range = (-floor(recon_matrix_size(1)/2):ceil(recon_matrix_size(1)/2)-1).' + floor(recon_matrix_size(1)*2/2) + 1;
    idx2_range = (-floor(recon_matrix_size(2)/2):ceil(recon_matrix_size(2)/2)-1).' + floor(recon_matrix_size(2)*2/2) + 1;
    idx3_range = (-floor(recon_matrix_size(3)/2):ceil(recon_matrix_size(3)/2)-1).' + floor(recon_matrix_size(3)*2/2) + 1;
    sens_nlinv_1Xgrid = sens_nlinv_2Xgrid(idx1_range, idx2_range, idx3_range, :);
    clear sens_nlinv_2Xgrid;
    fprintf('%s: done! (%6.4f/%6.4f sec)\n', datetime, toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % Normalize coil sensitivity maps estimated by NLINV
    %----------------------------------------------------------------------
    tstart = tic; fprintf('%s:(e=%d/%d) Normalizing coil sensitivity maps estimated by NLINV... \n', datetime, echo, nr_contrasts);
    sens_nlinv = bsxfun(@rdivide, sens_nlinv_1Xgrid, sqrt(sum(abs(sens_nlinv_1Xgrid).^2,4)));
    writecfl(fullfile(output_path, sens_nlinv_fiename), sens_nlinv);
    clear sens_nlinv_1Xgrid;
    fprintf('%s: done! (%6.4f/%6.4f sec)\n', datetime, toc(tstart), toc(start_time));

    %% Perform coil combination
    %----------------------------------------------------------------------
    % NUFFT reconstruction (ESPRiT)
    %----------------------------------------------------------------------
    cimg_cfl_file = fullfile(output_path, cimg_filename);
    sens_cfl_file = fullfile(output_path, sens_ecalib_fiename);
    if exist(cimg_cfl_file, 'file') && exist(sens_cfl_file, 'file')
        tstart = tic; fprintf('%s:(e=%d/%d) Performing coil combination (ESPIRiT)... \n', datetime, echo, nr_contrasts);
        cimg = readcfl(cimg_cfl_file);
        sens_ecalib = readcfl(sens_cfl_file);
        img_nufft_ecalib = sum(conj(sens_ecalib) .* cimg, 4);
        writecfl(fullfile(output_path, img_nufft_ecalib_filename), img_nufft_ecalib);
        fprintf('%s: done! (%6.4f/%6.4f sec)\n', datetime, toc(tstart), toc(start_time));
    end

    %----------------------------------------------------------------------
    % NUFFT reconstruction (NLINV)
    %----------------------------------------------------------------------
    cimg_cfl_file = fullfile(output_path, cimg_filename);
    sens_cfl_file = fullfile(output_path, sens_nlinv_fiename);
    if exist(cimg_cfl_file, 'file') && exist(sens_cfl_file, 'file')
        tstart = tic; fprintf('%s:(e=%d/%d) Performing coil combination (NLINV)... \n', datetime, echo, nr_contrasts);
        cimg = readcfl(cimg_cfl_file);
        sens_nlinv = readcfl(sens_cfl_file);
        img_nufft_nlinv = sum(conj(sens_nlinv) .* cimg, 4);
        writecfl(fullfile(output_path, img_nufft_nlinv_filename), img_nufft_nlinv);
        fprintf('%s: done! (%6.4f/%6.4f sec)\n', datetime, toc(tstart), toc(start_time));
    end

    %% Save as .dicom files
    %----------------------------------------------------------------------
    % NUFFT reconstruction (ESPRiT)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path, img_nufft_ecalib_filename);
    if exist(cfl_file, 'file')
        recon_type = sprintf('nufft_echo%d', echo);
        dicom_directory = fullfile(output_path, img_nufft_ecalib_filename);
        save_images_as_dicom_files(twix, recon_matrix_size, recon_resolution, recon_type, dicom_directory, cfl_file);
    end

    %----------------------------------------------------------------------
    % NUFFT reconstruction (NLINV)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path, img_nufft_nlinv_filename);
    if exist(cfl_file, 'file')
        recon_type = sprintf('nufft_echo%d', echo);
        dicom_directory = fullfile(output_path, sprintf('%s_dicom%d', img_nufft_nlinv_filename, dicom_matrix_size(1)));
        save_images_as_dicom_files(twix, recon_matrix_size, recon_resolution, recon_type, dicom_directory, cfl_file);
    end
end
