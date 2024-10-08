% demo_step6_bh_pics_recon.m
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
recon_matrix_size   = json.recon_conf.recon_matrix_size;   % reconstruction matrix size
dicom_matrix_size   = json.recon_conf.dicom_matrix_size;   % DICOM matrix size
recon_interp_factor = json.recon_conf.recon_interp_factor; % reconstruction interpolation factor: recon_resolution = encoded_resolution / recon_interp_factor;
cal_size            = json.recon_conf.cal_size;            % size of the calibration region
lambda              = json.recon_conf.lambda;              % l1 regularization parameter for pics in BART
max_iter            = json.recon_conf.max_iter;            % max. number of iterations for pics in BART
kappa               = json.recon_conf.kappa;               % kappa: kappa=0 => W=I

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

%% Perform PICS reconstruction per echo
for eco = 1:nr_echoes

    %% Define filenames
    traj_filename     = sprintf('traj_eco%d_%dx%dx%d', eco, traj_scale_factor(1), traj_scale_factor(2), traj_scale_factor(3));
    weights_filename  = sprintf('weights_eco%d_effMtx%d_kappa%g', eco, N1, kappa);
    ksp_filename      = sprintf('ksp_eco%d', eco);
    sens_fiename      = sprintf('sens_eco%d_%dx%dx%d', eco, N1, N2, N3);
    img_pics_filename = sprintf('img_pics_eco%d_kappa%g_l%g_i%d_%dx%dx%d', eco, kappa, lambda, max_iter, N1, N2, N3);

    traj_file     = strcat(bart_output_path, traj_filename);
    weights_file  = strcat(bart_output_path, weights_filename);
    ksp_file      = strcat(bart_output_path, ksp_filename);
    sens_file     = strcat(bart_output_path, sens_fiename);
    img_pics_file = strcat(bart_output_path, img_pics_filename);

    %% Read a .cfl file
    cfl_file = fullfile(output_path, sprintf('dcf_eco%d_effMtx%d', eco, N1));
    tstart = tic; fprintf('%s:(ECO=%d/%d) Reading a .cfl file: %s... ', datetime, eco, nr_echoes, cfl_file);
    dcf = readcfl(cfl_file);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %% Calculate W
    %----------------------------------------------------------------------
    % argmin ||sqrt(W) * (F * S * m - y)||_2^2 + Phi(m)
    % pics uses a normal operator to solve the least-squares problem,
    % we need the square root of our filter here:
    % weights = sqrt(W), where W = dcf^kappa
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path, weights_filename);
    tstart = tic; fprintf('%s:(ECO=%d/%d) Writing a .cfl file: %s... ', datetime, eco, nr_echoes, cfl_file);
    writecfl(cfl_file, sqrt((dcf).^kappa));
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
    clear dcf;

    %% Parallel imaging and compressed sensing reconstruction
    %----------------------------------------------------------------------
    % Usage: pics [-l ...] [-r f] [-R ...] [-c] [-s f] [-i d] [-t <file>] [-n]
    %             [-N] [-g] [--gpu-gridding] [-G d] [-p <file>] [--precond] [-I]
    %             [--fista] [--pridu] [-b d] [-e] [-W <file>] [-d d] [-u f]
    %             [-C d] [-f f] [-m] [-w f] [-S] [-L d] [-K] [-B <file>] [-P f]
    %             [-a] [-M] [-U,--lowmem] [--no-toeplitz] [--psf_export <file>]
    %             [--psf_import <file>] [--wavelet <string>]
    %             <kspace> <sensitivities> <output>
    %
    % Parallel-imaging compressed-sensing reconstruction.
    %
    % -l1/-l2              toggle l1-wavelet or l2 regularization.
    % -r lambda            regularization parameter
    % -R <T>:A:B:C         generalized regularization options (-Rh for help)
    % -c                   real-value constraint
    % -s step              iteration stepsize
    % -i iter              max. number of iterations
    % -t file              k-space trajectory
    % -n                   disable random wavelet cycle spinning
    % -N                   do fully overlapping LLR blocks
    % -g                   use GPU
    % --gpu-gridding       use GPU for gridding
    % -G gpun              use GPU device gpun
    % -p file              pattern or weights
    % --precond            interprete weights as preconditioner
    % -I                   select IST
    % --fista              select FISTA
    % --pridu              select primal dual
    % -b blk               Lowrank block size
    % -e                   Scale stepsize based on max. eigenvalue
    % -W <img>             Warm start with <img>
    % -d level             Debug level
    % -u rho               ADMM rho
    % -C iter              ADMM max. CG iterations
    % -f rfov              restrict FOV
    % -m                   select ADMM
    % -w f                 inverse scaling of the data
    % -S                   re-scale the image after reconstruction
    % -L flags             batch-mode
    % -K                   randshift for NUFFT
    % -B file              temporal (or other) basis
    % -P eps               Basis Pursuit formulation, || y- Ax ||_2 <= eps
    % -a                   select Primal Dual
    % -M                   Simultaneous Multi-Slice reconstruction
    % -U,--lowmem          Use low-mem mode of the nuFFT
    % --no-toeplitz        Turn off Toeplitz mode of nuFFT
    % --psf_export file    Export PSF to file
    % --psf_import file    Import PSF from file
    % --wavelet name       wavelet type (haar,dau2,cdf44)
    % -h                   help
    %
    %----------------------------------------------------------------------
    % Generalized regularization options (experimental)
    %
    % -R <T>:A:B:C    <T> is regularization type (single letter),
    %                 A is transform flags, B is joint threshold flags,
    %                 and C is regularization value. Specify any number
    %                 of regularization terms.
    %
    % -R Q:C          l2-norm in image domain
    % -R I:B:C        l1-norm in image domain
    % -R W:A:B:C      l1-wavelet
    % -R N:A:B:C      Normalized Iterative Hard Thresholding (NIHT), image domain
    %                 C is an integer percentage, i.e. from 0-100
    % -R H:A:B:C      NIHT, wavelet domain
    % -R F:A:B:C      l1-Fourier
    % -R T:A:B:C      total variation
    % -R T:7:0:.01    3D isotropic total variation with 0.01 regularization.
    % -R G:A:B:C      total generalized variation
    % -R C:A:B:C      infimal convolution TV
    % -R L:7:7:.02    Locally low rank with spatial decimation and 0.02 regularization.
    % -R M:7:7:.03    Multi-scale low rank with spatial decimation and 0.03 regularization.
    % -R TF:{graph_path}:lambda       TensorFlow loss
    %----------------------------------------------------------------------
    command = sprintf('%s pics -e -S -d5 -R W:7:0:%g -i %d -t %s -p %s --fista --wavelet dau2 %s %s %s', bart_command, lambda, max_iter, traj_file, weights_file, ksp_file, sens_file, img_pics_file);
    tstart = tic; fprintf('%s:(ECO=%d/%d)[BART] parallel imaging and compressed sensing reconstruction:\n%s\n', datetime, eco, nr_echoes, command);
    [status_pics,result_pics] = system(command);
    fprintf('%s: done! (%6.4f/%6.4f sec)\n', datetime, toc(tstart), toc(start_time));

    %% Save the output of the pics command as a .txt file
    txt_file = fullfile(output_path, sprintf('%s_debug.txt', img_pics_filename));
    [fid,message] = fopen(txt_file, 'w');
    fwrite(fid, result_pics, 'char');
    fclose(fid);

    %% Save as a .dicom file
    recon_type = sprintf('eco%d', eco);
    cfl_file = fullfile(output_path, img_pics_filename);
    dicom_directory = fullfile(output_path, sprintf('%s_dicom%d', img_pics_filename, dicom_matrix_size(1)));
    save_images_as_dicom_files(twix, dicom_matrix_size, recon_resolution, recon_type, dicom_directory, cfl_file);
end

%% Combine both echoes
img_pics_combined = complex(zeros(N1, N2, N3, 'single'));
for eco = 1:nr_echoes
    img_pics_filename = sprintf('img_pics_eco%d_kappa%g_l%g_i%d_%dx%dx%d', eco, kappa, lambda, max_iter, N1, N2, N3);
    img_pics_combined = img_pics_combined + readcfl(fullfile(output_path, img_pics_filename));
end

%% Define the full path of a .cfl file without file extension
img_pics_combined_filename = sprintf('img_pics_combined_kappa%g_l%g_i%d_%dx%dx%d', kappa, lambda, max_iter, N1, N2, N3);

%% Save as a .cfl file
cfl_file = fullfile(output_path, img_pics_combined_filename);
tstart = tic; fprintf('%s: Writing a .cfl file: %s... ', datetime, cfl_file);
writecfl(cfl_file, img_pics_combined);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%% Save as a .dicom file
recon_type = 'combined';
dicom_directory = fullfile(output_path, sprintf('%s_dicom%d',img_pics_combined_filename, dicom_matrix_size(1)));
save_images_as_dicom_files(twix, dicom_matrix_size, recon_resolution, recon_type, dicom_directory, cfl_file);
