% demo_step3_bh_estimate_dcf.m
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
    output_path = strrep(json.output_path, '/', '\');
else
    output_path = json.output_path;
end

%--------------------------------------------------------------------------
% Reconstruction parameters
%--------------------------------------------------------------------------
recon_matrix_size = json.recon_conf.recon_matrix_size; % reconstruction matrix size

%--------------------------------------------------------------------------
% Define the number of echoes
%--------------------------------------------------------------------------
nr_echoes = 2;

%% Estimate density compensation factors (1 x adc_samples_per_echo x nr_radial_views)
for eco = 1:nr_echoes
    %% Read a .cfl file
    cfl_file = fullfile(output_path, sprintf('traj_unscaled_eco%d', eco));
    tstart = tic; fprintf('%s:(ECO=%d/%d) Reading a .cfl file: %s... ', datetime, eco, nr_echoes, cfl_file);
    coords = double(readcfl(cfl_file)); % 3 x adc_samples_per_echo x nr_radial_views
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %% Estimate density compensation factors
    tstart = tic; fprintf('%s:(ECO=%d/%d) Estimating density compensation factors using sdc3_MAT.c... ', datetime, eco, nr_echoes);
    numIter = 25;                   % number of iterations
    effMtx  = recon_matrix_size(1); % the length of one side of the grid matrix
    verbose = 0;                    % 1:verbose 0:quiet
    osf     = 2;                    % the grid oversample factor
    DCF = sdc3_MAT(coords, numIter, effMtx, verbose, osf); % coords in [-0.5 0.5]
    dcf = reshape(DCF / max(DCF(:)), [1 size(coords, [2 3])]);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %% Write a .cfl file
    cfl_file = fullfile(output_path, sprintf('dcf_eco%d_effMtx%d', eco, recon_matrix_size(1)));
    tstart = tic; fprintf('%s:(ECO=%d/%d) Writing a .cfl file: %s... ', datetime, eco, nr_echoes, cfl_file);
    writecfl(cfl_file, dcf);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));    
end
