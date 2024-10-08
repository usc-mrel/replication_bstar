% demo_step4_fb_estimate_dcf.m
% Written by Nam Gyun Lee
% Email: namgyunl@usc.edu, ggang56@gmail.com (preferred)
% Started: 04/26/2024, Last modified: 09/18/2024

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
    output_parent_path = strrep(json.output_path, '/', '\');
else
    output_parent_path = json.output_path;
end

%--------------------------------------------------------------------------
% Reconstruction parameters
%--------------------------------------------------------------------------
recon_matrix_size = json.recon_conf.recon_matrix_size.'; % reconstruction matrix size
B                 = json.recon_conf.B;                   % number of frames (respiratory motion states)

%--------------------------------------------------------------------------
% Define the number of echoes
%--------------------------------------------------------------------------
nr_echoes = 2;

%% Define an output path
output_path = fullfile(output_parent_path, sprintf('%d_motion_states', B));

%% Estimate density compensation factors (1 x adc_samples_per_echo x nr_radial_views)
for b = 1%:B
    for eco = 1:nr_echoes
        %% Read a .cfl file
        cfl_file = fullfile(output_path, sprintf('btraj%d_eco%d_unscaled', b, eco));
        tstart = tic; fprintf('%s:(B=%2d/%2d)(ECO=%d/%d) Reading a .cfl file: %s... ', datetime, b, B, eco, nr_echoes, cfl_file);
        coords = double(readcfl(cfl_file)); % 3 x adc_samples_per_echo x nr_radial_views
        fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

        %% Estimate density compensation factors
        tstart = tic; fprintf('%s:(B=%2d/%2d)(ECO=%d/%d) Estimating density compensation factors using sdc3_MAT.c... ', datetime, b, B, eco, nr_echoes);
        numIter = 25;                   % number of iterations
        effMtx  = recon_matrix_size(1); % the length of one side of the grid matrix
        verbose = 0;                    % 1:verbose 0:quiet
        osf     = 2;                    % the grid oversample factor
        DCF = sdc3_MAT(coords, numIter, effMtx, verbose, osf); % coords in [-0.5 0.5]
        dcf = reshape(DCF / max(DCF(:)), [1 size(coords, [2 3])]);
        fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

        %% Write a .cfl file
        cfl_file = fullfile(output_path, sprintf('bdcf%d_eco%d_effMtx%d', b, eco, recon_matrix_size(1)));
        tstart = tic; fprintf('%s:(B=%2d/%2d)(ECO=%d/%d) Writing a .cfl file: %s... ', datetime, b, B, eco, nr_echoes, cfl_file);
        writecfl(cfl_file, dcf);
        fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
    end
end
