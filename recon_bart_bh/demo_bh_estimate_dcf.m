% demo_bh_estimate_dcf.m
% Written by Nam Gyun Lee
% Email: namgyunl@usc.edu, ggang56@gmail.com (preferred)
% Started: 08/02/2023, Last modified: 08/04/2023

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
    output_path = strrep(json.output_path, '/', '\');
else
    output_path = json.output_path;
end

%--------------------------------------------------------------------------
% Reconstruction parameters
%--------------------------------------------------------------------------
recon_matrix_size = json.recon_parameters.recon_matrix_size.'; % reconstruction matrix size

%--------------------------------------------------------------------------
% Sequence parameters
%--------------------------------------------------------------------------
nr_contrasts = json.seq.nr_contrasts;

%% Estimate density compensation factors (1 x Nk x NiNs)
dir_info = dir(fullfile(output_path, 'traj*.cfl'));
nr_files = length(dir_info);

for idx = 1:nr_files
    %----------------------------------------------------------------------
    % Get the respiratory state number and echo number
    %----------------------------------------------------------------------
    under_loc = strfind(dir_info(idx).name, '_');
    b = str2double(dir_info(idx).name(under_loc(1)-1));
    echo = str2double(dir_info(idx).name(under_loc(2)-1));

    %----------------------------------------------------------------------
    % Read a .cfl file
    %----------------------------------------------------------------------
    dot_loc = strfind(dir_info(idx).name, '.');
    cfl_file = fullfile(dir_info(idx).folder, dir_info(idx).name(1:dot_loc-1));
    tstart = tic; fprintf('%s: (%2d/%2d) Reading a .cfl file: %s... \n', datetime, idx, nr_files, cfl_file);
    coords = double(readcfl(cfl_file)); % 3 x Nk/2 x NiNs
    fprintf('%s: done! (%6.4f/%6.4f sec)\n', datetime, toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % Estimate density compensation factors
    %----------------------------------------------------------------------
    tstart = tic; fprintf('%s:(e=%d/%d) Estimating density compensation factors using sdc3_MAT.c... \n', datetime, echo, nr_contrasts);
    numIter = 25;                   % number of iterations
    effMtx  = recon_matrix_size(1); % the length of one side of the grid matrix
    verbose = 0;                    % 1:verbose 0:quiet
    osf     = 2;                    % the grid oversample factor
    DCF = sdc3_MAT(coords, numIter, effMtx, verbose, osf);
    dcf = reshape(DCF / max(DCF(:)), [1 size(coords, [2 3])]);
    fprintf('%s: done! (%6.4f/%6.4f sec)\n', datetime, toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % Write a .cfl file
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path, sprintf('dcf_echo%d', echo));
    tstart = tic; fprintf('%s: Writing a .cfl file: %s... \n', datetime, cfl_file);
    writecfl(cfl_file, dcf);
    fprintf('%s: done! (%6.4f/%6.4f sec)\n', datetime, toc(tstart), toc(start_time));
end
