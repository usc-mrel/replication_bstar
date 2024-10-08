% demo_step5_fb_prepare_bksp.m
% Written by Nam Gyun Lee
% Email: namgyunl@usc.edu, ggang56@gmail.com (preferred)
% Started: 08/02/2023, Last modified: 09/18/2024

%% Clean slate
close all; clearvars -except json_number nr_json_files json_files json_file; clc;

%% Start a stopwatch timer
start_time = tic;

%% Read a .json file
tstart = tic; fprintf('%s: Reading a .json file: %s... \n', datetime, json_file);
fid = fopen(json_file); 
json_txt = fread(fid, [1 inf], 'char=>char'); 
fclose(fid);
json = jsondecode(json_txt);
fprintf('%s: done! (%6.4f/%6.4f sec)\n', datetime, toc(tstart), toc(start_time));

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
B = json.recon_conf.B; % number of frames (respiratory motion states)

%--------------------------------------------------------------------------
% Respiratory filtering parameters
%--------------------------------------------------------------------------
n  = json.resp_conf.n;
fl = json.resp_conf.fl;
fh = json.resp_conf.fh;
fw = json.resp_conf.fw;

%--------------------------------------------------------------------------
% Define the number of echoes
%--------------------------------------------------------------------------
nr_echoes = 2;

%% Define an output path
output_path = fullfile(output_parent_path, sprintf('%d_motion_states', B));

%% Read a .cfl file
%--------------------------------------------------------------------------
% ksp (1 x adc_samples x nr_radial_views x nr_channels)
%--------------------------------------------------------------------------
cfl_file = fullfile(output_parent_path, 'ksp');
tstart = tic; fprintf('%s: Reading a .cfl file: %s... \n', datetime, cfl_file);
ksp = readcfl(cfl_file);
[adc_samples,nr_radial_views,nr_channels] = size(ksp, [2 3 4]);
fprintf('%s: done! (%6.4f/%6.4f sec)\n', datetime, toc(tstart), toc(start_time));

%--------------------------------------------------------------------------
% resp (nr_radial_views x 1)
%--------------------------------------------------------------------------
cfl_file = fullfile(output_parent_path, sprintf('resp_n%d_fl%4.2f_fh%4.2f_fw%4.2f', n, fl, fh, fw));
tstart = tic; fprintf('%s: Reading a .cfl file: %s... \n', datetime, cfl_file);
resp = real(readcfl(cfl_file));
fprintf('%s: done! (%6.4f/%6.4f sec)\n', datetime, toc(tstart), toc(start_time));

%% Perform data binning for respiratory motion
margin = 10;
bins = prctile(resp, linspace(0 + margin, 100 - margin, B + 1).');
bidx = cell(B,1);
for b = 1:B
    idx = find((resp >= bins(b)) & (resp < bins(b + 1)));
    bidx{b} = idx;
end

%% Write a .cfl file
adc_samples_per_echo = adc_samples / nr_echoes;

for b = 1%:B
    for eco = 1:nr_echoes
        %------------------------------------------------------------------
        % Calculate the index range of samples per echo
        %------------------------------------------------------------------
        echo_range = (1:adc_samples_per_echo) + (eco - 1) * adc_samples_per_echo;

        %------------------------------------------------------------------
        % bksp (1 x adc_samples_per_echo x # of TRs per motion state x nr_channels)
        %------------------------------------------------------------------
        cfl_file = fullfile(output_path, sprintf('bksp%d_eco%d', b, eco));
        tstart = tic; fprintf('%s:(B=%2d/%2d)(ECO=%d/%d) Writing a .cfl file: %s... \n', datetime, b, B, eco, nr_echoes, cfl_file);
        writecfl(cfl_file, ksp(:,echo_range,bidx{b},:));
        fprintf('%s: done! (%6.4f/%6.4f sec)\n', datetime, toc(tstart), toc(start_time));
    end
end
