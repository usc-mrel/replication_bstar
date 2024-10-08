% demo_step2_fb_estimate_resp.m
% Written by Nam Gyun Lee
% Email: namgyunl@usc.edu, ggang56@gmail.com (preferred)
% Started: 09/17/2024, Last modified: 09/18/2024

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
    ismrmrd_data_file = strrep(json.ismrmrd_data_file, '/', '\');
    ismrmrd_traj_file = strrep(json.ismrmrd_traj_file, '/', '\');
    output_path       = strrep(json.output_path, '/', '\');
else
    ismrmrd_data_file = json.ismrmrd_data_file;
    ismrmrd_traj_file = json.ismrmrd_traj_file;
    output_path       = json.output_path;
end

%--------------------------------------------------------------------------
% Respiratory filtering parameters
%--------------------------------------------------------------------------
n  = json.resp_conf.n;
fl = json.resp_conf.fl;
fh = json.resp_conf.fh;
fw = json.resp_conf.fw;

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

%% Get imaging parameters from an XML header
traj_header = ismrmrd.xml.deserialize(traj_dset.readxml);
data_header = ismrmrd.xml.deserialize(data_dset.readxml);

%% Get sequence parameters
%--------------------------------------------------------------------------
% self_navigation_flag
%--------------------------------------------------------------------------
self_navigation_flag = traj_header.encoding.trajectoryDescription.userParameterLong(4).value;

%--------------------------------------------------------------------------
% TR
%--------------------------------------------------------------------------
TR = data_header.sequenceParameters.TR * 1e-3; % [msec] * [sec/1e3msec] => [sec]

%% Read a .cfl file
if self_navigation_flag
    cfl_file = fullfile(output_path, 'ksp_nav');
else
    cfl_file = fullfile(output_path, 'k0');
end

tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
ksp = readcfl(cfl_file);
fprintf('%s: done! (%6.4f/%6.4f sec)\n', datetime, toc(tstart), toc(start_time));

%% Get the size of "ksp"
[discard_pre,nr_radial_views,nr_channels] = size(ksp, [2 3 4]);

%% Get the DC signals (nr_radial_views x nr_channels)
index_range = (7:discard_pre-5).';
dc = reshape(mean(ksp(1,index_range,:,:)), [nr_radial_views nr_channels]);

%% Estimate a respiratory signal from DC
tstart = tic; fprintf('%s: Estimating a respiratory signal from DC... \n', datetime);
resp = estimate_resp(dc, TR, n, fl, fh, fw, self_navigation_flag);
fprintf('%s: done! (%6.4f/%6.4f sec)\n', datetime, toc(tstart), toc(start_time));

%% Write a .cfl file
%--------------------------------------------------------------------------
% resp (nr_radial_views x 1)
%--------------------------------------------------------------------------
cfl_file = fullfile(output_path, sprintf('resp_n%d_fl%4.2f_fh%4.2f_fw%4.2f', n, fl, fh, fw));
tstart = tic; fprintf('%s: Writing a .cfl file: %s... ', datetime, cfl_file);
writecfl(cfl_file, resp);
fprintf('%s: done! (%6.4f/%6.4f sec)\n', datetime, toc(tstart), toc(start_time));

%% Calculate the spectrum of a respiratory signal
resp_freq = 1 / sqrt(nr_radial_views) * fftshift(fft(resp, [], 1), 1);

%% Display a respiratory signal and its spectrum
t = (0:nr_radial_views-1).' * TR; % time axis [sec]
df = 1 / (nr_radial_views * TR); % frequency resolution [Hz]
f = (-floor(nr_radial_views/2):ceil(nr_radial_views/2)-1).' * df; % frequency axis [Hz]

figure('Color', 'w', 'Position', [1 218 1918 773]);
ax1 = subplot(2,1,1);
plot(t, resp);
grid on;
grid minor;
xlabel('Time [sec]');
ylabel('Amplitude [a.u.]');
title(sprintf('Respiratory signal, n = %d, fl = %4.2f, fh = %4.2f, fw = %4.2f', n, fl, fh, fw));

subplot(2,1,2);
plot(f, abs(resp_freq));
grid on;
grid minor;
xlim([-5 5]);
ylim([0 0.5e2]);
xlabel('Frequency [Hz]');
ylabel('Magnitude [a.u.]');
title('Spectrum');

export_fig(fullfile(output_path, sprintf('resp_signal_spectrum_n%d_fl%4.2f_fh%4.2f_fw%4.2f_full', n, fl, fh, fw)), '-r300', '-tif');
xlim(ax1, [0 350]);
export_fig(fullfile(output_path, sprintf('resp_signal_spectrum_n%d_fl%4.2f_fh%4.2f_fw%4.2f', n, fl, fh, fw)), '-r300', '-tif');
