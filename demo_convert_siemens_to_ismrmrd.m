% demo_convert_siemens_to_ismrmrd.m
% Written by Nam Gyun Lee
% Email: namgyunl@usc.edu, ggang56@gmail.com (preferred)
% Started: 08/16/2022, Last modified: 04/22/2024

%% Clean slate
close all; clear all; clc;

%% Define a data path containing a .dat file (Siemens TWIX file)
data_path = 'D:\replication_bstar\data\vol0964_20240910\raw';

%% Start a stopwatch timer
start_time = tic;

%% Get directory information
dir_info = dir(fullfile(data_path, "*.dat"));

%% Convert TWIX format to ISMRMRD format
cd(data_path);
nr_twix_files = length(dir_info);

for idx = 1:nr_twix_files
    twix_filename = dir_info(idx).name;
    dot_loc = strfind(twix_filename, '.');

    %----------------------------------------------------------------------
    % Convert noise data
    %----------------------------------------------------------------------
    linux_command1 = sprintf('wsl siemens_to_ismrmrd -f %s -z 1 -o noise_%s.h5', twix_filename, twix_filename(1:dot_loc-1));
    tic; fprintf('(%2d/%2d): Running %s... ', idx, nr_twix_files, linux_command1);
    [status1,result1] = system(linux_command1);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc, toc(start_time));

    %----------------------------------------------------------------------
    % Convert imaging data
    %----------------------------------------------------------------------
    linux_command2 = sprintf('wsl siemens_to_ismrmrd -f %s -z 2 -o %s.h5', twix_filename, twix_filename(1:dot_loc-1));
    tic; fprintf('(%2d/%2d): Running %s... ', idx, nr_twix_files, linux_command2);
    [status2,result2] = system(linux_command2);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc, toc(start_time));
end
