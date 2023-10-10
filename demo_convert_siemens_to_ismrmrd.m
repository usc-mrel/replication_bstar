% demo_convert_siemens_to_ismrmrd.m
% Written by Nam Gyun Lee
% Email: namgyunl@usc.edu, ggang56@gmail.com (preferred)
% Started: 08/16/2022, Last modified: 08/16/2022

%% Clean slate
close all; clear all; clc;

%% Define directory containing TWIX files
data_directory = 'D:\data_pulseq_bstar\pulseq_bstar_acr_phantom_20230614\raw';

%% Get directory information
dir_info = dir(fullfile(data_directory, "*.dat"));

%% Convert TWIX to ISMRMRD format
start_time = tic;
cd(data_directory);
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
