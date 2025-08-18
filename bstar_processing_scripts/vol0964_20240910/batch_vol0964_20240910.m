% vol0964_20240910.m
% Written by Nam Gyun Lee
% Email: namgyunl@usc.edu, ggang56@gmail.com (preferred)
% Started: 09/12/2024, Last modified: 08/18/2025

%% Clean slate
close all; clear all; clc;

%% Set source paths
if ispc
    package_path = 'D:\replication_bstar';
    ismrmrd_path = 'D:\ismrmrd';
    pulseq_path  = 'D:\pulseq';
else
    package_path = '/server/sdata/nlee/replication_bstar';
    ismrmrd_path = '/server/sdata/nlee/ismrmrd';
    pulseq_path  = '/server/home/nlee/pulseq';   
end

%% Add source paths to search path
addpath(genpath(package_path));
addpath(genpath(ismrmrd_path));
addpath(genpath(pulseq_path));

%% Define a .json file
json_files{1} = 'D:\replication_bstar\data\vol0964_20240910\meas_MID01338_FID22747_bstar_exhale_43s.json';

%% Perform image reconstruction
nr_json_files = length(json_files);

for json_nr = 1:nr_json_files
    json_file = json_files{json_nr};

    %% Step 1: Prepare k-space data as a .cfl file
    demo_step1_bh_prepare_ksp;

    %% Step 2: Prepare k-space trajectories
    demo_step2_bh_prepare_traj;

    %% Step 3: Estimate density compensation factors
    demo_step3_bh_estimate_dcf;

    %% Step 4: Perform gridding recon
    demo_step4_bh_gridding_recon;

    %% Step 5: Estimate coil sensitivity maps
    demo_step5_bh_estimate_csm;

    %% Step 6: Perform PICS reconstruction
    demo_step6_bh_pics_recon;
end
