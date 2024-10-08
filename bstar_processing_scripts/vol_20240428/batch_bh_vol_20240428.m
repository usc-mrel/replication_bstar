% batch_vol_20240428.m
% Written by Nam Gyun Lee
% Email: namgyunl@usc.edu, ggang56@gmail.com (preferred)
% Started: 04/28/2024, Last modified: 04/28/2024

%% Clean slate
close all; clear all; clc;

%% Set source paths
package_path = 'D:\replication_bstar_v2';
ismrmrd_path = 'D:\ismrmrd';
pulseq_path = 'D:\pulseq';

%% Add source paths to search path
addpath(genpath(package_path));
addpath(genpath(ismrmrd_path));
addpath(genpath(pulseq_path));

%% Define a .json file
json_files{1} = 'D:\replication_bstar_v2\bstar_processing_scripts\vol_20240428\meas_MID00565_FID14187_pulseq_bstar_breathhold.json';

%% Perform image reconstruction
nr_json_files = length(json_files);

for json_nr = 1:nr_json_files
    json_file = json_files{json_nr};

    %% Step 1: Prepare k-space data as a .cfl file
    demo_bh_prepare_ksp;

    return
    
    %% Step 2: Prepare k-space trajectories
    demo_bh_prepare_traj;

    %% Step 3: Estimate density compensation factors
    demo_bh_estimate_dcf;

    %% Step 4: Perform gridding recon
    demo_bh_gridding_recon;

    %% Step 5: Estimate coil sensitivity maps
    demo_bh_estimate_csm;

    %% Step 6: Perform PICS reconstruction
    demo_bh_pics_recon;
end
