% batch_fb_vol_20240428.m
% Written by Nam Gyun Lee
% Email: namgyunl@usc.edu, ggang56@gmail.com (preferred)
% Started: 04/28/2024, Last modified: 04/28/2024

%% Clean slate
close all; clear all; clc;

%% Set source paths
if ispc
    package_path = 'D:\replication_bstar_v2';
    ismrmrd_path = 'D:\ismrmrd';
    pulseq_path  = 'D:\pulseq';
else
    package_path = '/server/sdata/nlee/replication_bstar_v2';
    ismrmrd_path = '/server/sdata/nlee/cartesian_maxgirf_epi/thirdparty/ismrmrd';
    pulseq_path  = '/server/home/nlee/pulseq';   
end

%% Add source paths to search path
addpath(genpath(package_path));
addpath(genpath(ismrmrd_path));
addpath(genpath(pulseq_path));

%% Define a .json file
json_files{1} = '/server/sdata/nlee/replication_bstar_v2/bstar_processing_scripts/vol_20240428/meas_MID00572_FID14194_pulseq_bstar_0_9mm_fb_server.json';

%% Perform image reconstruction
nr_json_files = length(json_files);

for json_nr = 1:nr_json_files
    json_file = json_files{json_nr};

    %% Step 1: Prepare k-space data as a .cfl file
    demo_step1_fb_prepare_ksp;

    %% Step 2: Estimate a respiratory signal for respiratory data binning
    demo_step2_fb_estimate_resp;

    %% Step 3: Prepare k-space trajectories per respiratory state
    demo_step3_fb_prepare_btraj;

    %% Step 4: Estimate density compensation factors
    demo_step4_fb_estimate_dcf;

    %% Step 5: Prepare k-space data per respiratory state
    demo_step5_fb_prepare_bksp;

    %% Step 6: Perform gridding recon
    %demo_bh_gridding_recon;

    %% Step 7: Estimate coil sensitivity maps
    demo_step7_fb_estimate_csm;

    %% Step 8: Perform PICS reconstruction
    demo_step8_fb_pics_recon;
end
