% batch_pulseq_bstar_acr_phantom_20230614.m
% Written by Nam Gyun Lee
% Email: namgyunl@usc.edu, ggang56@gmail.com (preferred)
% Started: 06/12/2023, Last modified: 06/12/2023

%% Clean slate
close all; clear all; clc;

%% Set source directories
package_directory = 'D:\replication_bstar';
ismrmrd_directory = 'D:\ismrmrd';
pulseq_directory = 'D:\pulseq\pulseq';

%% Add source directories to search path
addpath(genpath(package_directory));
addpath(genpath(ismrmrd_directory));
addpath(genpath(pulseq_directory));

%% Define a .json file
json_files{1} = 'D:\replicability_bstar\recon_bart_bh\pulseq_bstar_acr_phantom_20230614\meas_MID00184_FID66316_pulseq_bstar_1_38ms_1_61mm_i89_31k_SP.json';
json_files{2} = 'D:\replicability_bstar\recon_bart_bh\pulseq_bstar_acr_phantom_20230614\meas_MID00185_FID66317_pulseq_bstar_1_38ms_1_61mm_i88_31k_SP.json';
json_files{3} = 'D:\replicability_bstar\recon_bart_bh\pulseq_bstar_acr_phantom_20230614\meas_MID00186_FID66318_pulseq_bstar_1_38ms_1_61mm_i89_31k_WASP.json';
json_files{4} = 'D:\replicability_bstar\recon_bart_bh\pulseq_bstar_acr_phantom_20230614\meas_MID00187_FID66319_pulseq_bstar_1_38ms_1_61mm_i88_31k_WASP.json';

%% Perform bSTAR reconstruction
nr_json_files = length(json_files);

for idx = 1:nr_json_files
    json_file = json_files{idx};

    %----------------------------------------------------------------------
    % Step 1: Prepare k-space data as a .cfl file
    %----------------------------------------------------------------------
    demo_bh_prepare_ksp;

    %----------------------------------------------------------------------
    % Step 2: Prepare k-space trajectories
    %----------------------------------------------------------------------
    demo_bh_prepare_traj;

    %----------------------------------------------------------------------
    % Step 3: Estimate density compensation factors
    %----------------------------------------------------------------------
    demo_bh_estimate_dcf;

    %----------------------------------------------------------------------
    % Step 4: Gridding recon
    %----------------------------------------------------------------------
    %demo_bh_gridding_recon;

    %----------------------------------------------------------------------
    % Step 5: Estimate CSM
    %----------------------------------------------------------------------
    demo_bh_estimate_csm;

    %----------------------------------------------------------------------
    % Step 6: PICS reconstruction
    %----------------------------------------------------------------------
    demo_bh_pics_recon;
end
