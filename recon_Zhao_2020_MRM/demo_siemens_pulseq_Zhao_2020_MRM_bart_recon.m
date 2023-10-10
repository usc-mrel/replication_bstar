% demo_siemens_pulseq_Zhao_2020_MRM_bart_recon.m
% Written by Nam Gyun Lee
% Email: namgyunl@usc.edu, ggang56@gmail.com (preferred)
% Started: 10/02/2022, Last modified: 10/02/2022

%% Clean slate
close all; clear all; clc;

%% Set source directories
src_directory = 'D:\lowfield_bstar\src';
thirdparty_directory = 'D:\lowfield_bstar\pulseq_bstar\recon_Zhao_2020_MRM\thirdparty';
ismrmrd_directory = 'D:\ismrmrd';
pulseq_directory = 'D:\pulseq\pulseq';

%% Add source directories to search path
addpath(genpath(src_directory));
addpath(genpath(thirdparty_directory));
addpath(genpath(ismrmrd_directory));
addpath(genpath(pulseq_directory));

%% Add '$(TOOLBOX_PATH)/matlab' to the library path
bart_directory = '\\wsl$\Ubuntu-20.04\home\image\bart\matlab';
bart_path = '/home/image/bart';
setenv('TOOLBOX_PATH', bart_path);
addpath(bart_directory);

%% Define data fullpath
% siemens_twix_path  = 'D:\data_pulseq_bstar\pulseq_Zhao_ball1_20221002\trajectory\meas_MID00367_FID46436_pulseq_Zhao_v2_100kHz_RT140us_TT520us_delay500us_amp6_64x64.dat';
% ismrmrd_data_path  = 'D:\data_pulseq_bstar\pulseq_Zhao_ball1_20221002\trajectory\meas_MID00367_FID46436_pulseq_Zhao_v2_100kHz_RT140us_TT520us_delay500us_amp6_64x64.h5';
% ismrmrd_noise_path = 'D:\data_pulseq_bstar\pulseq_Zhao_ball1_20221002\trajectory\noise_meas_MID00367_FID46436_pulseq_Zhao_v2_100kHz_RT140us_TT520us_delay500us_amp6_64x64.h5';
% seq_path           = 'D:\data_pulseq_bstar\pulseq_Zhao_ball1_20221002\seq\Zhao_v2_100kHz_RT140us_TT520us_delay500us_amp6_Smax50_slice5.00mm_64x64.seq';
% output_filename    = sprintf('pulseq_Zhao_ball1_20221002');

output_filename = sprintf('pulseq_Zhao_ball2_20221002');
seq_scope_path  = 'D:\data_pulseq_bstar\pulseq_Zhao_ball2_20221002\seq\Zhao_v2_100kHz_RT140us_TT520us_delay500us_amp6_Smax50_slice5.00mm_64x64.seq';
seq_bstar_path  = 'D:\data_pulseq_bstar\pulseq_bstar_vol422_20220830\seq\bstar_ecg0_TR1.38ms_1.61mm_b1929_rf200_i89_40k_FA25.seq';

%--------------------------------------------------------------------------
% Coronal, R >> L (readout direction = 'z')
%--------------------------------------------------------------------------
siemens_twix_path1  = 'D:\data_pulseq_bstar\pulseq_Zhao_ball2_20221002\trajectory\meas_MID00384_FID46452_pulseq_Zhao_v2_cor_RL.dat';
ismrmrd_data_path1  = 'D:\data_pulseq_bstar\pulseq_Zhao_ball2_20221002\trajectory\meas_MID00384_FID46452_pulseq_Zhao_v2_cor_RL.h5';
ismrmrd_noise_path1 = 'D:\data_pulseq_bstar\pulseq_Zhao_ball2_20221002\trajectory\noise_meas_MID00384_FID46452_pulseq_Zhao_v2_cor_RL.h5';

%--------------------------------------------------------------------------
% Coronal, F >> H (readout direction = 'x')
%--------------------------------------------------------------------------
siemens_twix_path2  = 'D:\data_pulseq_bstar\pulseq_Zhao_ball2_20221002\trajectory\meas_MID00386_FID46454_pulseq_Zhao_v2_cor_FH.dat';
ismrmrd_data_path2  = 'D:\data_pulseq_bstar\pulseq_Zhao_ball2_20221002\trajectory\meas_MID00386_FID46454_pulseq_Zhao_v2_cor_FH.h5';
ismrmrd_noise_path2 = 'D:\data_pulseq_bstar\pulseq_Zhao_ball2_20221002\trajectory\noise_meas_MID00386_FID46454_pulseq_Zhao_v2_cor_FH.h5';

%--------------------------------------------------------------------------
% Transversal, A >> P (readout direction = 'x')
%--------------------------------------------------------------------------
siemens_twix_path3  = 'D:\data_pulseq_bstar\pulseq_Zhao_ball2_20221002\trajectory\meas_MID00382_FID46450_pulseq_Zhao_v2_tra_AP.dat';
ismrmrd_data_path3  = 'D:\data_pulseq_bstar\pulseq_Zhao_ball2_20221002\trajectory\meas_MID00382_FID46450_pulseq_Zhao_v2_tra_AP.h5';
ismrmrd_noise_path3 = 'D:\data_pulseq_bstar\pulseq_Zhao_ball2_20221002\trajectory\noise_meas_MID00382_FID46450_pulseq_Zhao_v2_tra_AP.h5';

%--------------------------------------------------------------------------
% Transversal, R >> L (readout direction = 'y')
%--------------------------------------------------------------------------
siemens_twix_path4  = 'D:\data_pulseq_bstar\pulseq_Zhao_ball2_20221002\trajectory\meas_MID00387_FID46455_pulseq_Zhao_v2_tra_RL.dat';
ismrmrd_data_path4  = 'D:\data_pulseq_bstar\pulseq_Zhao_ball2_20221002\trajectory\meas_MID00387_FID46455_pulseq_Zhao_v2_tra_RL.h5';
ismrmrd_noise_path4 = 'D:\data_pulseq_bstar\pulseq_Zhao_ball2_20221002\trajectory\noise_meas_MID00387_FID46455_pulseq_Zhao_v2_tra_RL.h5';

%--------------------------------------------------------------------------
% Sagittal, A >> P (readout direction = 'z')
%--------------------------------------------------------------------------
siemens_twix_path5  = 'D:\data_pulseq_bstar\pulseq_Zhao_ball2_20221002\trajectory\meas_MID00383_FID46451_pulseq_Zhao_v2_sag_AP.dat';
ismrmrd_data_path5  = 'D:\data_pulseq_bstar\pulseq_Zhao_ball2_20221002\trajectory\meas_MID00383_FID46451_pulseq_Zhao_v2_sag_AP.h5';
ismrmrd_noise_path5 = 'D:\data_pulseq_bstar\pulseq_Zhao_ball2_20221002\trajectory\noise_meas_MID00383_FID46451_pulseq_Zhao_v2_sag_AP.h5';

%--------------------------------------------------------------------------
% Sagittal, H >> F (readout direction = 'y')
%--------------------------------------------------------------------------
siemens_twix_path6  = 'D:\data_pulseq_bstar\pulseq_Zhao_ball2_20221002\trajectory\meas_MID00385_FID46453_pulseq_Zhao_v2_sag_HF.dat';
ismrmrd_data_path6  = 'D:\data_pulseq_bstar\pulseq_Zhao_ball2_20221002\trajectory\meas_MID00385_FID46453_pulseq_Zhao_v2_sag_HF.h5';
ismrmrd_noise_path6 = 'D:\data_pulseq_bstar\pulseq_Zhao_ball2_20221002\trajectory\noise_meas_MID00385_FID46453_pulseq_Zhao_v2_sag_HF.h5';

%% Make an output directory
output_path = fullfile(pwd, output_filename);
mkdir(output_path);

%% Process the first dataset (read_direction = 'z')
[trj_header1, trj1, trj_filename1, read_direction1, phase_direction1] = process_siemens_pulseq_Zhao_2020_MRM_bart_recon(siemens_twix_path1, ismrmrd_data_path1, ismrmrd_noise_path1, seq_scope_path, seq_bstar_path);

%% Process the second dataset (read_direction = 'x')
[trj_header2, trj2, trj_filename2, read_direction2, phase_direction2] = process_siemens_pulseq_Zhao_2020_MRM_bart_recon(siemens_twix_path2, ismrmrd_data_path2, ismrmrd_noise_path2, seq_scope_path, seq_bstar_path);

%% Process the third dataset (read_direction = 'x')
[trj_header3, trj3, trj_filename3, read_direction3, phase_direction3] = process_siemens_pulseq_Zhao_2020_MRM_bart_recon(siemens_twix_path3, ismrmrd_data_path3, ismrmrd_noise_path3, seq_scope_path, seq_bstar_path);

%% Process the fourth dataset (read_direction = 'y')
[trj_header4, trj4, trj_filename4, read_direction4, phase_direction4] = process_siemens_pulseq_Zhao_2020_MRM_bart_recon(siemens_twix_path4, ismrmrd_data_path4, ismrmrd_noise_path4, seq_scope_path, seq_bstar_path);

%% Process the fifth dataset (read_direction = 'z')
[trj_header5, trj5, trj_filename5, read_direction5, phase_direction5] = process_siemens_pulseq_Zhao_2020_MRM_bart_recon(siemens_twix_path5, ismrmrd_data_path5, ismrmrd_noise_path5, seq_scope_path, seq_bstar_path);

%% Process the sixth dataset (read_direction = 'y')
[trj_header6, trj6, trj_filename6, read_direction6, phase_direction6] = process_siemens_pulseq_Zhao_2020_MRM_bart_recon(siemens_twix_path6, ismrmrd_data_path6, ismrmrd_noise_path6, seq_scope_path, seq_bstar_path);

%% Create a .trj file
trj_header = trj_header1;
trj_header.lNumDirs = 6;
trj = (trj1 + trj2 + trj3 + trj4 + trj5 + trj6) / 2;
trj_path = fullfile(output_path, trj_filename1);
create_trj_file(trj_header, trj, trj_path);

%% Display a .trj file
display_trj_file;
