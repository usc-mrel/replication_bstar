% demo_video_S2.m
% Written by Nam Gyun Lee
% Email: namgyunl@usc.edu, ggang56@gmail.com (preferred)
% Started: 07/27/2022, Last modified: 08/20/2022

%% Clean slate
close all; clear all; clc;

%% Set source directories
thirdparty_directory = 'D:\replication_bstar\thirdparty';

%% Add source directories to search path
addpath(genpath(thirdparty_directory));

%% Define input directories
pulseq_bstar_dicom_directory = 'D:\data_pulseq_bstar\pulseq_bstar_vol490_20221104\pulseq_bstar_1_38ms_1_60mm_bh_17k_expiratory\img_pics_combined_kappa1_l0.0001_i30_360x360x360_res1.519_nlinv_a16_b16_i25_dicom256';
idea_bstar_dicom_directory = 'D:\data_bstar\bSTAR_USC_vol411_20220810\dicom\DCM_BH_1.6mm_1.38ms';

%% Read in a .dcm file for Pulseq bSTAR
dir_info = dir(fullfile(pulseq_bstar_dicom_directory, '*dcm'));
nr_dicom_files = length(dir_info);

%--------------------------------------------------------------------------
% Get DICOM information from the first dcm file
%--------------------------------------------------------------------------
dicom_file = fullfile(dir_info(1).folder, dir_info(1).name);
dicom_info = dicominfo(dicom_file);

N1 = dicom_info.Rows;
N2 = dicom_info.Columns;

%--------------------------------------------------------------------------
% Read all dcm files
%--------------------------------------------------------------------------
start_time = tic;
img_pulseq = zeros(N1, N2, nr_dicom_files, 'double');
for idx = 1:nr_dicom_files
    dicom_file = fullfile(dir_info(idx).folder, dir_info(idx).name);
    tstart = tic; fprintf('Reading a dicom file: %s... ', dicom_file);
    dicom_info = dicominfo(dicom_file);
    img_pulseq(:,:,idx) = double(dicomread(dicom_info));
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
end

%% Read in a .dcm file for IDEA bSTAR
dir_info = dir(fullfile(idea_bstar_dicom_directory, '*dcm'));
nr_dicom_files = length(dir_info);

%--------------------------------------------------------------------------
% Read all dcm files
%--------------------------------------------------------------------------
img_idea = zeros(N1, N2, nr_dicom_files, 'double');
for idx = 1:nr_dicom_files
    dicom_file = fullfile(dir_info(idx).folder, dir_info(idx).name);
    tstart = tic; fprintf('Reading a dicom file: %s... ', dicom_file);
    dicom_info = dicominfo(dicom_file);
    img_idea(:,:,idx) = double(dicomread(dicom_info));
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
end

%% Display images
FontSize = 25;

[N1,N2,N3] = size(img_pulseq);
c1 = floor(N1/2) + 1;
c2 = floor(N2/2) + 1;

scale_factor1 = abs(img_pulseq(161,c1,c2));
scale_factor2 = abs(img_idea(161,c1,c2));

index2_range = (1:N2).';
index3_range = (1:N3).';

figure('Color', 'k', 'Position', [330 11 1575 835]);
color_order = get(gca, 'ColorOrder');

video_fullpath = fullfile(pwd, sprintf('video_S2'));
%video_file = VideoWriter(video_fullpath, 'Uncompressed AVI');
video_file = VideoWriter(video_fullpath, 'MPEG-4');
%video_file = VideoWriter(video_fullpath, 'Motion JPEG AVI');
video_file.Quality = 100;

video_file.FrameRate = 3;
open(video_file);

for idx = 100:220
    % Coronal
    img1 = flip(rot90(squeeze(img_pulseq(idx,index2_range,index3_range)),-1),2);
    img2 = flip(rot90(squeeze(img_idea(idx,index2_range,index3_range)),-1),2);

    img1_scaled = img1 / scale_factor1;
    img2_scaled = img2 / scale_factor2;

    ax1 = subplot(1,2,1);
    imagesc(img1_scaled);
    axis image off;
    colormap(gray(256));
    clim([0 2.5]);
    text(N2/2, 0, {'Open-source bSTAR', 'Volunteer 1 (day 1)'}, 'Color', 'w', 'FontSize', FontSize, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

    ax2 = subplot(1,2,2);
    imagesc(img2_scaled);
    axis image off;
    colormap(gray(256));
    clim([0 2.5]);
    text(N2/2, 0, {'Original bSTAR', 'Volunteer 1 (day 2)'}, 'Color', 'w', 'FontSize', FontSize, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

    set(ax1, 'Position', [0.1300-0.1 0.1100-0.04 0.3347+0.15 0.8150]);
    set(ax2, 'Position', [0.5703-0.1 0.1100-0.04 0.3347+0.15 0.8150]);

    writeVideo(video_file, getframe(gcf));
end
close(video_file);