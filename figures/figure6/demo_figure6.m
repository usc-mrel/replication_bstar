% demo_figure6.m
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

%% Read in a .dcm file for open-source bSTAR
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

%% Read in a .dcm file for original bSTAR
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
FontSize = 20;

[N1,N2,N3] = size(img_pulseq);
c1 = floor(N1/2) + 1;
c2 = floor(N2/2) + 1;
c3 = floor(N2/2) + 1;

scale_factor1 = abs(img_pulseq(161,c2,c3));
scale_factor2 = abs(img_idea(161,c2,c3));

index1_range = (1:N1).';
index2_range = (1:N2).';
index3_range = (1:N3).';

img1_scaled = img_pulseq / scale_factor1;
img2_scaled = img_idea   / scale_factor2;

coronal_slice_nr1 = 256 - 116 + 1;
coronal_slice_nr2 = 256 - 116 + 1;

axial_slice_nr1 = 256 - 108 + 1;
axial_slice_nr2 = 256 - 116 + 1;

figure('Color', 'w', 'Position', [10 2 927 990]);
color_order = get(gca, 'ColorOrder');

climits1 = [0 2.5];
climits2 = [0 2.5];

%--------------------------------------------------------------------------
% Open-source bSTAR (coronal)
%--------------------------------------------------------------------------
ax1 = subplot(2,2,1);
imagesc(flip(rot90(squeeze(img1_scaled(coronal_slice_nr1,index2_range,index3_range)),-1),2));
axis image off;
colormap(gray(256));
clim(climits1);
text(N2 / 2, 0, {'Open-source bSTAR', 'Volunteer 1 (day 1)'}, 'FontSize', FontSize, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
text(5, N2 - 4, {'(A)'}, 'Rotation', 0, 'Color', 'w', 'FontSize', FontSize, 'FontWeight', 'Bold', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');

%--------------------------------------------------------------------------
% Original bSTAR (coronal)
%--------------------------------------------------------------------------
ax2 = subplot(2,2,2);
imagesc(flip(rot90(squeeze(img2_scaled(coronal_slice_nr2,index2_range,index3_range)),-1),2));
axis image off;
colormap(gray(256));
clim(climits2);
text(N2 / 2, 0, {'Original bSTAR', 'Volunteer 1 (day 2)'}, 'FontSize', FontSize, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
text(5, N2 - 4, {'(C)'}, 'Rotation', 0, 'Color', 'w', 'FontSize', FontSize, 'FontWeight', 'Bold', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');

%--------------------------------------------------------------------------
% Open-source bSTAR (axial)
%--------------------------------------------------------------------------
ax3 = subplot(2,2,3);
imagesc(squeeze(img1_scaled(index1_range,index2_range,axial_slice_nr1)));
axis image off;
colormap(gray(256));
clim(climits1);
text(5, N2 - 4, {'(B)'}, 'Rotation', 0, 'Color', 'w', 'FontSize', FontSize, 'FontWeight', 'Bold', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');

%--------------------------------------------------------------------------
% Original bSTAR (axial)
%--------------------------------------------------------------------------
ax4 = subplot(2,2,4);
imagesc(squeeze(img2_scaled(index1_range,index2_range,axial_slice_nr2)));
axis image off;
colormap(gray(256));
clim(climits2);
text(5, N2 - 4, {'(D)'}, 'Rotation', 0, 'Color', 'w', 'FontSize', FontSize, 'FontWeight', 'Bold', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');

set(ax1, 'Position', [0.1300-0.1 0.1100+0.2 0.3347+0.1 0.8150]);
set(ax2, 'Position', [0.5703-0.1 0.1100+0.2 0.3347+0.1 0.8150]);
set(ax3, 'Position', [0.1300-0.1 0.1100-0.2118 0.3347+0.1 0.8150]);
set(ax4, 'Position', [0.5703-0.1 0.1100-0.2118 0.3347+0.1 0.8150]);

export_fig('figure6', '-r700', '-tif');
