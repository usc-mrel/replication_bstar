% demo_figure4.m
% Written by Nam Gyun Lee
% Email: namgyunl@usc.edu, ggang56@gmail.com (preferred)
% Started: 07/27/2022, Last modified: 06/14/2023

%% Clean slate
close all; clear all; clc;

%% Set source directories
thirdparty_directory = 'D:\replication_bstar\thirdparty';

%% Add source directories to search path
addpath(genpath(thirdparty_directory));
start_time = tic;

%% Read a .cfl file
cfl_path1 = 'D:\data_pulseq_bstar\pulseq_bstar_acr_phantom_20230614\pulseq_bstar_1_38ms_1_61mm_i89_31k_WASP\img_pics_combined_kappa1_l0.0001_i30_360x360x360_res1.519_nlinv_a16_b16_i25';
cfl_path2 = 'D:\data_pulseq_bstar\pulseq_bstar_acr_phantom_20230614\pulseq_bstar_1_38ms_1_61mm_i88_31k_WASP\img_pics_combined_kappa1_l0.0001_i30_360x360x360_res1.519_nlinv_a16_b16_i25';
cfl_path3 = 'D:\data_pulseq_bstar\pulseq_bstar_acr_phantom_20230614\pulseq_bstar_1_38ms_1_61mm_i89_31k_SP\img_pics_combined_kappa1_l0.0001_i30_360x360x360_res1.519_nlinv_a16_b16_i25';
cfl_path4 = 'D:\data_pulseq_bstar\pulseq_bstar_acr_phantom_20230614\pulseq_bstar_1_38ms_1_61mm_i88_31k_SP\img_pics_combined_kappa1_l0.0001_i30_360x360x360_res1.519_nlinv_a16_b16_i25';

%--------------------------------------------------------------------------
% WASP 89 interleaves
%--------------------------------------------------------------------------
tstart = tic; fprintf('Reading a .cfl file: %s\n', cfl_path1);
img_wasp_89 = readcfl(cfl_path1);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%--------------------------------------------------------------------------
% WASP 88 interleaves
%--------------------------------------------------------------------------
tstart = tic; fprintf('Reading a .cfl file: %s\n', cfl_path2);
img_wasp_88 = readcfl(cfl_path2);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%--------------------------------------------------------------------------
% SP 89 interleaves
%--------------------------------------------------------------------------
tstart = tic; fprintf('Reading a .cfl file: %s\n', cfl_path3);
img_sp_89 = readcfl(cfl_path3);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%--------------------------------------------------------------------------
% SP 88 interleaves
%--------------------------------------------------------------------------
tstart = tic; fprintf('Reading a .cfl file: %s\n', cfl_path4);
img_sp_88 = readcfl(cfl_path4);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%% Crop images
[N1_osf,N2_osf,N3_osf] = size(img_wasp_89);

N1 = 256;
N2 = 256;
N3 = 256;

index1_range = (-floor(N1/2):ceil(N1/2)-1).' + floor(N1_osf/2) + 1;
index2_range = (-floor(N2/2):ceil(N2/2)-1).' + floor(N2_osf/2) + 1;
index3_range = (-floor(N3/2):ceil(N3/2)-1).' + floor(N3_osf/2) + 1;

im1 = img_wasp_89(index1_range,index2_range,index3_range);
im2 = img_wasp_88(index1_range,index2_range,index3_range);
im3 = img_sp_89(index1_range,index2_range,index3_range);
im4 = img_sp_88(index1_range,index2_range,index3_range);

N1_zoom = 220;
N2_zoom = 220;
N3_zoom = 220;

index1_range_zoom = (-floor(N1_zoom/2):ceil(N1_zoom/2)-1).' + floor(N1/2) + 1;
index2_range_zoom = (-floor(N2_zoom/2):ceil(N2_zoom/2)-1).' + floor(N2/2) + 1;
index3_range_zoom = (-floor(N3_zoom/2):ceil(N3_zoom/2)-1).' + floor(N3/2) + 1;

%% Create Figure 3
FontSize = 20;

coronal_slice_nr = 97;
sagittal_slice_nr = 132;

climits = [0 120];

figure('Color', 'w', 'Position', [-2 178 1191 800]);

%--------------------------------------------------------------------------
% Coronal (WASP 89)
%--------------------------------------------------------------------------
ax1 = subplot(2,4,1);
imagesc(abs(squeeze(im1(index1_range_zoom,coronal_slice_nr,index3_range_zoom+5))));
axis image off;
colormap(gray(256));
clim(climits);
text(N2_zoom, -38, {'WASP'}, 'FontSize', FontSize, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
text(N2_zoom / 2, 0, {'89 interleaves'}, 'FontSize', FontSize, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
text(5, N2_zoom - 4, {'(A)'}, 'Rotation', 0, 'Color', 'w', 'FontSize', FontSize, 'FontWeight', 'Bold', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');

%--------------------------------------------------------------------------
% Coronal (WASP 88)
%--------------------------------------------------------------------------
ax2 = subplot(2,4,2);
imagesc(abs(squeeze(im2(index1_range_zoom,coronal_slice_nr,index3_range_zoom+5))));
axis image off;
colormap(gray(256));
clim(climits);
text(N2_zoom / 2, 0, {'88 interleaves'}, 'FontSize', FontSize, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
text(5, N2_zoom - 4, {'(C)'}, 'Rotation', 0, 'Color', 'w', 'FontSize', FontSize, 'FontWeight', 'Bold', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');

%--------------------------------------------------------------------------
% Coronal (SP 89)
%--------------------------------------------------------------------------
ax3 = subplot(2,4,3);
imagesc(abs(squeeze(im3(index1_range_zoom,coronal_slice_nr,index3_range_zoom+5))));
axis image off;
colormap(gray(256));
clim(climits);
text(N2_zoom, -38, {'Spiral phyllotaxis'}, 'FontSize', FontSize, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
text(N2_zoom / 2, 0, {'89 interleaves'}, 'FontSize', FontSize, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
text(5, N2_zoom - 4, {'(E)'}, 'Rotation', 0, 'Color', 'w', 'FontSize', FontSize, 'FontWeight', 'Bold', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');

%--------------------------------------------------------------------------
% Coronal (SP 88)
%--------------------------------------------------------------------------
ax4 = subplot(2,4,4);
imagesc(abs(squeeze(im4(index1_range_zoom,coronal_slice_nr,index3_range_zoom+5))));
axis image off;
colormap(gray(256));
clim(climits);
text(N2_zoom / 2, 0, {'88 interleaves'}, 'FontSize', FontSize, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
text(5, N2_zoom - 4, {'(G)'}, 'Rotation', 0, 'Color', 'w', 'FontSize', FontSize, 'FontWeight', 'Bold', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');

%--------------------------------------------------------------------------
% Sagittal (WASP 89)
%--------------------------------------------------------------------------
ax5 = subplot(2,4,5);
imagesc(abs(squeeze(im1(index1_range_zoom,index2_range_zoom,sagittal_slice_nr))).');
axis image off;
colormap(gray(256));
clim(climits);
text(5, N2_zoom - 4, {'(B)'}, 'Rotation', 0, 'Color', 'w', 'FontSize', FontSize, 'FontWeight', 'Bold', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');

%--------------------------------------------------------------------------
% Sagittal (WASP 88)
%--------------------------------------------------------------------------
ax6 = subplot(2,4,6);
imagesc(abs(squeeze(im2(index1_range_zoom,index2_range_zoom,sagittal_slice_nr))).');
axis image off;
colormap(gray(256));
clim(climits);
text(5, N2_zoom - 4, {'(D)'}, 'Rotation', 0, 'Color', 'w', 'FontSize', FontSize, 'FontWeight', 'Bold', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');

%--------------------------------------------------------------------------
% Sagittal (SP 89)
%--------------------------------------------------------------------------
ax7 = subplot(2,4,7);
imagesc(abs(squeeze(im3(index1_range_zoom,index2_range_zoom,sagittal_slice_nr))).');
axis image off;
colormap(gray(256));
clim(climits);
text(5, N2_zoom - 4, {'(F)'}, 'Rotation', 0, 'Color', 'w', 'FontSize', FontSize, 'FontWeight', 'Bold', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');

%--------------------------------------------------------------------------
% Sagittal (SP 88)
%--------------------------------------------------------------------------
ax8 = subplot(2,4,8);
imagesc(abs(squeeze(im4(index1_range_zoom,index2_range_zoom,sagittal_slice_nr))).');
axis image off;
colormap(gray(256));
clim(climits);
text(5, N2_zoom - 4, {'(H)'}, 'Rotation', 0, 'Color', 'w', 'FontSize', FontSize, 'FontWeight', 'Bold', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');

set(ax1, 'Position', [0.1300-0.005*2+0.001*-2           0.5838-0.16+0.002*1 0.1566+0.05 0.3412+0.05]);
set(ax2, 'Position', [0.3361-0.005*1+0.001*-1           0.5838-0.16+0.002*1 0.1566+0.05 0.3412+0.05]);
set(ax3, 'Position', [0.5422-0.005*0+0.001*+0+0.0067*+1 0.5838-0.16+0.002*1 0.1566+0.05 0.3412+0.05]);
set(ax4, 'Position', [0.7484+0.005*1+0.001*+1+0.0067*+1 0.5838-0.16+0.002*1 0.1566+0.05 0.3412+0.05]);
set(ax5, 'Position', [0.1300-0.005*2+0.001*-2           0.1100              0.1566+0.05 0.3412+0.05]);
set(ax6, 'Position', [0.3361-0.005*1+0.001*-1           0.1100              0.1566+0.05 0.3412+0.05]);
set(ax7, 'Position', [0.5422-0.005*0+0.001*+0+0.0067*+1 0.1100              0.1566+0.05 0.3412+0.05]);
set(ax8, 'Position', [0.7484+0.005*1+0.001*+1+0.0067*+1 0.1100              0.1566+0.05 0.3412+0.05]);

% Create line
annotation(gcf, 'line', [0.5428 0.5428], [0.82025 0.15125]); % 0.82025 = height

export_fig('figure4', '-r500', '-tif');
