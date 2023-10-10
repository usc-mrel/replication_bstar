% demo_figure3.m
% Written by Nam Gyun Lee
% Email: namgyunl@usc.edu, ggang56@gmail.com (preferred)
% Started: 07/27/2022, Last modified: 06/13/2023

%% Clean slate
close all; clear all; clc;

%% Set source directories
thirdparty_directory = 'D:\replication_bstar\thirdparty';

%% Add source directories to search path
addpath(genpath(thirdparty_directory));
start_time = tic;

%% Read a .cfl file
cfl_path1 = 'D:\data_pulseq_bstar\pulseq_bstar_acr_phantom_20230614\pulseq_bstar_1_38ms_1_61mm_i89_31k_SP\img_pics_echo1_kappa1_l0.0001_i30_360x360x360_res1.519_nlinv_a16_b16_i25';
cfl_path2 = 'D:\data_pulseq_bstar\pulseq_bstar_acr_phantom_20230614\pulseq_bstar_1_38ms_1_61mm_i89_31k_SP\img_pics_echo2_kappa1_l0.0001_i30_360x360x360_res1.519_nlinv_a16_b16_i25';

%--------------------------------------------------------------------------
% SP 89 interleaves (echo 1)
%--------------------------------------------------------------------------
tstart = tic; fprintf('Reading a .cfl file: %s\n', cfl_path1);
img_pics_echo1 = readcfl(cfl_path1);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%--------------------------------------------------------------------------
% SP 89 interleaves (echo 2)
%--------------------------------------------------------------------------
tstart = tic; fprintf('Reading a .cfl file: %s\n', cfl_path2);
img_pics_echo2 = readcfl(cfl_path2);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%% Crop images
[N1_osf,N2_osf,N3_osf] = size(img_pics_echo1);

N1 = 256;
N2 = 256;
N3 = 256;

index1_range = (-floor(N1/2):ceil(N1/2)-1).' + floor(N1_osf/2) + 1;
index2_range = (-floor(N2/2):ceil(N2/2)-1).' + floor(N2_osf/2) + 1;
index3_range = (-floor(N3/2):ceil(N3/2)-1).' + floor(N3_osf/2) + 1;

im1 = img_pics_echo1(index1_range,index2_range,index3_range);
im2 = img_pics_echo2(index1_range,index2_range,index3_range);

N1_zoom = 220;
N2_zoom = 220;
N3_zoom = 220;

index1_range_zoom = (-floor(N1_zoom/2):ceil(N1_zoom/2)-1).' + floor(N1/2) + 1;
index2_range_zoom = (-floor(N2_zoom/2):ceil(N2_zoom/2)-1).' + floor(N2/2) + 1;
index3_range_zoom = (-floor(N3_zoom/2):ceil(N3_zoom/2)-1).' + floor(N3/2) + 1;

%% Create Figure 2
FontSize = 20;

coronal_slice_nr = 95;
sagittal_slice_nr = 128;
axial_slice_nr = 152;

climits = [0 60];

figure('Color', 'w', 'Position', [-2 178 1191 800]);

%--------------------------------------------------------------------------
% Coronal (echo 1)
%--------------------------------------------------------------------------
ax1 = subplot(2,3,1);
imagesc(abs(squeeze(im1(index1_range_zoom,coronal_slice_nr,index3_range_zoom))));
axis image off;
colormap(gray(256));
clim(climits);
text(0, N1_zoom / 2, {'Echo 1'}, 'Rotation', 90, 'FontSize', FontSize, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
text(5, N2_zoom - 4, {'(A)'}, 'Rotation', 0, 'Color', 'w', 'FontSize', FontSize, 'FontWeight', 'Bold', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');

%--------------------------------------------------------------------------
% Sagittal (echo 1)
%--------------------------------------------------------------------------
ax2 = subplot(2,3,2);
imagesc(abs(squeeze(im1(index1_range_zoom,index2_range_zoom,sagittal_slice_nr))).');
axis image off;
colormap(gray(256));
clim(climits);
text(5, N2_zoom - 4, {'(B)'}, 'Rotation', 0, 'Color', 'w', 'FontSize', FontSize, 'FontWeight', 'Bold', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');

%--------------------------------------------------------------------------
% Axial (echo 1)
%--------------------------------------------------------------------------
ax3 = subplot(2,3,3);
imagesc(abs(squeeze(im1(axial_slice_nr,index2_range_zoom,index3_range_zoom))));
axis image off;
colormap(gray(256));
clim(climits);
text(5, N2_zoom - 4, {'(C)'}, 'Rotation', 0, 'Color', 'w', 'FontSize', FontSize, 'FontWeight', 'Bold', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');

%--------------------------------------------------------------------------
% Coronal (echo 2)
%--------------------------------------------------------------------------
ax4 = subplot(2,3,4);
imagesc(abs(squeeze(im2(index1_range_zoom,coronal_slice_nr,index3_range_zoom))));
axis image off;
colormap(gray(256));
clim(climits);
text(0, N1_zoom / 2, {'Echo 2'}, 'Rotation', 90, 'FontSize', FontSize, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
text(5, N2_zoom - 4, {'(D)'}, 'Rotation', 0, 'Color', 'w', 'FontSize', FontSize, 'FontWeight', 'Bold', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');

%--------------------------------------------------------------------------
% Sagittal (echo 2)
%--------------------------------------------------------------------------
ax5 = subplot(2,3,5);
imagesc(abs(squeeze(im2(index1_range_zoom,index2_range_zoom,sagittal_slice_nr))).');
axis image off;
clim(climits);
text(5, N2_zoom - 4, {'(E)'}, 'Rotation', 0, 'Color', 'w', 'FontSize', FontSize, 'FontWeight', 'Bold', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');

%--------------------------------------------------------------------------
% Axial (echo 2)
%--------------------------------------------------------------------------
ax6 = subplot(2,3,6);
imagesc(abs(squeeze(im2(axial_slice_nr,index2_range_zoom,index3_range_zoom))));
axis image off;
clim(climits);
text(5, N2_zoom - 4, {'(F)'}, 'Rotation', 0, 'Color', 'w', 'FontSize', FontSize, 'FontWeight', 'Bold', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');

set(ax1, 'Position', [0.1300-0.03*2+0.004*-1 0.5838-0.078+0.0055*1 0.2134+0.05 0.3412+0.05]);
set(ax2, 'Position', [0.3361-0.03*0+0.004*0  0.5838-0.078+0.0055*1 0.2134+0.05 0.3412+0.05]);
set(ax3, 'Position', [0.5422+0.03*2+0.004*1  0.5838-0.078+0.0055*1 0.2134+0.05 0.3412+0.05]);
set(ax4, 'Position', [0.1300-0.03*2+0.004*-1 0.1100                0.2134+0.05 0.3412+0.05]);
set(ax5, 'Position', [0.3361-0.03*0+0.004*0  0.1100                0.2134+0.05 0.3412+0.05]);
set(ax6, 'Position', [0.5422+0.03*2+0.004*1  0.1100                0.2134+0.05 0.3412+0.05]);


% Create arrow
annotation(gcf, 'arrow', [0.245172124265323 0.225020990764064], [0.18125 0.2325], ...
    'Color', [251 176 64] / 255, 'LineWidth', 4, 'HeadWidth', 12, 'HeadStyle', 'plain', 'HeadLength', 12);

% Create arrow
annotation(gcf, 'arrow', [0.377833753148615 0.403862300587739], [0.45125 0.41375], ...
    'Color', [251 176 64] / 255, 'LineWidth', 4, 'HeadWidth', 12, 'HeadStyle', 'plain', 'HeadLength', 12);

% Create arrow
annotation(gcf, 'arrow', [0.55919395465995 0.533165407220822], [0.45125 0.41375],...
    'Color', [251 176 64] / 255, 'LineWidth', 4, 'HeadWidth', 12, 'HeadStyle', 'plain', 'HeadLength', 12);

% Create arrow
annotation(gcf, 'arrow', [0.638958858102432 0.664987405541557], [0.45125 0.41375],...
    'Color', [251 176 64] / 255, 'LineWidth', 4, 'HeadWidth', 12, 'HeadStyle', 'plain', 'HeadLength', 12);

% Create arrow
annotation(gcf, 'arrow', [0.84466834592779 0.818639798488662], [0.45125 0.41375],...
    'Color', [251 176 64] / 255, 'LineWidth', 4, 'HeadWidth', 12, 'HeadStyle', 'plain', 'HeadLength', 12);

export_fig('figure3', '-r500', '-tif');
