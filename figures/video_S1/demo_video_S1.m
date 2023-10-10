% demo_video_S1.m
% Written by Nam Gyun Lee
% Email: namgyunl@usc.edu, ggang56@gmail.com (preferred)
% Started: 11/07/2022, Last modified: 06/29/2023

%% Clean slate
close all; clear all; clc;

%% Set source directories
thirdparty_directory = 'D:\replication_bstar\thirdparty';

%% Add source directories to search path
addpath(genpath(thirdparty_directory));
start_time = tic;

%% Read a .cfl file
cfl_path1 = 'D:\data_pulseq_bstar\pulseq_bstar_vol422_20220830\pulseq_bstar_1_38ms_1_60mm_bh_40k_expiratory\img_pics_echo1_kappa1_l0.0001_i30_360x360x360_res1.519_nlinv_a16_b16_i25';
cfl_path2 = 'D:\data_pulseq_bstar\pulseq_bstar_vol422_20220830\pulseq_bstar_1_38ms_1_60mm_bh_40k_expiratory\img_pics_echo2_kappa1_l0.0001_i30_360x360x360_res1.519_nlinv_a16_b16_i25';
cfl_path3 = 'D:\data_pulseq_bstar\pulseq_bstar_vol422_20220830\pulseq_bstar_1_38ms_1_60mm_bh_40k_expiratory\img_pics_combined_kappa1_l0.0001_i30_360x360x360_res1.519_nlinv_a16_b16_i25';

%--------------------------------------------------------------------------
% Echo 1
%--------------------------------------------------------------------------
tstart = tic; fprintf('Reading a .cfl file: %s\n', cfl_path1);
img1 = readcfl(cfl_path1);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%--------------------------------------------------------------------------
% Echo 2
%--------------------------------------------------------------------------
tstart = tic; fprintf('Reading a .cfl file: %s\n', cfl_path2);
img2 = readcfl(cfl_path2);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%--------------------------------------------------------------------------
% Echo combined
%--------------------------------------------------------------------------
tstart = tic; fprintf('Reading a .cfl file: %s\n', cfl_path3);
img3 = readcfl(cfl_path3);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%% Calculate the range of indices
[N1,N2,N3] = size(img1);

c1 = floor(N1/2) + 1;
c2 = floor(N2/2) + 1;
c3 = floor(N3/2) + 1;

scale_factor1 = abs(img1(c1,c2,c3));
scale_factor2 = abs(img2(c1,c2,c3));
scale_factor3 = abs(img3(c1,c2,c3));

img1_scaled = img1 / scale_factor1;
img2_scaled = img2 / scale_factor2;
img3_scaled = img3 / scale_factor3;

%% Calculate the range of indices
[N1_osf,N2_osf,N3_osf] = size(img1);

N1_zoom = 200;
N2_zoom = 140;
N3_zoom = 212;

index1_range_zoom = (-floor(N1_zoom/2):ceil(N1_zoom/2)-1).' + floor(N1/2) + 1;
index2_range_zoom = (-floor(N2_zoom/2):ceil(N2_zoom/2)-1).' + floor(N2/2) + 1;
index3_range_zoom = (-floor(N3_zoom/2):ceil(N3_zoom/2)-1).' + floor(N3/2) + 1;

%% Display Figure 4
FontSize = 20;
climits = [0 2.5];

figure('Color', 'w', 'Position', [-2 178 1918 800]);
color_order = get(gca, 'colororder');

video_fullpath = fullfile(pwd, sprintf('video_S1'));
video_file = VideoWriter(video_fullpath, 'MPEG-4');
video_file.Quality = 100;

video_file.FrameRate = 3;
open(video_file);

for idx = 40:220
    %----------------------------------------------------------------------
    % Sagittal view
    %----------------------------------------------------------------------
    sagittal_slice_nr = idx + 40;
    vertical_block = 100 * ones(N1_zoom, 1, 'double');
    sagittal_montage = cat(2, img1_scaled(index1_range_zoom + 30, index2_range_zoom + 30, sagittal_slice_nr), vertical_block, ...
                              img2_scaled(index1_range_zoom + 30, index2_range_zoom + 30, sagittal_slice_nr), vertical_block, ...
                              img3_scaled(index1_range_zoom + 30, index2_range_zoom + 30, sagittal_slice_nr));

    %----------------------------------------------------------------------
    % Axial view
    %----------------------------------------------------------------------
    axial_slice_nr = idx + 80;
    horizontal_block = 100 * ones(1, 350, 'double');
    axial_idx2_range = (-floor(140/2):ceil(140/2)-1).' + floor(360/2) + 1;
    axial_idx3_range = (-floor(350/2):ceil(350/2)-1).' + floor(360/2) + 1;
    axial_montage = cat(1, squeeze(img1_scaled(axial_slice_nr, axial_idx2_range + 25, axial_idx3_range + 5)), horizontal_block, ...
                           squeeze(img2_scaled(axial_slice_nr, axial_idx2_range + 25, axial_idx3_range + 5)), horizontal_block, ...
                           squeeze(img3_scaled(axial_slice_nr, axial_idx2_range + 25, axial_idx3_range + 5)));

    %----------------------------------------------------------------------
    % Coronal view
    %----------------------------------------------------------------------
    coronal_slice_nr = idx + 100;
    coronal_idx1_range = (-floor(240/2):ceil(240/2)-1).' + floor(360/2) + 1;
    coronal_idx3_range = (-floor(280/2):ceil(280/2)-1).' + floor(360/2) + 1;
    vertical_block = 100 * ones(length(coronal_idx1_range), 1, 'double');
    coronal_montage = cat(2, squeeze(img1_scaled(coronal_idx1_range + 20, coronal_slice_nr, coronal_idx3_range + 15)), vertical_block, ...
                             squeeze(img2_scaled(coronal_idx1_range + 20, coronal_slice_nr, coronal_idx3_range + 15)), vertical_block, ...
                             squeeze(img3_scaled(coronal_idx1_range + 20, coronal_slice_nr, coronal_idx3_range + 15)));

    ax1 = subplot(2,2,1);
    imagesc(abs(sagittal_montage));
    axis image off;
    colormap(gray(256));
    clim(climits);
    text(140 + 140*0, 0, 'Echo 1'  , 'FontSize', FontSize-2, 'Color', color_order(3,:), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');
    text(140 + 140*1, 0, 'Echo 2'  , 'FontSize', FontSize-2, 'Color', color_order(3,:), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');
    text(140 + 140*2, 0, 'Combined', 'FontSize', FontSize-2, 'Color', color_order(3,:), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');
    text(3 + 140*0    , 0, {'(A)'}, 'Rotation', 0, 'Color', 'w', 'FontSize', FontSize, 'FontWeight', 'Bold', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');
    text(3 + 140*1 + 1, 0, {'(B)'}, 'Rotation', 0, 'Color', 'w', 'FontSize', FontSize, 'FontWeight', 'Bold', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');
    text(3 + 140*2 + 2, 0, {'(C)'}, 'Rotation', 0, 'Color', 'w', 'FontSize', FontSize, 'FontWeight', 'Bold', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');

    ax2 = subplot(2,2,[2,4]);
    imagesc(abs(axial_montage));
    axis image off;
    colormap(gray(256));
    clim(climits);
    text(3, 140*0    , {'(G)'}, 'Rotation', 0, 'Color', 'w', 'FontSize', FontSize, 'FontWeight', 'Bold', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');
    text(3, 140*1 + 1, {'(H)'}, 'Rotation', 0, 'Color', 'w', 'FontSize', FontSize, 'FontWeight', 'Bold', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');
    text(3, 140*2 + 2, {'(I)'}, 'Rotation', 0, 'Color', 'w', 'FontSize', FontSize, 'FontWeight', 'Bold', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');
    text(350-2, 140*0 + 0, 'Echo 1'  , 'FontSize', FontSize-2, 'Color', color_order(3,:), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');
    text(350-2, 140*1 + 1, 'Echo 2'  , 'FontSize', FontSize-2, 'Color', color_order(3,:), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');
    text(350-2, 140*2 + 2, 'Combined', 'FontSize', FontSize-2, 'Color', color_order(3,:), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');

    ax3 = subplot(2,2,3);
    imagesc(abs(coronal_montage));
    axis image off;
    colormap(gray(256));
    clim(climits);
    text(5 + 280*0    , 0, {'(D)'}, 'Rotation', 0, 'Color', 'w', 'FontSize', FontSize, 'FontWeight', 'Bold', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');
    text(5 + 280*1 + 1, 0, {'(E)'}, 'Rotation', 0, 'Color', 'w', 'FontSize', FontSize, 'FontWeight', 'Bold', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');
    text(5 + 280*2 + 2, 0, {'(F)'}, 'Rotation', 0, 'Color', 'w', 'FontSize', FontSize, 'FontWeight', 'Bold', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');
    text(280 + 280*0, 0, 'Echo 1'  , 'FontSize', FontSize-2, 'Color', color_order(3,:), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');
    text(280 + 280*1, 0, 'Echo 2'  , 'FontSize', FontSize-2, 'Color', color_order(3,:), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');
    text(280 + 280*2, 0, 'Combined', 'FontSize', FontSize-2, 'Color', color_order(3,:), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');

    set(ax1, 'Position', [0.1300-0.005 0.5838-0.16 0.3347+0.16 0.3412+0.16]);
    set(ax3, 'Position', [0.1300+0.0215 0.1100-0.073 0.3347+0.107 0.3412+0.107]);
    set(ax2, 'Position', [0.5703 0.1100 0.3347 0.8150]);

    frame = getframe(gcf);
    frame.cdata = frame.cdata(59:714,290:1688,:);
    writeVideo(video_file, frame);
end
close(video_file);
