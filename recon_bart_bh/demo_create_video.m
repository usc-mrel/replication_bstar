% demo_create_video.m
% Written by Nam Gyun Lee
% Email: namgyunl@usc.edu, ggang56@gmail.com (preferred)
% Started: 04/27/2024, Last modified: 04/27/2024

%% Clean slate
close all; clear all; clc;

%% Define the full path of a .json file
json_file = 'D:\replication_bstar_v2\recon_bart_bh\meas_MID00189_FID13599_bstar_340mm_1_60mm_true0_TR1_54ms_rf200_i89_31k.json';

%% Define an output filename
img1_filename = 'img_rss_eco1_360x360x360';
img2_filename = 'img_rss_eco2_360x360x360';

img1_filename = 'img_pics_eco1_kappa1_l0.0001_i30_360x360x360';
img2_filename = 'img_pics_eco2_kappa1_l0.0001_i30_360x360x360';

%% Start a stopwatch timer
start_time = tic;

%% Read a .json file
tstart = tic; fprintf('%s: Reading a .json file: %s... ', datetime, json_file);
fid = fopen(json_file);
json_txt = fread(fid, [1 inf], 'char=>char');
fclose(fid);
json = jsondecode(json_txt);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%% Get parameters from a .json file
%--------------------------------------------------------------------------
% Define the full path of a filename
%--------------------------------------------------------------------------
if ispc
    output_path = strrep(json.output_path, '/', '\');
else
    output_path = json.output_path;
end

%--------------------------------------------------------------------------
% Reconstruction parameters
%--------------------------------------------------------------------------
recon_matrix_size = json.recon_conf.recon_matrix_size.';   % reconstruction matrix size

N1 = recon_matrix_size(1);
N2 = recon_matrix_size(2);
N3 = recon_matrix_size(3);

%% Read a .cfl file
%--------------------------------------------------------------------------
% x
%--------------------------------------------------------------------------
cfl_file = fullfile(output_path, 'x');
x = real(readcfl(cfl_file));

%--------------------------------------------------------------------------
% y
%--------------------------------------------------------------------------
cfl_file = fullfile(output_path, 'y');
y = real(readcfl(cfl_file));

%--------------------------------------------------------------------------
% z
%--------------------------------------------------------------------------
cfl_file = fullfile(output_path, 'z');
z = real(readcfl(cfl_file));

%--------------------------------------------------------------------------
% img1
%--------------------------------------------------------------------------
cfl_file = fullfile(output_path, img1_filename);
img1 = readcfl(cfl_file);

%--------------------------------------------------------------------------
% img2
%--------------------------------------------------------------------------
cfl_file = fullfile(output_path, img2_filename);
img2 = readcfl(cfl_file);

%% Display reconstructed images in the DCS
X = reshape(x, [N1 N2 N3]);
Y = reshape(y, [N1 N2 N3]);
Z = reshape(z, [N1 N2 N3]);

c1 = floor(N1/2) + 1;
c2 = floor(N2/2) + 1;
c3 = floor(N3/2) + 1;

FontSize = 14;

figure('Color', 'w', 'Position', [36 27 1466 951]);
ax1 = subplot(2,3,1);
ax2 = subplot(2,3,2);
ax3 = subplot(2,3,3);
ax4 = subplot(2,3,4);
ax5 = subplot(2,3,5);
ax6 = subplot(2,3,6);

video_file1 = fullfile(output_path, img1_filename);
video_file2 = fullfile(output_path, sprintf('%s_zoom', img1_filename));

video_obj1 = VideoWriter(video_file1, 'MPEG-4');
video_obj2 = VideoWriter(video_file2, 'MPEG-4');

video_obj1.Quality = 100;
video_obj2.Quality = 100;
video_obj1.FrameRate = 3;
video_obj2.FrameRate = 3;

open(video_obj1);
open(video_obj2);

for idx = 1:N3
    %----------------------------------------------------------------------
    % Echo 1, sagittal slice
    %----------------------------------------------------------------------
    surf(ax1, X(:,:,idx) * 1e3, Y(:,:,idx) * 1e3, Z(:,:,idx) * 1e3, abs(img1(:,:,idx)), 'EdgeColor', 'none');
    set(ax1, 'Box', 'on', 'TickLabelInterpreter', 'latex', 'FontSize', FontSize);
    axis(ax1, 'image');
    xlabel(ax1, 'x [mm]', 'Interpreter', 'latex', 'FontSize', FontSize);
    ylabel(ax1, 'y [mm]', 'Interpreter', 'latex', 'FontSize', FontSize);
    zlabel(ax1, 'z [mm]', 'Interpreter', 'latex', 'FontSize', FontSize);
    title(ax1, {sprintf('Echo 1, sagittal slice %d', idx), sprintf('x = %3.1f mm', X(1,1,idx) * 1e3)}, 'FontWeight', 'normal');
    set(ax1, 'YDir', 'reverse');
    set(ax1, 'ZDir', 'reverse');
    colormap(ax1, gray(256));
    view(ax1, 90, 0);

    %----------------------------------------------------------------------
    % Echo 1, axial slice
    %----------------------------------------------------------------------
    surf(ax2, reshape(X(idx,:,:) * 1e3, [N2 N3]), reshape(Y(idx,:,:) * 1e3, [N2 N3]), reshape(Z(idx,:,:) * 1e3, [N2 N3]), reshape(abs(img1(idx,:,:)), [N2 N3]), 'EdgeColor', 'none');
    set(ax2, 'Box', 'on', 'TickLabelInterpreter', 'latex', 'FontSize', FontSize);
    axis(ax2, 'image');
    xlabel(ax2, 'x [mm]', 'Interpreter', 'latex', 'FontSize', FontSize);
    ylabel(ax2, 'y [mm]', 'Interpreter', 'latex', 'FontSize', FontSize);
    zlabel(ax2, 'z [mm]', 'Interpreter', 'latex', 'FontSize', FontSize);
    title(ax2, {sprintf('Echo 1, axial slice %d', idx), sprintf('z = %3.1f mm', Z(idx,1,1) * 1e3)}, 'FontWeight', 'normal');
    colormap(ax2, gray(256));
    view(ax2, 0, 90);

    %----------------------------------------------------------------------
    % Echo 1, coronal slice
    %----------------------------------------------------------------------
    surf(ax3, reshape(X(:,idx,:) * 1e3, [N1 N3]), reshape(Y(:,idx,:) * 1e3, [N1 N3]), reshape(Z(:,idx,:) * 1e3, [N1 N3]), reshape(abs(img1(:,idx,:)), [N1 N3]), 'EdgeColor', 'none');
    set(ax3, 'Box', 'on', 'TickLabelInterpreter', 'latex', 'FontSize', FontSize);
    axis(ax3, 'image');
    xlabel(ax3, 'x [mm]', 'Interpreter', 'latex', 'FontSize', FontSize);
    ylabel(ax3, 'y [mm]', 'Interpreter', 'latex', 'FontSize', FontSize);
    zlabel(ax3, 'z [mm]', 'Interpreter', 'latex', 'FontSize', FontSize);
    title(ax3, {sprintf('Echo 1, coronal slice %d', idx), sprintf('y = %3.1f mm', Y(1,idx,1) * 1e3)}, 'FontWeight', 'normal');
    set(ax3, 'ZDir', 'reverse');
    colormap(ax3, gray(256));
    view(ax3, 0, 0);

    %----------------------------------------------------------------------
    % Echo 2, sagittal slice
    %----------------------------------------------------------------------
    surf(ax4, X(:,:,idx) * 1e3, Y(:,:,idx) * 1e3, Z(:,:,idx) * 1e3, abs(img2(:,:,idx)), 'EdgeColor', 'none');
    set(ax4, 'Box', 'on', 'TickLabelInterpreter', 'latex', 'FontSize', FontSize);
    axis(ax4, 'image');
    xlabel(ax4, 'x [mm]', 'Interpreter', 'latex', 'FontSize', FontSize);
    ylabel(ax4, 'y [mm]', 'Interpreter', 'latex', 'FontSize', FontSize);
    zlabel(ax4, 'z [mm]', 'Interpreter', 'latex', 'FontSize', FontSize);
    title(ax4, {sprintf('Echo 2, sagittal slice %d', idx), sprintf('x = %3.1f mm', X(1,1,idx) * 1e3)}, 'FontWeight', 'normal');
    set(ax4, 'YDir', 'reverse');
    set(ax4, 'ZDir', 'reverse');
    colormap(ax4, gray(256));
    view(ax4, 90, 0);

    %----------------------------------------------------------------------
    % Echo 2, axial slice
    %----------------------------------------------------------------------
    surf(ax5, reshape(X(idx,:,:) * 1e3, [N2 N3]), reshape(Y(idx,:,:) * 1e3, [N2 N3]), reshape(Z(idx,:,:) * 1e3, [N2 N3]), reshape(abs(img2(idx,:,:)), [N2 N3]), 'EdgeColor', 'none');
    set(ax5, 'Box', 'on', 'TickLabelInterpreter', 'latex', 'FontSize', FontSize);
    axis(ax5, 'image');
    xlabel(ax5, 'x [mm]', 'Interpreter', 'latex', 'FontSize', FontSize);
    ylabel(ax5, 'y [mm]', 'Interpreter', 'latex', 'FontSize', FontSize);
    zlabel(ax5, 'z [mm]', 'Interpreter', 'latex', 'FontSize', FontSize);
    title(ax5, {sprintf('Echo 2, axial slice %d', idx), sprintf('z = %3.1f mm', Z(idx,1,1) * 1e3)}, 'FontWeight', 'normal');
    colormap(ax5, gray(256));
    view(ax5, 0, 90);

    %----------------------------------------------------------------------
    % Echo 2, coronal slice
    %----------------------------------------------------------------------
    surf(ax6, reshape(X(:,idx,:) * 1e3, [N1 N3]), reshape(Y(:,idx,:) * 1e3, [N1 N3]), reshape(Z(:,idx,:) * 1e3, [N1 N3]), reshape(abs(img2(:,idx,:)), [N1 N3]), 'EdgeColor', 'none');
    set(ax6, 'Box', 'on', 'TickLabelInterpreter', 'latex', 'FontSize', FontSize);
    axis(ax6, 'image');
    xlabel(ax6, 'x [mm]', 'Interpreter', 'latex', 'FontSize', FontSize);
    ylabel(ax6, 'y [mm]', 'Interpreter', 'latex', 'FontSize', FontSize);
    zlabel(ax6, 'z [mm]', 'Interpreter', 'latex', 'FontSize', FontSize);
    title(ax6, {sprintf('Echo 2, coronal slice %d', idx), sprintf('y = %3.1f mm', Y(1,idx,1) * 1e3)}, 'FontWeight', 'normal');
    set(ax6, 'ZDir', 'reverse');
    colormap(ax6, gray(256));
    view(ax6, 0, 0);

    %----------------------------------------------------------------------
    % Full FOV
    %----------------------------------------------------------------------
    frame = getframe(gcf);
    frame.cdata = frame.cdata(10:940,120:1360,:);
    writeVideo(video_obj1, frame);

    %----------------------------------------------------------------------
    % Zoom
    %----------------------------------------------------------------------
    ylim(ax1, [-120 120]);
    zlim(ax1, [-120 120]);

    xlim(ax2, [-120 120]);
    ylim(ax2, [-120 120]);

    xlim(ax3, [-120 120]);
    zlim(ax3, [-120 120]);

    ylim(ax4, [-120 120]);
    zlim(ax4, [-120 120]);

    xlim(ax5, [-120 120]);
    ylim(ax5, [-120 120]);

    xlim(ax6, [-120 120]);
    zlim(ax6, [-120 120]);

    frame = getframe(gcf);
    frame.cdata = frame.cdata(10:940,120:1360,:);
    writeVideo(video_obj2, frame);
end

close(video_obj1);
close(video_obj2);
