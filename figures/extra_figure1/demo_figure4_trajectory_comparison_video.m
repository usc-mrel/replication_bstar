% demo_figure4_trajectory_comparison_video.m
% Written by Nam Gyun Lee
% Email: namgyunl@usc.edu, ggang56@gmail.com (preferred)
% Started: 07/24/2023, Last modified: 07/24/2023

%% Clean slate
close all; clear all; clc;

%% Set source directories
package_directory = 'D:\replication_bstar';

%% Add source directories to search path
addpath(genpath(package_directory));
start_time = tic;

%% Calculate a wobbling Archimedean spiral pole (WASP) pattern (WASP 89)
%--------------------------------------------------------------------------
% WASP: 350 half-projections per interleaf, 89 interleaves, 31150
%--------------------------------------------------------------------------
nr_readouts    = 350;    % number of half-projections per interleaf
nr_interleaves = 89;     % number of interleaves
rpa_range      = [0 10]; % [degree]
[phi_wasp89, theta_wasp89] = calculate_wasp_pattern(nr_readouts, nr_interleaves, rpa_range);

%--------------------------------------------------------------------------
% Calculate spatial coordinates of 3D k-space trajectories
%--------------------------------------------------------------------------
x_wasp89 = sin(theta_wasp89) .* cos(phi_wasp89);
y_wasp89 = sin(theta_wasp89) .* sin(phi_wasp89);
z_wasp89 = cos(theta_wasp89);

x_wasp89 = reshape(x_wasp89, [nr_readouts nr_interleaves]);
y_wasp89 = reshape(y_wasp89, [nr_readouts nr_interleaves]);
z_wasp89 = reshape(z_wasp89, [nr_readouts nr_interleaves]);

%% Calculate a wobbling Archimedean spiral pole (WASP) pattern (WASP 88)
%--------------------------------------------------------------------------
% WASP: 354 half-projections per interleaf, 88 interleaves, 31152
%--------------------------------------------------------------------------
nr_readouts    = 354;    % number of half-projections per interleaf
nr_interleaves = 88;     % number of interleaves
rpa_range      = [0 10]; % [degree]
[phi_wasp88, theta_wasp88] = calculate_wasp_pattern(nr_readouts, nr_interleaves, rpa_range);

%--------------------------------------------------------------------------
% Calculate spatial coordinates of 3D k-space trajectories
%--------------------------------------------------------------------------
x_wasp88 = sin(theta_wasp88) .* cos(phi_wasp88);
y_wasp88 = sin(theta_wasp88) .* sin(phi_wasp88);
z_wasp88 = cos(theta_wasp88);

x_wasp88 = reshape(x_wasp88, [nr_readouts nr_interleaves]);
y_wasp88 = reshape(y_wasp88, [nr_readouts nr_interleaves]);
z_wasp88 = reshape(z_wasp88, [nr_readouts nr_interleaves]);

%% Calculate a spiral phyllotaxis pattern (SP 89)
%--------------------------------------------------------------------------
% SP: 350 half-projections per interleaf, 89 interleaves, 31150
%--------------------------------------------------------------------------
nr_readouts    = 350;    % number of half-projections per interleaf
nr_interleaves = 89;     % number of interleaves
[phi_sp89, theta_sp89] = calculate_spiral_phyllotaxis_pattern(nr_readouts, nr_interleaves, 0);

%--------------------------------------------------------------------------
% Calculate spatial coordinates of 3D k-space trajectories
%--------------------------------------------------------------------------
x_sp89 = sin(theta_sp89) .* cos(phi_sp89);
y_sp89 = sin(theta_sp89) .* sin(phi_sp89);
z_sp89 = cos(theta_sp89);

x_sp89 = reshape(x_sp89, [nr_readouts nr_interleaves]);
y_sp89 = reshape(y_sp89, [nr_readouts nr_interleaves]);
z_sp89 = reshape(z_sp89, [nr_readouts nr_interleaves]);

%% Calculate a spiral phyllotaxis pattern (SP 88)
%--------------------------------------------------------------------------
% SP: 354 half-projections per interleaf, 88 interleaves, 31152
%--------------------------------------------------------------------------
nr_readouts    = 354;    % number of half-projections per interleaf
nr_interleaves = 88;     % number of interleaves
[phi_sp88, theta_sp88] = calculate_spiral_phyllotaxis_pattern(nr_readouts, nr_interleaves, 0);

%--------------------------------------------------------------------------
% Calculate spatial coordinates of 3D k-space trajectories
%--------------------------------------------------------------------------
x_sp88 = sin(theta_sp88) .* cos(phi_sp88);
y_sp88 = sin(theta_sp88) .* sin(phi_sp88);
z_sp88 = cos(theta_sp88);

x_sp88 = reshape(x_sp88, [nr_readouts nr_interleaves]);
y_sp88 = reshape(y_sp88, [nr_readouts nr_interleaves]);
z_sp88 = reshape(z_sp88, [nr_readouts nr_interleaves]);

%% Display 3D k-space trajectories
FontSize = 20;
MarkerSize = 12;

figure('Color', 'w', 'Position', [-4 50 1924 468*2]);
ColorOrder = get(gca, 'ColorOrder');

ax1 = subplot(2,4,1);
ax2 = subplot(2,4,2);
ax3 = subplot(2,4,3);
ax4 = subplot(2,4,4);
ax5 = subplot(2,4,5); hold on;
ax6 = subplot(2,4,6); hold on;
ax7 = subplot(2,4,7); hold on;
ax8 = subplot(2,4,8); hold on;

az = -45;
el = 22;

video_fullpath = fullfile(pwd, sprintf('extra_video1'));
video_file = VideoWriter(video_fullpath, 'MPEG-4');
video_file.Quality = 100;

video_file.FrameRate = 1;
open(video_file);

for idx = 1:88
    %----------------------------------------------------------------------
    % WASP 89
    %----------------------------------------------------------------------
    plot3(ax1, x_wasp89(:,idx), y_wasp89(:,idx), z_wasp89(:,idx), '.-', 'Color', ColorOrder(mod(idx-1,7)+1,:), 'MarkerSize', MarkerSize);
    grid(ax1, 'On');
    axis(ax1, 'image');
    xlim(ax1, [-1 1]);
    ylim(ax1, [-1 1]);
    zlim(ax1, [-1 1]);
    view(ax1, az, el);
    set(ax1, 'Box', 'off', 'FontSize', 12);
    hXLabel1 = xlabel(ax1, 'k_x', 'FontSize', FontSize);
    hYLabel1 = ylabel(ax1, 'k_y', 'FontSize', FontSize);
    hZLabel1 = zlabel(ax1, 'k_z', 'FontSize', FontSize, 'Rotation', 0);
    hText1 = text(ax1, 0, 0, {'89 interleaves', sprintf('%d half-spokes per interleaf', 350), sprintf('Interleaf number = %d', idx)}, 'FontSize', FontSize, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
    hText5 = text(ax1, 0, 0, {'WASP'}, 'FontSize', FontSize, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');

    %----------------------------------------------------------------------
    % WASP 88
    %----------------------------------------------------------------------
    plot3(ax2, x_wasp88(:,idx), y_wasp88(:,idx), z_wasp88(:,idx), '.-', 'Color', ColorOrder(mod(idx-1,7)+1,:), 'MarkerSize', MarkerSize);
    grid(ax2, 'On');
    axis(ax2, 'image');
    xlim(ax2, [-1 1]);
    ylim(ax2, [-1 1]);
    zlim(ax2, [-1 1]);
    view(ax2, az, el);
    set(ax2, 'Box', 'off', 'FontSize', 12);
    hXLabel2 = xlabel(ax2, 'k_x', 'FontSize', FontSize);
    hYLabel2 = ylabel(ax2, 'k_y', 'FontSize', FontSize);
    hZLabel2 = zlabel(ax2, 'k_z', 'FontSize', FontSize, 'Rotation', 0);
    hText2 = text(ax2, 0, 0, {'88 interleaves', sprintf('%d half-spokes per interleaf', 354), sprintf('Interleaf number = %d', idx)}, 'FontSize', FontSize, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');

    %----------------------------------------------------------------------
    % SP 89
    %----------------------------------------------------------------------
    plot3(ax3, x_sp89(:,idx), y_sp89(:,idx), z_sp89(:,idx), '.-', 'Color', ColorOrder(mod(idx-1,7)+1,:), 'MarkerSize', MarkerSize);
    grid(ax3, 'On');
    axis(ax3, 'image');
    xlim(ax3, [-1 1]);
    ylim(ax3, [-1 1]);
    zlim(ax3, [-1 1]);
    view(ax3, az, el);
    set(ax3, 'Box', 'off', 'FontSize', 12);
    hXLabel3 = xlabel(ax3, 'k_x', 'FontSize', FontSize);
    hYLabel3 = ylabel(ax3, 'k_y', 'FontSize', FontSize);
    hZLabel3 = zlabel(ax3, 'k_z', 'FontSize', FontSize, 'Rotation', 0);
    hText3 = text(ax3, 0, 0, {'89 interleaves', sprintf('%d half-spokes per interleaf', 350), sprintf('Interleaf number = %d', idx)}, 'FontSize', FontSize, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
    hText6 = text(ax3, 0, 0, {'Spiral phyllotaxis'}, 'FontSize', FontSize, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');

    %----------------------------------------------------------------------
    % SP 88
    %----------------------------------------------------------------------
    plot3(ax4, x_sp88(:,idx), y_sp88(:,idx), z_sp88(:,idx), '.-', 'Color', ColorOrder(mod(idx-1,7)+1,:), 'MarkerSize', MarkerSize);
    grid(ax4, 'On');
    axis(ax4, 'image');
    xlim(ax4, [-1 1]);
    ylim(ax4, [-1 1]);
    zlim(ax4, [-1 1]);
    view(ax4, az, el);
    set(ax4, 'Box', 'off', 'FontSize', 12);
    hXLabel4 = xlabel(ax4, 'k_x', 'FontSize', FontSize);
    hYLabel4 = ylabel(ax4, 'k_y', 'FontSize', FontSize);
    hZLabel4 = zlabel(ax4, 'k_z', 'FontSize', FontSize, 'Rotation', 0);
    hText4 = text(ax4, 0, 0, {'88 interleaves', sprintf('%d half-spokes per interleaf', 354), sprintf('Interleaf number = %d', idx)}, 'FontSize', FontSize, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');

    %----------------------------------------------------------------------
    % WASP 89 (accumulated)
    %----------------------------------------------------------------------
    plot3(ax5, x_wasp89(:,idx), y_wasp89(:,idx), z_wasp89(:,idx), '.-', 'Color', ColorOrder(mod(idx-1,7)+1,:), 'MarkerSize', MarkerSize);
    grid(ax5, 'On');
    axis(ax5, 'image');
    xlim(ax5, [-1 1]);
    ylim(ax5, [-1 1]);
    zlim(ax5, [-1 1]);
    view(ax5, az, el);
    set(ax5, 'Box', 'off', 'FontSize', 12);
    hXLabel5 = xlabel(ax5, 'k_x', 'FontSize', FontSize);
    hYLabel5 = ylabel(ax5, 'k_y', 'FontSize', FontSize);
    hZLabel5 = zlabel(ax5, 'k_z', 'FontSize', FontSize, 'Rotation', 0);
    title(ax5, {sprintf('Interleaves from 1 to %d', idx)}, 'FontWeight', 'normal', 'FontSize', 18);

    %----------------------------------------------------------------------
    % WASP 88 (accumulated)
    %----------------------------------------------------------------------
    plot3(ax6, x_wasp88(:,idx), y_wasp88(:,idx), z_wasp88(:,idx), '.-', 'Color', ColorOrder(mod(idx-1,7)+1,:), 'MarkerSize', MarkerSize);
    grid(ax6, 'On');
    axis(ax6, 'image');
    xlim(ax6, [-1 1]);
    ylim(ax6, [-1 1]);
    zlim(ax6, [-1 1]);
    view(ax6, az, el);
    set(ax6, 'Box', 'off', 'FontSize', 12);
    hXLabel6 = xlabel(ax6, 'k_x', 'FontSize', FontSize);
    hYLabel6 = ylabel(ax6, 'k_y', 'FontSize', FontSize);
    hZLabel6 = zlabel(ax6, 'k_z', 'FontSize', FontSize, 'Rotation', 0);
    title(ax6, {sprintf('Interleaves from 1 to %d', idx)}, 'FontWeight', 'normal', 'FontSize', 18);

    %----------------------------------------------------------------------
    % SP 89 (accumulated)
    %----------------------------------------------------------------------
    plot3(ax7, x_sp89(:,idx), y_sp89(:,idx), z_sp89(:,idx), '.-', 'Color', ColorOrder(mod(idx-1,7)+1,:), 'MarkerSize', MarkerSize);
    grid(ax7, 'On');
    axis(ax7, 'image');
    xlim(ax7, [-1 1]);
    ylim(ax7, [-1 1]);
    zlim(ax7, [-1 1]);
    view(ax7, az, el);
    set(ax7, 'Box', 'off', 'FontSize', 12);
    hXLabel7 = xlabel(ax7, 'k_x', 'FontSize', FontSize);
    hYLabel7 = ylabel(ax7, 'k_y', 'FontSize', FontSize);
    hZLabel7 = zlabel(ax7, 'k_z', 'FontSize', FontSize, 'Rotation', 0);
    title(ax7, {sprintf('Interleaves from 1 to %d', idx)}, 'FontWeight', 'normal', 'FontSize', 18);

    %----------------------------------------------------------------------
    % SP 88 (accumulated)
    %----------------------------------------------------------------------
    plot3(ax8, x_sp88(:,idx), y_sp88(:,idx), z_sp88(:,idx), '.-', 'Color', ColorOrder(mod(idx-1,7)+1,:), 'MarkerSize', MarkerSize);
    grid(ax8, 'On');
    axis(ax8, 'image');
    xlim(ax8, [-1 1]);
    ylim(ax8, [-1 1]);
    zlim(ax8, [-1 1]);
    view(ax8, az, el);
    set(ax8, 'Box', 'off', 'FontSize', 12);
    hXLabel8 = xlabel(ax8, 'k_x', 'FontSize', FontSize);
    hYLabel8 = ylabel(ax8, 'k_y', 'FontSize', FontSize);
    hZLabel8 = zlabel(ax8, 'k_z', 'FontSize', FontSize, 'Rotation', 0);
    title(ax8, {sprintf('Interleaves from 1 to %d', idx)}, 'FontWeight', 'normal', 'FontSize', 18);

    set([hXLabel1 hXLabel2 hXLabel3 hXLabel4 hXLabel5 hXLabel6 hXLabel7 hXLabel8], 'Position', [0.2605 -1.2348 -1.0478]);
    set([hYLabel1 hYLabel2 hYLabel3 hYLabel4 hYLabel5 hYLabel6 hYLabel7 hYLabel8], 'Position', [-1.2685 0.2069 -0.9988]);

    set(ax1, 'Position', [0.1300 0.5838-0.05*2 0.1550 0.3412]);
    set(ax2, 'Position', [0.3361 0.5838-0.05*2 0.1550 0.3412]);
    set(ax3, 'Position', [0.5422 0.5838-0.05*2 0.1550 0.3412]);
    set(ax4, 'Position', [0.7484 0.5838-0.05*2 0.1550 0.3412]);
    set(ax5, 'Position', [0.1300 0.1100-0.05 0.1550 0.3412]);
    set(ax6, 'Position', [0.3361 0.1100-0.05 0.1550 0.3412]);
    set(ax7, 'Position', [0.5422 0.1100-0.05 0.1550 0.3412]);
    set(ax8, 'Position', [0.7484 0.1100-0.05 0.1550 0.3412]);

    set([hText1 hText2 hText3 hText4], 'Position', [0.3989 0.3855 1.3728]);
    set([hText5 hText6], 'Position', [1.8450 -0.8310 2.4000]);

    frame = getframe(gcf);
    frame.cdata = frame.cdata(28:900,190:1750,:);
    writeVideo(video_file, frame);
end
close(video_file);
