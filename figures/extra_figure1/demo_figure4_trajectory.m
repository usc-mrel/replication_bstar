% demo_figure4_trajectory.m
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

figure('Color', 'w', 'Position', [-6 520 1924 468]);
ColorOrder = get(gca, 'ColorOrder');

cmap = hsv(89);

ax1 = subplot(1,4,1);
ax2 = subplot(1,4,2);
ax3 = subplot(1,4,3);
ax4 = subplot(1,4,4);

az = -45;
el = 22;

idx = 1;

%--------------------------------------------------------------------------
% WASP 89
%--------------------------------------------------------------------------
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
hText1 = text(ax1, 0, 0, {'89 interleaves', sprintf('%d half-spokes per interleaf', 350)}, 'FontSize', FontSize, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
hText5 = text(ax1, 0, 0, {'WASP'}, 'FontSize', FontSize, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');

%--------------------------------------------------------------------------
% WASP 88
%--------------------------------------------------------------------------
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
hText2 = text(ax2, 0, 0, {'88 interleaves', sprintf('%d half-spokes per interleaf', 354)}, 'FontSize', FontSize, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');

%--------------------------------------------------------------------------
% SP 89
%--------------------------------------------------------------------------
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
hText3 = text(ax3, 0, 0, {'89 interleaves', sprintf('%d half-spokes per interleaf', 350)}, 'FontSize', FontSize, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
hText6 = text(ax3, 0, 0, {'Spiral phyllotaxis'}, 'FontSize', FontSize, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');

%--------------------------------------------------------------------------
% SP 88
%--------------------------------------------------------------------------
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
hText4 = text(ax4, 0, 0, {'88 interleaves', sprintf('%d half-spokes per interleaf', 354)}, 'FontSize', FontSize, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');

set([hXLabel1 hXLabel2 hXLabel3 hXLabel4], 'Position', [0.2605 -1.2348 -1.0478]);
set([hYLabel1 hYLabel2 hYLabel3 hYLabel4], 'Position', [-1.2685 0.2069 -0.9988]);

set([hText1 hText2 hText3 hText4], 'Position', [0.3989 0.3855 1.3728]);
set([hText5 hText6], 'Position', [1.8450 -0.8310 2.1000]);

set(ax1, 'Position', [0.1300 0.1100-0.1 0.1550 0.8150]);
set(ax2, 'Position', [0.3361 0.1100-0.1 0.1550 0.8150]);
set(ax3, 'Position', [0.5422 0.1100-0.1 0.1550 0.8150]);
set(ax4, 'Position', [0.7484 0.1100-0.1 0.1550 0.8150]);

export_fig('extra_figure1', '-r500', '-tif');
