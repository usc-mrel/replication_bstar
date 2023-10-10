% input_pulseq_bstar_protocol_20220727_1_38ms_1_60mm.m
% Written by Nam Gyun Lee
% Email: namgyunl@usc.edu, ggang56@gmail.com (preferred)
% Started: 08/21/2022, Last modified: 09/24/2022

%% Resolution
%--------------------------------------------------------------------------
% Common
%--------------------------------------------------------------------------
fov_read        = 340e-3; % FoV read [m]
base_resolution = 144;    % Base resolution
nr_radial_views = 31150;  % Radial views (i.e., total number of readouts)

% WASP: 354 half-projections per interleaf,    88 interleaves,  31152 radial views
% WASP: 350 half-projections per interleaf,    89 interleaves,  31150 radial views
% WASP: 600 half-projections per interleaf,   600 interleaves, 360000 radial views (free-breathing)

%% Contrast
%--------------------------------------------------------------------------
% Common
%--------------------------------------------------------------------------
flip_angle = 25; % Flip angle [deg]

%% Sequence
% bSTAR_1.84ms_1.01mm
% 208
% 1335 Hz/Px
% 140/750/1497.6/1.01
% 137.06

% bSTAR_1.70ms_1.35mm
% 160
% 1488 Hz/Px
% 140/680/1344.0/1.35
% 117.48

% bSTAR_1.55ms_1.56mm
% 144
% 1653 Hz/Px
% 140/610/1209.6/1.55
% 117.48

% bSTAR_1.38ms_1.60mm
% 144
% 1929 Hz/Px
% 140/520/1036.8/1.62
% 137.06

%--------------------------------------------------------------------------
% Part 1
%--------------------------------------------------------------------------
bandwidth = 1929; % Bandwidth [Hz/Px]

%--------------------------------------------------------------------------
% Part 2
%--------------------------------------------------------------------------
grad_mode = 'Fast';

%--------------------------------------------------------------------------
% Special
%--------------------------------------------------------------------------
rf_length        = 200e-6;        % RF length [sec]
rf_dummies       = 100;           % RF dummies
rf_phase         = 180;           % RF phase [deg]
trajectory       = 'Phyllotaxis';        % Trajectory: Archimedean (NOT USED!), WASP, CardioWASP (NOT USED!), Phyllotaxis
f_mod            = 3.6;           % F-mod (WASP) (NOT USED!)
tilt             = 0.7;           % Tilt (WASP) (NOT USED!)
nr_interleaves   = 89;            % number of interleaves
gradient_factor  = 0.80;          % Gradient factor (NOT USED!)
amplitude_factor = 2.00;          % Amplitude factor (NOT USED!)
spoiler_factor   = 0.50;          % Spoiler factor (NOT USED!)

%% Self-navigation
flag_self = 0; % self-navigation (SI projections): 0=no, 1=yes

%% Define misc parameters
readout_os_factor = 2; % readout oversampling factor

%--------------------------------------------------------------------------
% Geometry (Common)
%--------------------------------------------------------------------------
slices_per_slab = base_resolution; % Slices per slab

%--------------------------------------------------------------------------
% Sequence (Part 1)
%--------------------------------------------------------------------------
nr_echoes    = 2;          % Contrasts
readout_mode = 'Bipolar';  % Readout mode

%% Define main field strength [T]
B0 = 0.55; 

%% Define an output directory
output_directory = fullfile(pwd, mfilename);
