% input_pulseq_bstar_protocol_20220727_1_38ms_1_60mm.m
% Written by Nam Gyun Lee
% Email: namgyunl@usc.edu, ggang56@gmail.com (preferred)
% Started: 08/21/2022, Last modified: 09/24/2022

%% Resolution
%--------------------------------------------------------------------------
% Common
%--------------------------------------------------------------------------
fov_read        = 340e-3;     % FoV read [m]
base_resolution = 224;        % Base resolution
nr_radial_views = 600 * 600;  % Radial views (i.e., total number of readouts)

%% Contrast
%--------------------------------------------------------------------------
% Common
%--------------------------------------------------------------------------
flip_angle = 25; % Flip angle [deg]

%% Sequence
% bSTAR_2.14ms_0.90mm
% 224
% 1116 Hz/Px
% 140/900/1792.0/0.90mm
% 123.35
% 17.27 mT/m

%--------------------------------------------------------------------------
% Part 1
%--------------------------------------------------------------------------
bandwidth = 1116; % Bandwidth [Hz/Px]

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
trajectory       = 'WASP';        % Trajectory: Archimedean, WASP, CardioWASP, Phyllotaxis
f_mod            = 3.6;           % F-mod (WASP) (NOT USED!)
tilt             = 0.7;           % Tilt (WASP) (NOT USED!)
nr_interleaves   = 600;           % number of interleaves
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
