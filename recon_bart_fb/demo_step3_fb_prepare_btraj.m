% demo_step3_fb_prepare_btraj.m
% Written by Nam Gyun Lee
% Email: namgyunl@usc.edu, ggang56@gmail.com (preferred)
% Started: 04/25/2024, Last modified: 09/18/2024

%% Clean slate
close all; clearvars -except json_number nr_json_files json_files json_file; clc;

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
    siemens_twix_file  = strrep(json.siemens_twix_file, '/', '\');
    ismrmrd_data_file  = strrep(json.ismrmrd_data_file, '/', '\');
    ismrmrd_traj_file  = strrep(json.ismrmrd_traj_file, '/', '\');
    trj_file           = strrep(json.trj_file, '/', '\');
    output_parent_path = strrep(json.output_path, '/', '\');
else
    siemens_twix_file  = json.siemens_twix_file;
    ismrmrd_data_file  = json.ismrmrd_data_file;
    ismrmrd_traj_file  = json.ismrmrd_traj_file;
    trj_file           = json.trj_file;
    output_parent_path = json.output_path;
end

%--------------------------------------------------------------------------
% Reconstruction parameters
%--------------------------------------------------------------------------
recon_matrix_size   = json.recon_conf.recon_matrix_size.';   % reconstruction matrix size
dicom_matrix_size   = json.recon_conf.dicom_matrix_size.';   % DICOM matrix size
recon_interp_factor = json.recon_conf.recon_interp_factor.'; % reconstruction interpolation factor: recon_resolution = encoded_resolution / recon_interp_factor;
cal_size            = json.recon_conf.cal_size.';            % size of the calibration region
lambda              = json.recon_conf.lambda;                % l1 regularization parameter for pics in BART
max_iter            = json.recon_conf.max_iter;              % max. number of iterations for pics in BART
B                   = json.recon_conf.B;                     % number of frames (respiratory motion states)
kappa               = json.recon_conf.kappa;                 % kappa: kappa=0 => W=I

%--------------------------------------------------------------------------
% Respiratory filtering parameters
%--------------------------------------------------------------------------
n  = json.resp_conf.n;
fl = json.resp_conf.fl;
fh = json.resp_conf.fh;
fw = json.resp_conf.fw;

%--------------------------------------------------------------------------
% Define the number of echoes
%--------------------------------------------------------------------------
nr_echoes = 2;

%% Define an output path
output_path = fullfile(output_parent_path, sprintf('%d_motion_states', B));

%% Make an output path
mkdir(output_path);

%% Read an ISMRMRD file
%--------------------------------------------------------------------------
% k-space trajectories
%--------------------------------------------------------------------------
tstart = tic; fprintf('%s: Reading an ISMRMRD file: %s... ', datetime, ismrmrd_traj_file);
if exist(ismrmrd_traj_file, 'file')
    traj_dset = ismrmrd.Dataset(ismrmrd_traj_file, 'dataset');
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
else
    error('File %s does not exist.  Please generate it.' , ismrmrd_traj_file);
end

%--------------------------------------------------------------------------
% k-space data
%--------------------------------------------------------------------------
tstart = tic; fprintf('%s: Reading an ISMRMRD file: %s... ', datetime, ismrmrd_data_file);
if exist(ismrmrd_data_file, 'file')
    data_dset = ismrmrd.Dataset(ismrmrd_data_file, 'dataset');
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
else
    error('File %s does not exist.  Please generate it.' , ismrmrd_data_file);
end

%% Read acquisitions
tstart = tic; fprintf('%s: Reading acquisitions (traj)... ', datetime);
raw_traj = traj_dset.readAcquisition(); % read all the acquisitions
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%--------------------------------------------------------------------------
% Get data type
%--------------------------------------------------------------------------
acq_is_navigation_data = raw_traj.head.flagIsSet('ACQ_IS_NAVIGATION_DATA'); % FID navigator
acq_user1              = raw_traj.head.flagIsSet('ACQ_USER1'); % self-navigation

%% Get k-space trajectories for each data type
nav_traj  = raw_traj.select(find(acq_is_navigation_data));
self_traj = raw_traj.select(find(acq_user1));
img_traj  = raw_traj.select(find(~(acq_is_navigation_data | acq_user1)));

%% Get imaging parameters from an XML header
traj_header = ismrmrd.xml.deserialize(traj_dset.readxml);
data_header = ismrmrd.xml.deserialize(data_dset.readxml);

%--------------------------------------------------------------------------
% measurement information
%--------------------------------------------------------------------------
patient_position = data_header.measurementInformation.patientPosition;

%--------------------------------------------------------------------------
% Encoding Space (Nkx, Nky, Nkz)
%--------------------------------------------------------------------------
encoded_fov(1) = traj_header.encoding.encodedSpace.fieldOfView_mm.x * 1e-3; % [m] RO
encoded_fov(2) = traj_header.encoding.encodedSpace.fieldOfView_mm.y * 1e-3; % [m] PE
encoded_fov(3) = traj_header.encoding.encodedSpace.fieldOfView_mm.z * 1e-3; % [m] SL

Nkx = traj_header.encoding.encodedSpace.matrixSize.x; % number of readout samples in k-space
Nky = traj_header.encoding.encodedSpace.matrixSize.y; % number of phase-encoding steps in k-space
Nkz = traj_header.encoding.encodedSpace.matrixSize.z; % number of slice-encoding steps in k-space

encoded_resolution = encoded_fov ./ [Nkx Nky Nkz]; % [m]

recon_resolution = encoded_resolution;

%--------------------------------------------------------------------------
% self_navigation_flag
%--------------------------------------------------------------------------
self_navigation_flag = traj_header.encoding.trajectoryDescription.userParameterLong(4).value;

%--------------------------------------------------------------------------
% TR
%--------------------------------------------------------------------------
TR = data_header.sequenceParameters.TR * 1e-3; % [msec] * [sec/1e3msec] => [sec]

%% Read a Siemens .dat file
fprintf('%s: Reading a Siemens .dat file: %s\n', datetime, siemens_twix_file);
twix = mapVBVD(siemens_twix_file);
if length(twix) > 1
    twix = twix{end};
end

%% Get a slice normal vector from Siemens TWIX format
%--------------------------------------------------------------------------
% dNormalSag: Sagittal component of a slice normal vector (in the PCS)
%--------------------------------------------------------------------------
if isfield(twix.hdr.MeasYaps.sSliceArray.asSlice{1}.sNormal, 'dSag')
    dNormalSag = twix.hdr.MeasYaps.sSliceArray.asSlice{1}.sNormal.dSag;
else
    dNormalSag = 0;
end

%--------------------------------------------------------------------------
% dNormalCor: Coronal component of a slice normal vector (in the PCS)
%--------------------------------------------------------------------------
if isfield(twix.hdr.MeasYaps.sSliceArray.asSlice{1}.sNormal, 'dCor')
    dNormalCor = twix.hdr.MeasYaps.sSliceArray.asSlice{1}.sNormal.dCor;
else
    dNormalCor = 0;
end

%--------------------------------------------------------------------------
% dNormalTra: Transverse component of a slice normal vector (in the PCS)
%--------------------------------------------------------------------------
if isfield(twix.hdr.MeasYaps.sSliceArray.asSlice{1}.sNormal, 'dTra')
    dNormalTra = twix.hdr.MeasYaps.sSliceArray.asSlice{1}.sNormal.dTra;
else
    dNormalTra = 0;
end

%--------------------------------------------------------------------------
% dRotAngle: Slice rotation angle ("swap Fre/Pha")
%--------------------------------------------------------------------------
if isfield(twix.hdr.MeasYaps.sSliceArray.asSlice{1}, 'dInPlaneRot')
    dRotAngle = twix.hdr.MeasYaps.sSliceArray.asSlice{1}.dInPlaneRot; % [rad]
else
    dRotAngle = 0; % [rad]
end

%% Get a slice offset in the PCS from Siemens TWIX format
if isfield(twix.hdr.MeasYaps.sSliceArray.asSlice{1}, 'sPosition')
    if isfield(twix.hdr.MeasYaps.sSliceArray.asSlice{1}.sPosition, 'dSag')
        sag_offset_twix = twix.hdr.MeasYaps.sSliceArray.asSlice{1}.sPosition.dSag; % [mm]
    else
        sag_offset_twix = 0; % [mm]
    end
    if isfield(twix.hdr.MeasYaps.sSliceArray.asSlice{1}.sPosition, 'dCor')
        cor_offset_twix = twix.hdr.MeasYaps.sSliceArray.asSlice{1}.sPosition.dCor; % [mm]
    else
        cor_offset_twix = 0; % [mm]
    end
    if isfield(twix.hdr.MeasYaps.sSliceArray.asSlice{1}.sPosition, 'dTra')
        tra_offset_twix = twix.hdr.MeasYaps.sSliceArray.asSlice{1}.sPosition.dTra; % [mm]
    else
        tra_offset_twix = 0; % [mm]
    end
else
    sag_offset_twix = 0; % [mm]
    cor_offset_twix = 0; % [mm]
    tra_offset_twix = 0; % [mm]
end

%% Set a slice offset in the PCS [m]
pcs_offset = [sag_offset_twix; cor_offset_twix; tra_offset_twix] * 1e-3; % [mm] * [m/1e3mm] => [m]

%% Calculate a transformation matrix from the RCS to the GCS [r,c,s] <=> [PE,RO,SL]
R_rcs2gcs = [0    1    0 ; % [PE]   [0 1 0] * [r]
             1    0    0 ; % [RO] = [1 0 0] * [c]
             0    0    1]; % [SL]   [0 0 1] * [s]

%% Calculate a rotation matrix from the GCS to the PCS
[R_gcs2pcs, phase_sign, read_sign, main_orientation] = siemens_calculate_transform_gcs_to_pcs(dNormalSag, dNormalCor, dNormalTra, dRotAngle);

%% Calculate a rotation matrix from the PCS to the DCS
R_pcs2dcs = siemens_calculate_transform_pcs_to_dcs(patient_position);

%% Calculate a rotation matrix from the GCS to the DCS
R_gcs2dcs = R_pcs2dcs * R_gcs2pcs;

%% Calculate a rotation matrix from the RCS to the DCS
R_rcs2dcs = R_pcs2dcs * R_gcs2pcs * R_rcs2gcs;

%% Calculate a slice offset in the DCS [m]
dcs_offset = R_pcs2dcs * pcs_offset; % 3 x 1

%% Display slice information
fprintf('======================= SLICE INFORMATION ========================\n');
fprintf('main_orientation = %d (SAGITTAL/CORONAL/TRANSVERSAL = 0/1/2)\n', main_orientation);
fprintf('dNormalSag = %+g \ndNormalCor = %+g \ndNormalTra = %+g \ndRotAngle = %g [rad]\n', dNormalSag, dNormalCor, dNormalTra, dRotAngle);
fprintf('phase_sign = %+g, read_sign = %+g\n', phase_sign, read_sign);
fprintf('---------------------- From Siemens TWIX format ------------------\n');
fprintf('                   [sag]   %10.5f [mm]\n', sag_offset_twix);
fprintf('slice offset     : [cor] = %10.5f [mm]\n', cor_offset_twix);
fprintf('                   [tra]   %10.5f [mm]\n', tra_offset_twix);
fprintf('---------------------- From Siemens TWIX format ------------------\n');
fprintf('                   [sag]   [%10.5f %10.5f %10.5f][PE]\n', R_gcs2pcs(1,1), R_gcs2pcs(1,2), R_gcs2pcs(1,3));
fprintf('R_gcs2pcs        : [cor] = [%10.5f %10.5f %10.5f][RO]\n', R_gcs2pcs(2,1), R_gcs2pcs(2,2), R_gcs2pcs(2,3));
fprintf('                   [tra]   [%10.5f %10.5f %10.5f][SL]\n', R_gcs2pcs(3,1), R_gcs2pcs(3,2), R_gcs2pcs(3,3));
fprintf('------------------------------------------------------------------\n');
fprintf('                   [ x ]   [%10.5f %10.5f %10.5f][sag]\n', R_pcs2dcs(1,1), R_pcs2dcs(1,2), R_pcs2dcs(1,3));
fprintf('R_pcs2dcs        : [ y ] = [%10.5f %10.5f %10.5f][cor]\n', R_pcs2dcs(2,1), R_pcs2dcs(2,2), R_pcs2dcs(2,3));
fprintf('                   [ z ]   [%10.5f %10.5f %10.5f][tra]\n', R_pcs2dcs(3,1), R_pcs2dcs(3,2), R_pcs2dcs(3,3));
fprintf('------------------------------------------------------------------\n');
fprintf('                   [ x ]   [%10.5f %10.5f %10.5f][PE]\n', R_gcs2dcs(1,1), R_gcs2dcs(1,2), R_gcs2dcs(1,3));
fprintf('R_gcs2dcs        : [ y ] = [%10.5f %10.5f %10.5f][RO]\n', R_gcs2dcs(2,1), R_gcs2dcs(2,2), R_gcs2dcs(2,3));
fprintf('                   [ z ]   [%10.5f %10.5f %10.5f][SL]\n', R_gcs2dcs(3,1), R_gcs2dcs(3,2), R_gcs2dcs(3,3));
fprintf('==================================================================\n');

%% Parse an ISMRMRD header
adc_samples  = double(max(img_traj.head.number_of_samples));
discard_pre  = double(max(img_traj.head.discard_pre));
discard_post = double(max(img_traj.head.discard_post));

%--------------------------------------------------------------------------
% Number of half-spokes per interleaf
%--------------------------------------------------------------------------
nr_readouts = traj_header.encoding.encodingLimits.kspace_encoding_step_1.maximum + 1;

%--------------------------------------------------------------------------
% Number of interleaves
%--------------------------------------------------------------------------
nr_interleaves = traj_header.encoding.encodingLimits.segment.maximum + 1;

%--------------------------------------------------------------------------
% Define the total number of half-spokes
%--------------------------------------------------------------------------
nr_radial_views = nr_readouts * nr_interleaves;

%--------------------------------------------------------------------------
% Set N1, N2, and N3
%--------------------------------------------------------------------------
N1 = recon_matrix_size(1); % number of samples in the row (r) direction of the RCS
N2 = recon_matrix_size(2); % number of samples in the col (c) direction of the RCS
N3 = recon_matrix_size(3); % number of samples in the slice (s) direction of the RCS

%--------------------------------------------------------------------------
% Calculate the total number of voxels
%--------------------------------------------------------------------------
N = N1 * N2 * N3;

%% Calculate a scaling matrix [m]
scaling_matrix = diag(recon_resolution); % [m]

%% Calculate spatial coordinates in the RCS [m]
[I1,I2,I3] = ndgrid((1:N1).', (1:N2).', (1:N3).');
r_rcs = scaling_matrix * cat(1, I1(:).' - (floor(N1/2) + 1), I2(:).' - (floor(N2/2) + 1), I3(:).' - (floor(N3/2) + 1)); % 3 x N

%% Calculate spatial coordinates in the DCS [m]
r_dcs = repmat(dcs_offset, [1 N]) + R_pcs2dcs * R_gcs2pcs * R_rcs2gcs * r_rcs; % 3 x N

%% Write a .cfl file
%--------------------------------------------------------------------------
% x (N1 x N2 x N3)
%--------------------------------------------------------------------------
cfl_file = fullfile(output_parent_path, 'x');
tstart = tic; fprintf('%s: Writing a .cfl file: %s... ', datetime, cfl_file);
writecfl(cfl_file, reshape(r_dcs(1,:), [N1 N2 N3]));
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%--------------------------------------------------------------------------
% y (N1 x N2 x N3)
%--------------------------------------------------------------------------
cfl_file = fullfile(output_parent_path, 'y');
tstart = tic; fprintf('%s: Writing a .cfl file: %s... ', datetime, cfl_file);
writecfl(cfl_file, reshape(r_dcs(2,:), [N1 N2 N3]));
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%--------------------------------------------------------------------------
% z (N1 x N2 x N3)
%--------------------------------------------------------------------------
cfl_file = fullfile(output_parent_path, 'z');
tstart = tic; fprintf('%s: Writing a .cfl file: %s... ', datetime, cfl_file);
writecfl(cfl_file, reshape(r_dcs(3,:), [N1 N2 N3]));
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%% Read a .trj file [rad/m]
%--------------------------------------------------------------------------
% trj(:,1)  - +X, echo 1 | trj(:,7)  - +X, echo 2
% trj(:,2)  - -X, echo 1 | trj(:,8)  - -X, echo 2
% trj(:,3)  - +Y, echo 1 | trj(:,9)  - +Y, echo 2
% trj(:,4)  - -Y, echo 1 | trj(:,10) - -Y, echo 2
% trj(:,5)  - +Z, echo 1 | trj(:,11) - +Z, echo 2
% trj(:,6)  - -Z, echo 1 | trj(:,12) - -Z, echo 2
%--------------------------------------------------------------------------
[trj_header,trj] = read_trj_file(trj_file); % [rad/m]

%% Get the parameterization (phi,theta) of 3D radial k-space trajectories
phi   = zeros(nr_radial_views, 1, 'double');
theta = zeros(nr_radial_views, 1, 'double');

for i = 1:nr_radial_views
    phi(i) = img_traj.data{i}(1);
    theta(i) = img_traj.data{i}(2);
end

%% Calculate measured k-space trajectories in the DCS [rad/m] [x,y,z]
k_dcs = zeros(3, adc_samples, nr_radial_views, 'double');

for i = 1:nr_radial_views
    %----------------------------------------------------------------------
    % Calculate a rotation matrix
    % theta: polar angle     (angle from RO to SL)
    % phi  : azimuthal angle (angle from SL to PE)
    %
    %         ^ RO(2)                   ^ RO(2)                  ^ RO(2)
    %         |                  theta  |                        |
    %         |                       \ |                        | /|
    %         |                        \|                        |/ |
    %         +---------> PE(1) =>      +---------> PE(1)        +--|------> PE(1)
    %        /                         /                        / \ |
    %       /                         /                        /phi\|
    %      v SL(3)                   v SL(3)                  v SL(3)
    %
    % All rotation matrices are right-handed.
    % 1. Apply rotation about the PE direction by theta
    % 2. Apply rotation about the RO direction by phi
    %
    % R_PE(theta) = [1     0           0     ], R_RO(phi) = [ cos(phi) 0 sin(phi)]
    %               [0 cos(theta) -sin(theta)]              [    0     1     0   ]
    %               [0 sin(theta)  cos(theta)]              [-sin(phi) 0 cos(phi)]
    %----------------------------------------------------------------------
    cos_theta = cos(theta(i));
    sin_theta = sin(theta(i));
    cos_phi   = cos(phi(i));
    sin_phi   = sin(phi(i));

    %----------------------------------------------------------------------
    % Calculate the rotation matrix about the PE direction by theta
    %----------------------------------------------------------------------
    R_PE = [1        0             0     ;
            0    cos_theta    -sin_theta ;
            0    sin_theta     cos_theta];

    %----------------------------------------------------------------------
    % Calculate the rotation matrix about the RO direction by phi
    %----------------------------------------------------------------------
    R_RO = [ cos_phi    0    sin_phi ;
                0       1       0    ;
            -sin_phi    0    cos_phi];

    %----------------------------------------------------------------------
    % Define an initial vector (along the RO direction, PRS)
    %----------------------------------------------------------------------
    initial_vector = [0; 1; 0]; % [unitless]

    %----------------------------------------------------------------------
    % Calculate the Cartesian coordinates of a rotated vector in the DCS
    %----------------------------------------------------------------------
    xyz = R_gcs2dcs * R_RO * R_PE * initial_vector;

    %----------------------------------------------------------------------
    % Calculate measured k-space trajectories in the DCS [rad/m]
    % trj(:,1)  - +X, echo 1 | trj(:,7)  - +X, echo 2
    % trj(:,2)  - -X, echo 1 | trj(:,8)  - -X, echo 2
    % trj(:,3)  - +Y, echo 1 | trj(:,9)  - +Y, echo 2
    % trj(:,4)  - -Y, echo 1 | trj(:,10) - -Y, echo 2
    % trj(:,5)  - +Z, echo 1 | trj(:,11) - +Z, echo 2
    % trj(:,6)  - -Z, echo 1 | trj(:,12) - -Z, echo 2
    %----------------------------------------------------------------------
    if xyz(1) > 0
        k_dcs(1,:,i) = xyz(1) * cat(1, trj(:,1), trj(:,7)); % +X (echo 1), +X (echo 2)
    else
        k_dcs(1,:,i) = xyz(1) * cat(1, trj(:,2), trj(:,8)); % -X (echo 1), -X (echo 2)
    end

    if xyz(2) > 0
        k_dcs(2,:,i) = xyz(2) * cat(1, trj(:,3), trj(:,9));  % +Y (echo 1), +Y (echo 2)
    else
        k_dcs(2,:,i) = xyz(2) * cat(1, trj(:,4), trj(:,10)); % -Y (echo 1), -Y (echo 2)
    end

    if xyz(3) > 0
        k_dcs(3,:,i) = xyz(3) * cat(1, trj(:,5), trj(:,11)); % +Z (echo 1), +Z (echo 2)
    else
        k_dcs(3,:,i) = xyz(3) * cat(1, trj(:,6), trj(:,12)); % -Z (echo 1), -Z (echo 2)
    end
end

%% Calculate measured k-space trajectories in the GCS [rad/m] [PE,RO,SL]
k_gcs = zeros(3, adc_samples, nr_radial_views, 'double');
for i = 1:nr_radial_views
    k_gcs(:,:,i) = R_gcs2dcs.' * k_dcs(:,:,i); % 3 x adc_samples
end

%% Change the sign of k-space trajectories in the GCS [rad/m] [PE,RO,SL]
k_gcs(1,:,:) = phase_sign * k_gcs(1,:,:); % PE (ku)
k_gcs(2,:,:) = read_sign  * k_gcs(2,:,:); % RO (kv)

%% Calculate measured k-space trajectories in the RCS [rad/m] [r,c,s]
k_rcs = zeros(3, adc_samples, nr_radial_views, 'double');
for i = 1:nr_radial_views
    k_rcs(:,:,i) = R_rcs2gcs.' * k_gcs(:,:,i); % 3 x adc_samples
end

%% Calculate the maximum k-space value [rad/m]
krmax = (2 * pi / recon_resolution(1)) / 2; % [rad/m]

%% Calculate "unscaled" k-space trajectories for BART (3 x adc_samples x nr_radial_views)
%--------------------------------------------------------------------------
% A minus is required to perform forward FFT to move from k-space to image space
% BART requires "traj" scaled to [-0.5,0.5] * N
% Here, we calculate "unscaled" traj scaled to [-0.5,0.5]
%--------------------------------------------------------------------------
traj_unscaled = zeros(3, adc_samples, nr_radial_views, 'single');
traj_unscaled(1,:,:) = -k_rcs(1,:,:) / (2 * krmax);
traj_unscaled(2,:,:) = -k_rcs(2,:,:) / (2 * krmax); % AP direction for SAGITTAL with PE along AP
traj_unscaled(3,:,:) = -k_rcs(3,:,:) / (2 * krmax);

%% Read a .cfl file
cfl_file = fullfile(output_parent_path, sprintf('resp_n%d_fl%4.2f_fh%4.2f_fw%4.2f', n, fl, fh, fw));
tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
resp = real(readcfl(cfl_file));
fprintf('%s: done! (%6.4f/%6.4f sec)\n', datetime, toc(tstart), toc(start_time));

%% Perform data binning for respiratory motion
margin = 10;
bins = prctile(resp, linspace(0 + margin, 100 - margin, B + 1).');
bidx = cell(B,1);
for b = 1:B
    idx = find((resp >= bins(b)) & (resp < bins(b + 1)));
    bidx{b} = idx;
end

% %%
% nr_radial_views = size(resp);
% TR = 2.3e-3;
% t = (0:nr_radial_views-1).' * TR;
% 
% diff_resp = diff(cat(1, resp, resp(end)));
% positive_slope = (diff_resp > 0);
% nagative_slope = (diff_resp <= 0);
% 
% figure;
% hold on;
% plot(t(positive_slope), resp(positive_slope), '.');
% plot(t(nagative_slope), resp(nagative_slope), '.');
% 
% 
% 
% %%
% resp1 = resp(bidx{1});
% 
% %%
% figure;
% plot(resp(bidx{1}))

%% Write a .cfl file
adc_samples_per_echo = adc_samples / nr_echoes;

for b = 1:B
    for eco = 1:nr_echoes
        %------------------------------------------------------------------
        % Calculate the index range of samples per echo
        %------------------------------------------------------------------
        echo_range = (1:adc_samples_per_echo) + (eco - 1) * adc_samples_per_echo;

        %------------------------------------------------------------------
        % btraj_unscaled (3 x adc_samples_per_echo x # of TRs per motion state)
        %------------------------------------------------------------------
        cfl_file = fullfile(output_path, sprintf('btraj%d_eco%d_unscaled', b, eco));
        tstart = tic; fprintf('%s:(B=%2d/%2d)(ECO=%d/%d) Writing a .cfl file: %s... ', datetime, b, B, eco, nr_echoes, cfl_file);
        writecfl(cfl_file, traj_unscaled(:,echo_range,bidx{b}));
        fprintf('%s: done! (%6.4f/%6.4f sec)\n', datetime, toc(tstart), toc(start_time));
    end
end

%% Display results
if 1
t = (1:length(resp)).' * TR;
figure('Color', 'w', 'Position', [2 486 1670 492]);
hold on;
plot(t, resp, 'LineWidth', 1);
grid on; grid minor;
set(gca, 'Box', 'On', 'FontSize', 12);
for b = 1:B
    idx = bidx{b};
    plot(t(idx), resp(idx), '.');
end
xlim([0 t(end)]);
xlim([0 200]);
ylim([0.2 1]);
xlabel('Time [sec]', 'FontSize', 12);
ylabel('Amplitude [a.u.]', 'FontSize', 12);
drawnow;
title(sprintf('Respiratory data binning, %d repsiratory states', B), 'Interpreter', 'latex', 'FontSize', 20);
export_fig(fullfile(output_path, sprintf('resp_n%d_fl%4.2f_fh%4.2f_fw%4.2f_B%d', n, fl, fh, fw, B)), '-r300', '-tif');
close gcf;
end
