% demo_bh_prepare_traj.m
% Written by Nam Gyun Lee
% Email: namgyunl@usc.edu, ggang56@gmail.com (preferred)
% Started: 07/27/2022, Last modified: 08/04/2023

%% Clean slate
close all; clearvars -except json_nr json_files json_file; clc;

%% Read a .json file
start_time = tic;
tstart = tic; fprintf('%s: Reading a .json file: %s... \n', datetime, json_file);
fid = fopen(json_file); 
json_txt = fread(fid, [1 inf], 'char=>char'); 
fclose(fid);
json = jsondecode(json_txt);
fprintf('%s: done! (%6.4f/%6.4f sec)\n', datetime, toc(tstart), toc(start_time));

%--------------------------------------------------------------------------
% Define the full path of a filename
%--------------------------------------------------------------------------
if ispc
    siemens_twix_file  = strrep(json.siemens_twix_file, '/', '\');
    ismrmrd_data_file  = strrep(json.ismrmrd_data_file, '/', '\');
    ismrmrd_noise_file = strrep(json.ismrmrd_noise_file, '/', '\');
    trj_file           = strrep(json.trj_file, '/', '\');
    output_path        = strrep(json.output_path, '/', '\');
else
    siemens_twix_file  = json.siemens_twix_file;
    ismrmrd_data_file  = json.ismrmrd_data_file;
    ismrmrd_noise_file = json.ismrmrd_noise_file;
    trj_file           = json.trj_file;
    output_path        = json.output_path;
end

%--------------------------------------------------------------------------
% Define the BART directory
%--------------------------------------------------------------------------
bart_path = json.bart_path;

%--------------------------------------------------------------------------
% Reconstruction parameters
%--------------------------------------------------------------------------
recon_matrix_size   = json.recon_parameters.recon_matrix_size.'; % reconstruction matrix size
dicom_matrix_size   = json.recon_parameters.dicom_matrix_size.'; % DICOM matrix size
recon_interp_factor = json.recon_parameters.recon_interp_factor; % reconstruction interpolation factor: recon_resolution = encoded_resolution / recon_interp_factor;
cal_size            = json.recon_parameters.cal_size.';          % size of the calibration region
lambda              = json.recon_parameters.lambda;              % l1 regularization parameter for pics in BART
max_iter            = json.recon_parameters.max_iter;            % max. number of iterations for pics in BART
kappa               = json.recon_parameters.kappa;               % kappa: kappa=0 => W=I

%--------------------------------------------------------------------------
% Sequence parameters
%--------------------------------------------------------------------------
flag_self         = json.seq.flag_self;
nr_contrasts      = json.seq.nr_contrasts;
nr_interleaves    = json.seq.nr_interleaves;
nr_readouts       = json.seq.nr_readouts;
readout_os_factor = json.seq.readout_os_factor;
resolution        = json.seq.resolution;
TR                = json.seq.TR;
trajectory        = json.seq.trajectory;

%% Make an output directory
mkdir(output_path);

%% Set up a BART command
%--------------------------------------------------------------------------
% Define a BART command
%--------------------------------------------------------------------------
if ispc
    command_prefix = 'wsl';
else
    command_prefix = '';
end
bart_command = sprintf('%s %s/bart', command_prefix, bart_path);

%--------------------------------------------------------------------------
% Translate from a Windows path to a WSL path 
%--------------------------------------------------------------------------
if ispc
    bart_output_path = strrep(output_path, '\', '/');
    bart_output_path = sprintf('/mnt/%s/%s/', lower(bart_output_path(1)), bart_output_path(4:end));
else
    bart_output_path = sprintf('%s/', output_path);
end

%% Read k-space data (ISMRMRD format)
tstart = tic; fprintf('%s: Reading an ISMRMRD file: %s... \n', datetime, ismrmrd_data_file);
if exist(ismrmrd_data_file, 'file')
    dset = ismrmrd.Dataset(ismrmrd_data_file, 'dataset');
    fprintf('%s: done! (%6.4f/%6.4f sec)\n', datetime, toc(tstart), toc(start_time));
else
    error('File %s does not exist.  Please generate it.' , ismrmrd_data_file);
end

%% Get imaging parameters from the XML header
header = ismrmrd.xml.deserialize(dset.readxml);

%--------------------------------------------------------------------------
% measurement information
%--------------------------------------------------------------------------
patient_position = header.measurementInformation.patientPosition;

%--------------------------------------------------------------------------
% Encoding
%--------------------------------------------------------------------------
encoded_fov(1) = header.encoding.encodedSpace.fieldOfView_mm.x * 1e-3; % RO
encoded_fov(2) = header.encoding.encodedSpace.fieldOfView_mm.y * 1e-3; % PE
encoded_fov(3) = header.encoding.encodedSpace.fieldOfView_mm.z * 1e-3; % SL

%% Parse the ISMRMRD header
tstart = tic; fprintf('%s: Parsing the ISMRMRD header... \n', datetime);
raw_data = dset.readAcquisition(1);
fprintf('%s: done! (%6.4f/%6.4f sec)\n', datetime, toc(tstart), toc(start_time));

%--------------------------------------------------------------------------
% ISMRMRD header
%--------------------------------------------------------------------------
% uint16_t version;                                    /**< First unsigned int indicates the version */
% uint64_t flags;                                      /**< bit field with flags */
% uint32_t measurement_uid;                            /**< Unique ID for the measurement */
% uint32_t scan_counter;                               /**< Current acquisition number in the measurement */
% uint32_t acquisition_time_stamp;                     /**< Acquisition clock */
% uint32_t physiology_time_stamp[ISMRMRD_PHYS_STAMPS]; /**< Physiology time stamps, e.g. ecg, breating, etc. */
% uint16_t number_of_samples;                          /**< Number of samples acquired */
% uint16_t available_channels;                         /**< Available coils */
% uint16_t active_channels;                            /**< Active coils on current acquisiton */
% uint64_t channel_mask[ISMRMRD_CHANNEL_MASKS];        /**< Mask to indicate which channels are active. Support for 1024 channels */
% uint16_t discard_pre;                                /**< Samples to be discarded at the beginning of  acquisition */
% uint16_t discard_post;                               /**< Samples to be discarded at the end of acquisition */
% uint16_t center_sample;                              /**< Sample at the center of k-space */
% uint16_t encoding_space_ref;                         /**< Reference to an encoding space, typically only one per acquisition */
% uint16_t trajectory_dimensions;                      /**< Indicates the dimensionality of the trajectory vector (0 means no trajectory) */
% float sample_time_us;                                /**< Time between samples in micro seconds, sampling BW */
% float position[3];                                   /**< Three-dimensional spatial offsets from isocenter */
% float read_dir[3];                                   /**< Directional cosines of the readout/frequency encoding */
% float phase_dir[3];                                  /**< Directional cosines of the phase */
% float slice_dir[3];                                  /**< Directional cosines of the slice direction */
% float patient_table_position[3];                     /**< Patient table off-center */
% ISMRMRD_EncodingCounters idx;                        /**< Encoding loop counters, see above */
% int32_t user_int[ISMRMRD_USER_INTS];                 /**< Free user parameters */
% float user_float[ISMRMRD_USER_FLOATS];               /**< Free user parameters */
%--------------------------------------------------------------------------
% Where EncodingCounters are defined as:
% uint16_t kspace_encode_step_1;    /**< e.g. phase encoding line number */
% uint16_t kspace_encode_step_2;    /**< e.g. partition encoding number */
% uint16_t average;                 /**< e.g. signal average number */
% uint16_t slice;                   /**< e.g. imaging slice number */
% uint16_t contrast;                /**< e.g. echo number in multi-echo */
% uint16_t phase;                   /**< e.g. cardiac phase number */
% uint16_t repetition;              /**< e.g. dynamic number for dynamic scanning */
% uint16_t set;                     /**< e.g. flow encoding set */
% uint16_t segment;                 /**< e.g. segment number for segmented acquisition */
% uint16_t user[ISMRMRD_USER_INTS]; /**< Free user parameters */
%--------------------------------------------------------------------------
number_of_samples  = double(max(raw_data.head.number_of_samples));
discard_pre        = double(max(raw_data.head.discard_pre));
discard_post       = double(max(raw_data.head.discard_post));
center_sample      = double(max(raw_data.head.center_sample));
nr_channels        = double(max(raw_data.head.active_channels));
nr_phase_encodings = double(max(raw_data.head.idx.kspace_encode_step_1)) + 1; % nr_interleaves for spiral imaging
nr_slice_encodings = double(max(raw_data.head.idx.kspace_encode_step_2)) + 1;
nr_averages        = double(max(raw_data.head.idx.average)) + 1;
nr_slices          = double(max(raw_data.head.idx.slice)) + 1;
nr_phases          = double(max(raw_data.head.idx.phase)) + 1;
nr_repetitions     = double(max(raw_data.head.idx.repetition)) + 1;
nr_sets            = double(max(raw_data.head.idx.set)) + 1;
nr_segments        = double(max(raw_data.head.idx.segment)) + 1;

%--------------------------------------------------------------------------
% Get the dimensionality of the trajectory vector (0 means no trajectory)
%--------------------------------------------------------------------------
trajectory_dimensions = double(max(raw_data.head.trajectory_dimensions));

%% Overwrite imaging parameters from ISMRMRD format with parameters from a .seq file
encoded_resolution = ones(1,3) * resolution;

encoded_fov = encoded_fov * readout_os_factor; % [RO,PE,SL] [m]

recon_resolution = encoded_resolution ./ recon_interp_factor;

recon_fov = recon_matrix_size .* recon_resolution; % [m]

%% Display the ISMRMRD header
fprintf('========================= ISMRMRD header =========================\n');
fprintf('encoded_fov        = %8.4f %8.4f %8.4f [mm]\n', encoded_fov(1) * 1e3, encoded_fov(2) * 1e3, encoded_fov(3) * 1e3);
fprintf('encoded_resolution = %8.4f %8.4f %8.4f [mm]\n', encoded_resolution(1) * 1e3, encoded_resolution(2) * 1e3, encoded_resolution(3) * 1e3);
fprintf('------------------------------------------------------------------\n');
fprintf('recon_fov          = %8.4f %8.4f %8.4f [mm]\n', recon_fov(1) * 1e3, recon_fov(2) * 1e3, recon_fov(3) * 1e3);
fprintf('recon_matrix_size  = %d      %d        %d\n', recon_matrix_size(1), recon_matrix_size(2), recon_matrix_size(3));
fprintf('recon_resolution   = %8.4f %8.4f %8.4f [mm]\n', recon_resolution(1) * 1e3, recon_resolution(2) * 1e3, recon_resolution(3) * 1e3);
fprintf('------------------------------------------------------------------\n');
fprintf('trajectory         = %s\n', header.encoding.trajectory);
fprintf('number_of_samples  = %d\n', number_of_samples);
fprintf('discard_pre        = %d\n', discard_pre);
fprintf('discard_post       = %d\n', discard_post);
fprintf('center_sample      = %d\n', center_sample);
fprintf('nr_channels        = %d\n', nr_channels);
fprintf('nr_phase_encodings = %d\n', nr_phase_encodings);
fprintf('nr_slice_encodings = %d\n', nr_slice_encodings);
fprintf('nr_averages        = %d\n', nr_averages);
fprintf('nr_slices          = %d\n', nr_slices);
fprintf('nr_contrasts       = %d\n', nr_contrasts);
fprintf('nr_phases          = %d\n', nr_phases);
fprintf('nr_repetitions     = %d\n', nr_repetitions);
fprintf('nr_sets            = %d\n', nr_sets);
fprintf('nr_segments        = %d\n', nr_segments);
fprintf('==================================================================\n');

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

%% Determine the main orientation of an imaging stack
main_orientation = fGSLClassOri(dNormalSag, dNormalCor, dNormalTra);

%% Get a slice offset in the PCS from Siemens TWIX format
if isfield(twix.hdr.MeasYaps.sSliceArray.asSlice{1}, 'sPosition')
    if isfield(twix.hdr.MeasYaps.sSliceArray.asSlice{1}.sPosition, 'dSag')
        sag_offset = twix.hdr.MeasYaps.sSliceArray.asSlice{1}.sPosition.dSag; % [mm]
    else
        sag_offset = 0; % [mm]
    end
    if isfield(twix.hdr.MeasYaps.sSliceArray.asSlice{1}.sPosition, 'dCor')
        cor_offset = twix.hdr.MeasYaps.sSliceArray.asSlice{1}.sPosition.dCor; % [mm]
    else
        cor_offset = 0; % [mm]
    end
    if isfield(twix.hdr.MeasYaps.sSliceArray.asSlice{1}.sPosition, 'dTra')
        tra_offset = twix.hdr.MeasYaps.sSliceArray.asSlice{1}.sPosition.dTra; % [mm]
    else
        tra_offset = 0; % [mm]
    end
else
    sag_offset = 0; % [mm]
    cor_offset = 0; % [mm]
    tra_offset = 0; % [mm]
end
pcs_offset = [sag_offset; cor_offset; tra_offset] * 1e-3; % [mm] * [m/1e3mm] => [m]

%% Get a slice offset in the PCS from ISMRMRD format
sag_offset_ismrmrd = double(raw_data.head.position(1,1)); % [mm]
cor_offset_ismrmrd = double(raw_data.head.position(2,1)); % [mm]
tra_offset_ismrmrd = double(raw_data.head.position(3,1)); % [mm]
pcs_offset_ismrmrd = [sag_offset_ismrmrd; cor_offset_ismrmrd; tra_offset_ismrmrd] * 1e-3; % [mm] * [m/1e3mm] => [m]

%% Get a rotation matrix from the GCS to the PCS (ISMRMRD format)
phase_dir = double(raw_data.head.phase_dir(:,1));
read_dir  = double(raw_data.head.read_dir(:,1));
slice_dir = double(raw_data.head.slice_dir(:,1));
R_gcs2pcs_ismrmrd = [phase_dir read_dir slice_dir];

%% Calculate a transformation matrix from the RCS to the GCS [r,c,s] <=> [PE,RO,SL]
R_rcs2gcs = [0    1    0 ; % [PE]   [0 1 0] * [r]
             1    0    0 ; % [RO] = [1 0 0] * [c]
             0    0    1]; % [SL]   [0 0 1] * [s]

%% Calculate a rotation matrix from the GCS to the PCS
[R_gcs2pcs,phase_sign,read_sign] = siemens_calculate_transform_gcs_to_pcs(dNormalSag, dNormalCor, dNormalTra, dRotAngle);

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
fprintf('                   [sag]   %10.5f [mm]\n', sag_offset);
fprintf('slice offset     : [cor] = %10.5f [mm]\n', cor_offset);
fprintf('                   [tra]   %10.5f [mm]\n', tra_offset);
fprintf('---------------------- From ISMRMRD format -----------------------\n');
fprintf('                   [sag]   %10.5f [mm]\n', sag_offset_ismrmrd);
fprintf('slice offset     : [cor] = %10.5f [mm]\n', cor_offset_ismrmrd);
fprintf('                   [tra]   %10.5f [mm]\n', tra_offset_ismrmrd);
fprintf('---------------------- From Siemens TWIX format ------------------\n');
fprintf('                   [sag]   [%10.5f %10.5f %10.5f][PE]\n', R_gcs2pcs(1,1), R_gcs2pcs(1,2), R_gcs2pcs(1,3));
fprintf('R_gcs2pcs        : [cor] = [%10.5f %10.5f %10.5f][RO]\n', R_gcs2pcs(2,1), R_gcs2pcs(2,2), R_gcs2pcs(2,3));
fprintf('                   [tra]   [%10.5f %10.5f %10.5f][SL]\n', R_gcs2pcs(3,1), R_gcs2pcs(3,2), R_gcs2pcs(3,3));
fprintf('---------------------- From ISMRMRD format (incorrect!)-----------\n');
fprintf('                   [sag]   [%10.5f %10.5f %10.5f][PE]\n', R_gcs2pcs_ismrmrd(1,1), R_gcs2pcs_ismrmrd(1,2), R_gcs2pcs_ismrmrd(1,3));
fprintf('R_gcs2pcs        : [cor] = [%10.5f %10.5f %10.5f][RO]\n', R_gcs2pcs_ismrmrd(2,1), R_gcs2pcs_ismrmrd(2,2), R_gcs2pcs_ismrmrd(2,3));
fprintf('                   [tra]   [%10.5f %10.5f %10.5f][SL]\n', R_gcs2pcs_ismrmrd(3,1), R_gcs2pcs_ismrmrd(3,2), R_gcs2pcs_ismrmrd(3,3));
fprintf('------------------------------------------------------------------\n');
fprintf('                   [ x ]   [%10.5f %10.5f %10.5f][sag]\n', R_pcs2dcs(1,1), R_pcs2dcs(1,2), R_pcs2dcs(1,3));
fprintf('R_pcs2dcs        : [ y ] = [%10.5f %10.5f %10.5f][cor]\n', R_pcs2dcs(2,1), R_pcs2dcs(2,2), R_pcs2dcs(2,3));
fprintf('                   [ z ]   [%10.5f %10.5f %10.5f][tra]\n', R_pcs2dcs(3,1), R_pcs2dcs(3,2), R_pcs2dcs(3,3));
fprintf('------------------------------------------------------------------\n');
fprintf('                   [ x ]   [%10.5f %10.5f %10.5f][PE]\n', R_gcs2dcs(1,1), R_gcs2dcs(1,2), R_gcs2dcs(1,3));
fprintf('R_gcs2dcs        : [ y ] = [%10.5f %10.5f %10.5f][RO]\n', R_gcs2dcs(2,1), R_gcs2dcs(2,2), R_gcs2dcs(2,3));
fprintf('                   [ z ]   [%10.5f %10.5f %10.5f][SL]\n', R_gcs2dcs(3,1), R_gcs2dcs(3,2), R_gcs2dcs(3,3));
fprintf('==================================================================\n');

%% Define parameters for non-Cartesian reconstruction
%--------------------------------------------------------------------------
% Calculate the number of samples per echo
%--------------------------------------------------------------------------
Nk = (number_of_samples - discard_pre - discard_post);

%--------------------------------------------------------------------------
% Calculate the number of half-spokes per interleave (without SI projections)
%--------------------------------------------------------------------------
Ns = nr_readouts;

%--------------------------------------------------------------------------
% Calculate the number of projections per interleave (spokes + SI projection)
%--------------------------------------------------------------------------
if flag_self
    Np = nr_readouts + 1;
else
    Np = nr_readouts;
end

%--------------------------------------------------------------------------
% Calculate the number of interleaves
%--------------------------------------------------------------------------
Ni = nr_interleaves;

%--------------------------------------------------------------------------
% Calculate the total number of half-spokes (with and without SI projections)
%--------------------------------------------------------------------------
NiNs = Ni * Ns;
NiNp = Ni * Np;

%--------------------------------------------------------------------------
% Calculate the number of channels
%--------------------------------------------------------------------------
Nc = nr_channels;

%--------------------------------------------------------------------------
% Set N1, N2, and N3
%--------------------------------------------------------------------------
N1 = recon_matrix_size(1); % number of samples in the row (r) direction of the RCS
N2 = recon_matrix_size(2); % number of samples in the col (c) direction of the RCS
N3 = recon_matrix_size(3); % number of samples in the slice (s) direction of the RCS
N = N1 * N2 * N3;

%% Calculate a scaling matrix
scaling_matrix = diag(recon_resolution); % [m]

%% Calculate spatial coordinates in the DCS [m]
%--------------------------------------------------------------------------
% Calculate spatial coordinates in the RCS [m]
%--------------------------------------------------------------------------
[I1,I2,I3] = ndgrid((1:N1).', (1:N2).', (1:N3).');
r_rcs = (scaling_matrix * cat(2, I1(:) - (floor(N1/2) + 1), I2(:) - (floor(N2/2) + 1), I3(:) - (floor(N3/2) + 1)).').'; % N x 3

%--------------------------------------------------------------------------
% Calculate spatial coordinates in the DCS [m]
%--------------------------------------------------------------------------
r_dcs = (repmat(dcs_offset, [1 N]) + R_rcs2dcs * r_rcs.').'; % N x 3
x = r_dcs(:,1); % N x 1 [m]
y = r_dcs(:,2); % N x 1 [m]
z = r_dcs(:,3); % N x 1 [m]

%% Read a .trj file
%--------------------------------------------------------------------------
% trj(:,1)  - +X, echo 1 | trj(:,7)  - +X, echo 2
% trj(:,2)  - -X, echo 1 | trj(:,8)  - -X, echo 2
% trj(:,3)  - +Y, echo 1 | trj(:,9)  - +Y, echo 2
% trj(:,4)  - -Y, echo 1 | trj(:,10) - -Y, echo 2
% trj(:,5)  - +Z, echo 1 | trj(:,11) - +Z, echo 2
% trj(:,6)  - -Z, echo 1 | trj(:,12) - -Z, echo 2
%--------------------------------------------------------------------------
[trj_header,trj] = read_trj_file(trj_file); % [cycle/m]

%% Calculate the (phi,theta) parameterization of 3D radial k-space trajectories
%-------------------------------------------------------------------------- 
% Calculate a spiral phyllotaxis pattern (half-spokes, Delacoste 2018 MRM)
%-------------------------------------------------------------------------- 
if strcmp(trajectory, 'Phyllotaxis')
    [phi, theta] = calculate_spiral_phyllotaxis_pattern(nr_readouts, nr_interleaves, flag_self);

%--------------------------------------------------------------------------
% Calculate a wobbling Archimedean spiral pole (WASP) pattern
%-------------------------------------------------------------------------- 
elseif strcmp(trajectory, 'WASP')
    rpa_range = [0 10]; % [degree]
    [phi, theta] = calculate_wasp_pattern(nr_readouts, nr_interleaves, rpa_range);
end

%% Calculate measured k-space trajectories in the DCS [cycle/m] [x,y,z]
k_dcs = zeros(3, Nk, NiNp, 'double');

for i = 1:NiNp
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
    % Calculate measured k-space trajectories in the DCS [cycle/m]
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

%% Calculate measured k-space trajectories in the GCS [cycle/m] [PE,RO,SL]
k_gcs = zeros(3, Nk, NiNp, 'double');
for i = 1:NiNp
    k_gcs(:,:,i) = R_gcs2dcs.' * k_dcs(:,:,i); % 3 x Nk
end

%% Change the sign of k-space trajectories in the GCS [cycle/m] [PE,RO,SL]
k_gcs(1,:,:) = phase_sign * k_gcs(1,:,:); % PE (ku)
k_gcs(2,:,:) = read_sign  * k_gcs(2,:,:); % RO (kv)

%% Calculate measured k-space trajectories in the RCS [cycle/m] [r,c,s]
k_rcs = zeros(3, Nk, NiNp, 'double');
for i = 1:NiNp
    k_rcs(:,:,i) = R_rcs2gcs.' * k_gcs(:,:,i); % 3 x Nk
end

%% Calculate measured k-space trajectories only for imaging data in the RCS [cycle/m] [r,c,s]
if flag_self
    start_index = 1;
else
    start_index = 0;
end

imaging_data_indices = zeros(NiNs, 1, 'double');
for idx = 1:Ni
    imaging_data_indices((1:Ns).' + (idx - 1) * Ns) = start_index + (1:Ns).' + (idx - 1) * Np;
end
k_rcs_imaging_data = k_rcs(:,:,imaging_data_indices); % 3 x Nk x NiNs

%% Calculate "unscaled" k-space trajectories for BART (3 x Nk x NiNs)
%--------------------------------------------------------------------------
% A minus is required to perform forward FFT to move from k-space to image space
% BART requires "traj" scaled to [-0.5,0.5] * N
% Here, we calculate "unscaled" traj scaled to [-0.5,0.5]
%--------------------------------------------------------------------------
krmax = 1 / (2 * resolution); % [cycle/m]
traj_unscaled = zeros(3, Nk, NiNs, 'single');
traj_unscaled(1,:,:) = -k_rcs_imaging_data(1,:,:) / (2 * krmax);
traj_unscaled(2,:,:) = -k_rcs_imaging_data(2,:,:) / (2 * krmax); % AP direction
traj_unscaled(3,:,:) = -k_rcs_imaging_data(3,:,:) / (2 * krmax);

%% Write a .cfl file
for echo = 1:nr_contrasts
    echo_range = (1:Nk/nr_contrasts) + (echo - 1) * Nk / nr_contrasts;
    %----------------------------------------------------------------------
    % btraj_unscaled (3 x Nk/nr_contrasts x NiNs)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path, sprintf('traj_echo%d_unscaled', echo));
    tstart = tic; fprintf('%s:(e=%d/%d) Writing a .cfl file: %s... \n', datetime, echo, nr_contrasts, cfl_file);
    writecfl(cfl_file, traj_unscaled(:,echo_range,:));
    fprintf('%s: done! (%6.4f/%6.4f sec)\n', datetime, toc(tstart), toc(start_time));
end
