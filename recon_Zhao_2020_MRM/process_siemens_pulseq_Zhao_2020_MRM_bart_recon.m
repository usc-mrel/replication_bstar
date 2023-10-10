function [trj_header, trj, trj_filename, read_direction, phase_direction] = process_siemens_pulseq_Zhao_2020_MRM_bart_recon(siemens_twix_path, ismrmrd_data_path, ismrmrd_noise_path, seq_scope_path, seq_bstar_path)
% Written by Nam Gyun Lee
% Email: namgyunl@usc.edu, ggang56@gmail.com (preferred)
% Started: 10/02/2022, Last modified: 10/19/2022


%% Read a .seq file
%--------------------------------------------------------------------------
% scope
%--------------------------------------------------------------------------
seq_scope = mr.Sequence;
seq_scope.read(seq_scope_path);

%--------------------------------------------------------------------------
% bstar
%--------------------------------------------------------------------------
seq_bstar = mr.Sequence;
seq_bstar.read(seq_bstar_path);

%% Get sequence parameters from .seq file
%--------------------------------------------------------------------------
% scope
%--------------------------------------------------------------------------
adc_grad_delay            = seq_scope.getDefinition('ADCGradDelay');        % ADC-Grad Delay [sec]
amplitude_scope           = seq_scope.getDefinition('Amplitude');           % Amplitude [mT/m]
nr_echoes                 = seq_scope.getDefinition('Contrasts');           % Contrasts
fov                       = seq_scope.getDefinition('FOV');                 % FOV [m]
max_slewrate              = seq_scope.getDefinition('MaxSlewRate');         % Max SlewRate [mT/m/ms]
nr_phase_encoding_steps_1 = seq_scope.getDefinition('PhaseEncodingSteps1');
nr_phase_encoding_steps_2 = seq_scope.getDefinition('PhaseEncodingSteps2');
rf_duration               = seq_scope.getDefinition('RFDuration');          % RF Duration [sec]
ramp_time                 = seq_scope.getDefinition('RampTime');            % RampTime [sec]
readout_os_factor         = seq_scope.getDefinition('ReadoutOSFactor');
real_dwell_time_scope     = seq_scope.getDefinition('RealDwellTime');       % RealDwellTime [sec]
total_time                = seq_scope.getDefinition('TotalTime');           % TotalTime [sec]

%--------------------------------------------------------------------------
% bstar
%--------------------------------------------------------------------------
adc_raster_time           = seq_bstar.getDefinition('AdcRasterTime');  % AdcRasterTime [sec]
amplitude_bstar           = seq_bstar.getDefinition('Amplitude');      % Amplitude [mT/m/ms]
base_resolution           = seq_bstar.getDefinition('BaseResolution'); % Base Resolution
real_dwell_time_bstar1    = seq_bstar.getDefinition('RealDwellTime1'); % RealDwellTime1 [sec]
real_dwell_time_bstar2    = seq_bstar.getDefinition('RealDwellTime2'); % RealDwellTime2 [sec]
resolution                = seq_bstar.getDefinition('Resolution');     % Resolution [m]
real_dwell_time_bstar     = real_dwell_time_bstar1;

fprintf('======================scope special card ========================\n');
fprintf('RampTime       = %4.0f [us]      RF Duration     = %4.0f [us]\n', ramp_time * 1e6, rf_duration * 1e6);
fprintf('TotalTime      = %4.0f [us]      Slice Thickness = %4.2f [mm]\n', total_time * 1e6, 0.005 * 1e3);
fprintf('ADC-Grad Delay = %4.0f [us]      Slice Shift     = %4.0f [mm]\n', adc_grad_delay * 1e6, 0 * 1e3);
fprintf('Amplitude      = %4.1f [mT/m]\n', amplitude_scope);
fprintf('[\bDwell Time     = %4.0f [ns]]\b\n', real_dwell_time_scope * 1e9);
fprintf('Max SlewRate   = %4.0f [mT/m/ms]\n', max_slewrate);
fprintf('-----------------------------------------------------------------\n');
fprintf('scope: real dwell time          = %7.2f [ns]\n', real_dwell_time_scope * 1e9);
fprintf('bSTAR: real dwell time (echo 1) = %7.2f [ns]\n', real_dwell_time_bstar1 * 1e9);
fprintf('bSTAR: real dwell time (echo 2) = %7.2f [ns]\n', real_dwell_time_bstar2 * 1e9);
fprintf('=================================================================\n');

%% Calculate the bandwidth [Hz/Px]
bandwidth_scope  = 1 / real_dwell_time_scope  / (2 * base_resolution); % [Hz/Px]
bandwidth_bstar1 = 1 / real_dwell_time_bstar1 / (2 * base_resolution); % [Hz/Px]
bandwidth_bstar2 = 1 / real_dwell_time_bstar2 / (2 * base_resolution); % [Hz/Px]
bandwidth = bandwidth_bstar1;
fprintf('Bandwidth (scope)  = %4.0f [Hz/Px]\n', bandwidth_scope);
fprintf('Bandwidth (echo 1) = %4.0f [Hz/Px]\n', bandwidth_bstar1);
fprintf('Bandwidth (echo 2) = %4.0f [Hz/Px]\n', bandwidth_bstar2);

%% Read k-space data (ISMRMRD format)
start_time = tic;
tic; fprintf('Reading an ISMRMRD file: %s... ', ismrmrd_data_path);
if exist(ismrmrd_data_path, 'file')
    dset = ismrmrd.Dataset(ismrmrd_data_path, 'dataset');
    fprintf('done! (%6.4f/%6.4f sec)\n', toc, toc(start_time));
else
    error('File %s does not exist.  Please generate it.' , ismrmrd_data_path);
end

%% Get imaging parameters from the XML header
header = ismrmrd.xml.deserialize(dset.readxml);

%--------------------------------------------------------------------------
% measurement information
%--------------------------------------------------------------------------
patient_position = header.measurementInformation.patientPosition;

%--------------------------------------------------------------------------
% Sequence parameters
%--------------------------------------------------------------------------
TR         = header.sequenceParameters.TR * 1e-3;     % [msec] * [sec/1e3msec] => [sec]
TE         = header.sequenceParameters.TE * 1e-3;     % [msec] * [sec/1e3msec] => [sec]
flip_angle = header.sequenceParameters.flipAngle_deg; % [degrees]

%--------------------------------------------------------------------------
% Encoding
%--------------------------------------------------------------------------
encoded_fov(1) = header.encoding.encodedSpace.fieldOfView_mm.x; % RO
encoded_fov(2) = header.encoding.encodedSpace.fieldOfView_mm.y; % PE
encoded_fov(3) = header.encoding.encodedSpace.fieldOfView_mm.z; % SL

recon_fov(1) = header.encoding.reconSpace.fieldOfView_mm.x; % RO
recon_fov(2) = header.encoding.reconSpace.fieldOfView_mm.y; % PE
recon_fov(3) = header.encoding.reconSpace.fieldOfView_mm.z; % SL

Nkx = header.encoding.encodedSpace.matrixSize.x; % number of readout samples in k-space
Nky = header.encoding.encodedSpace.matrixSize.y; % number of phase encodes in k-space
Nkz = header.encoding.encodedSpace.matrixSize.z; % number of slice encodes in k-space
Nx  = header.encoding.reconSpace.matrixSize.x;   % number of samples in image-space (RO)
Ny  = header.encoding.reconSpace.matrixSize.y;   % number of samples in image-space (PE)
Nz  = header.encoding.reconSpace.matrixSize.z;   % number of samples in image-space (SL)
Nc  = header.acquisitionSystemInformation.receiverChannels;

encoded_resolution = encoded_fov ./ [Nkx Nky Nkz]; % [mm]
recon_resolution = recon_fov ./ [Nx Ny Nz]; % [mm]

%% Parse the ISMRMRD header
tic; fprintf('Parsing the ISMRMRD header... ');
raw_data = dset.readAcquisition(); % read all the acquisitions
fprintf('done! (%6.4f/%6.4f sec)\n', toc, toc(start_time));

%--------------------------------------------------------------------------
% Get data type
%--------------------------------------------------------------------------
is_noise                     = raw_data.head.flagIsSet('ACQ_IS_NOISE_MEASUREMENT');
is_parallel_calibration      = raw_data.head.flagIsSet('ACQ_IS_PARALLEL_CALIBRATION');
is_navigation                = raw_data.head.flagIsSet('ACQ_IS_NAVIGATION_DATA');
is_phasecorr                 = raw_data.head.flagIsSet('ACQ_IS_PHASECORR_DATA');
is_hpfeedback                = raw_data.head.flagIsSet('ACQ_IS_HPFEEDBACK_DATA');
is_dummyscan                 = raw_data.head.flagIsSet('ACQ_IS_DUMMYSCAN_DATA');
is_rtfeedback                = raw_data.head.flagIsSet('ACQ_IS_RTFEEDBACK_DATA');
is_surfacecoilcorrectionscan = raw_data.head.flagIsSet('ACQ_IS_SURFACECOILCORRECTIONSCAN_DATA');

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
number_of_samples  = double(max(raw_data.head.number_of_samples(~is_noise)));
discard_pre        = double(max(raw_data.head.discard_pre));
discard_post       = double(max(raw_data.head.discard_post));
center_sample      = double(max(raw_data.head.center_sample));
nr_channels        = double(max(raw_data.head.active_channels));
nr_phase_encodings = double(max(raw_data.head.idx.kspace_encode_step_1)) + 1; % nr_interleaves for spiral imaging
nr_slice_encodings = double(max(raw_data.head.idx.kspace_encode_step_2)) + 1;
nr_averages        = double(max(raw_data.head.idx.average)) + 1;
nr_slices          = double(max(raw_data.head.idx.slice)) + 1;
nr_contrasts       = double(max(raw_data.head.idx.contrast)) + 1;
nr_phases          = double(max(raw_data.head.idx.phase)) + 1;
nr_repetitions     = double(max(raw_data.head.idx.repetition)) + 1;
nr_sets            = double(max(raw_data.head.idx.set)) + 1;
nr_segments        = double(max(raw_data.head.idx.segment)) + 1;

%--------------------------------------------------------------------------
% Get the dimensionality of the trajectory vector (0 means no trajectory)
%--------------------------------------------------------------------------
trajectory_dimensions = double(max(raw_data.head.trajectory_dimensions));

%--------------------------------------------------------------------------
% Get the dwell time in [sec]
%--------------------------------------------------------------------------
dt = double(max(raw_data.head.sample_time_us)) * 1e-6; % [usec] * [sec/1e-6 usec] => [sec]

%--------------------------------------------------------------------------
% Calculate the readout duration [sec]
%--------------------------------------------------------------------------
T = (number_of_samples - discard_pre - discard_post) * dt; % readout duration [sec]

%% Display ISMRMRD header
fprintf('========================= ISMRMRD header =========================\n');
fprintf('encoded_fov        = %8.4f %8.4f %8.4f\n', encoded_fov(1), encoded_fov(2), encoded_fov(3));
fprintf('Nkx Nky Nkz        = %d      %d        %d\n', Nkx, Nky, Nkz);
fprintf('encoded_resolution = %8.4f %8.4f %8.4f\n', encoded_resolution(1), encoded_resolution(2), encoded_resolution(3));
fprintf('------------------------------------------------------------------\n');
fprintf('recon_fov          = %8.4f %8.4f %8.4f\n', recon_fov(1), recon_fov(2), recon_fov(3));
fprintf('Nx Ny Nz           = %d      %d        %d\n', Nx, Ny, Nz);
fprintf('recon_resolution   = %8.4f %8.4f %8.4f\n', recon_resolution(1), recon_resolution(2), recon_resolution(3));
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
fprintf('dt                 = %5.2f [usec]\n', dt * 1e6);
fprintf('readout duration   = %5.2f [msec]\n', T * 1e3);
fprintf('==================================================================\n');

%% Calculate the receiver noise matrix
[Psi,inv_L] = calculate_receiver_noise_matrix(ismrmrd_noise_path);

%% Read a Siemens .dat file
fprintf('Reading a Siemens .dat file: %s\n', siemens_twix_path);
twix = mapVBVD(siemens_twix_path);
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

%% Calculate a transformation matrix from the RCS to the GCS [r,c,s] <=> [PE,RO,SL]
R_rcs2gcs = [0    1    0 ; % [PE]   [0 1 0] * [r]
             1    0    0 ; % [RO] = [1 0 0] * [c]
             0    0    1]; % [SL]   [0 0 1] * [s]

%% Calculate a rotation matrix from the GCS to the PCS
[R_gcs2pcs,phase_sign,read_sign] = siemens_calculate_matrix_gcs_to_pcs(dNormalSag, dNormalCor, dNormalTra, dRotAngle);

%% Calculate a rotation matrix from the PCS to the DCS
R_pcs2dcs = siemens_calculate_matrix_pcs_to_dcs(patient_position);

%% Calculate a rotation matrix from the GCS to the DCS
R_gcs2dcs = R_pcs2dcs * R_gcs2pcs;

%% Calculate a rotation matrix from the RCS to the DCS
R_rcs2dcs = R_pcs2dcs * R_gcs2pcs * R_rcs2gcs;

%% Calculate a slice offset in the DCS [m]
dcs_offset = R_pcs2dcs * pcs_offset; % 3 x 1

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

%% Calculate a rotation matrix from the GCS to the DCS (ISMRMRD format)
R_gcs2dcs_ismrmrd = R_pcs2dcs * R_gcs2pcs_ismrmrd;

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
fprintf('---------------------- From ISMRMRD format ------------------------\n');
fprintf('                   [sag]   [%10.5f %10.5f %10.5f][PE]\n', R_gcs2pcs_ismrmrd(1,1), R_gcs2pcs_ismrmrd(1,2), R_gcs2pcs_ismrmrd(1,3));
fprintf('R_gcs2pcs_ismrmrd: [cor] = [%10.5f %10.5f %10.5f][RO]\n', R_gcs2pcs_ismrmrd(2,1), R_gcs2pcs_ismrmrd(2,2), R_gcs2pcs_ismrmrd(2,3));
fprintf('                   [tra]   [%10.5f %10.5f %10.5f][SL]\n', R_gcs2pcs_ismrmrd(3,1), R_gcs2pcs_ismrmrd(3,2), R_gcs2pcs_ismrmrd(3,3));
fprintf('------------------------------------------------------------------\n');
fprintf('                   [ x ]   [%10.5f %10.5f %10.5f][sag]\n', R_pcs2dcs(1,1), R_pcs2dcs(1,2), R_pcs2dcs(1,3));
fprintf('R_pcs2dcs        : [ y ] = [%10.5f %10.5f %10.5f][cor]\n', R_pcs2dcs(2,1), R_pcs2dcs(2,2), R_pcs2dcs(2,3));
fprintf('                   [ z ]   [%10.5f %10.5f %10.5f][tra]\n', R_pcs2dcs(3,1), R_pcs2dcs(3,2), R_pcs2dcs(3,3));
fprintf('---------------------- From Siemens TWIX format ------------------\n');
fprintf('                   [ x ]   [%10.5f %10.5f %10.5f][PE]\n', R_gcs2dcs(1,1), R_gcs2dcs(1,2), R_gcs2dcs(1,3));
fprintf('R_gcs2dcs        : [ y ] = [%10.5f %10.5f %10.5f][RO]\n', R_gcs2dcs(2,1), R_gcs2dcs(2,2), R_gcs2dcs(2,3));
fprintf('                   [ z ]   [%10.5f %10.5f %10.5f][SL]\n', R_gcs2dcs(3,1), R_gcs2dcs(3,2), R_gcs2dcs(3,3));
fprintf('---------------------- From ISMRMRD format ------------------------\n');
fprintf('                   [ x ]   [%10.5f %10.5f %10.5f][PE]\n', R_gcs2dcs_ismrmrd(1,1), R_gcs2dcs_ismrmrd(1,2), R_gcs2dcs_ismrmrd(1,3));
fprintf('R_gcs2dcs_ismrmrd: [ y ] = [%10.5f %10.5f %10.5f][RO]\n', R_gcs2dcs_ismrmrd(2,1), R_gcs2dcs_ismrmrd(2,2), R_gcs2dcs_ismrmrd(2,3));
fprintf('                   [ z ]   [%10.5f %10.5f %10.5f][SL]\n', R_gcs2dcs_ismrmrd(3,1), R_gcs2dcs_ismrmrd(3,2), R_gcs2dcs_ismrmrd(3,3));
fprintf('------------------------------------------------------------------\n');
fprintf('Use rotation matrices from TWIX format to calculate spatial coordinates\n');
fprintf('Use "R_gcs2dcs_ismrmrd" to determine the sign of a gradient\n');
fprintf('==================================================================\n');

%% Define reconstruction parameters
%--------------------------------------------------------------------------
% Set the number of phase-encoding steps in the RO direction
%--------------------------------------------------------------------------
Nkx = nr_phase_encoding_steps_1;

%--------------------------------------------------------------------------
% Set the number of phase-encoding steps in the PE direction
%--------------------------------------------------------------------------
Nky = nr_phase_encoding_steps_2;

%--------------------------------------------------------------------------
% Set the number of samples
%--------------------------------------------------------------------------
Nt = (number_of_samples - discard_pre - discard_post);

%--------------------------------------------------------------------------
% Set the number of samples in image-space (RO)
%--------------------------------------------------------------------------
Nx = Nkx;

%--------------------------------------------------------------------------
% Set the number of samples in image-space (PE)
%--------------------------------------------------------------------------
Ny = Nky;

%--------------------------------------------------------------------------
% Calculate the encoded resolution [m]
%--------------------------------------------------------------------------
encoded_resolution = fov ./ [Nkx Nky 1].'; % [mm]

%--------------------------------------------------------------------------
% Calculate the reconstruction matrix size
%--------------------------------------------------------------------------
recon_matrix_size = [Nkx Nky 1].';
N1 = recon_matrix_size(1); % number of samples in the row (r) direction of the RCS
N2 = recon_matrix_size(2); % number of samples in the col (c) direction of the RCS
N3 = recon_matrix_size(3); % number of samples in the slice (s) direction of the RCS
N = N1 * N2 * N3;

%--------------------------------------------------------------------------
% Calculate the recon resolution
%--------------------------------------------------------------------------
recon_resolution = fov ./ recon_matrix_size; % [mm]

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

%% Update ISMRMRD parameters
header.encoding.encodingLimits.kspace_encoding_step_1.center  = Nkx / 2;
header.encoding.encodingLimits.kspace_encoding_step_1.maximum = Nkx - 1;
header.encoding.encodingLimits.kspace_encoding_step_2.center  = Nky / 2;
header.encoding.encodingLimits.kspace_encoding_step_2.maximum = Nky - 1;

%% Get information about k-space encoding (counting from 0)
kspace_encoding_step_1_center = header.encoding.encodingLimits.kspace_encoding_step_1.center;
kspace_encoding_step_2_center = header.encoding.encodingLimits.kspace_encoding_step_2.center;
kspace_encoding_step_1_maximum = header.encoding.encodingLimits.kspace_encoding_step_1.maximum;
kspace_encoding_step_2_maximum = header.encoding.encodingLimits.kspace_encoding_step_2.maximum;

%% Get k-space data (Nkx x Nky x 1 x Nc x ones(1,6) x Nt x nr_repetitions)
nr_encodings = Nkx * Nky;
ksp = complex(zeros([Nkx Nky 1 Nc ones(1,6) Nt nr_repetitions], 'single'));

for idx2 = 1:nr_repetitions
    offset = (idx2 - 1) * nr_encodings;
    for idx1 = 1:nr_encodings
        tstart = tic; fprintf('(%d/%d): Reading k-space imaging data (%d/%d)... ', idx2, nr_repetitions, idx1, nr_encodings);
        %------------------------------------------------------------------
        % Calculate the current index
        %------------------------------------------------------------------
        index = idx1 + offset;

        %------------------------------------------------------------------
        % Calculate the (1st,2nd) matrix index of a profile
        %------------------------------------------------------------------
        kspace_encode_step_1 = double(raw_data.head.idx.kspace_encode_step_1(index));
        kspace_encode_step_2 = double(raw_data.head.idx.kspace_encode_step_2(index));
        k1_index = kspace_encode_step_1 - (kspace_encoding_step_1_maximum - kspace_encoding_step_1_center + 1) + floor(Nkx/2) + 1;
        k2_index = kspace_encode_step_2 - (kspace_encoding_step_2_maximum - kspace_encoding_step_2_center + 1) + floor(Nky/2) + 1;

        %------------------------------------------------------------------
        % Prewhiten k-space data
        %------------------------------------------------------------------
        profile = raw_data.data{index}; % number_of_samples x nr_channels
        profile = (inv_L * profile.').';

        %------------------------------------------------------------------
        % Accumulate k-space
        %------------------------------------------------------------------
        ksp(k1_index,k2_index,1,:,1,1,1,1,1,1,:,idx2) = reshape(profile.', [1 1 1 Nc ones(1,6) Nt]);
        fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
    end
end

%% Flip k-space
tstart = tic; fprintf('Flipping k-space... ');
if read_sign == -1
    ksp = flip(ksp,1);
end
if phase_sign == -1
    ksp = flip(ksp,2);
end
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%% Perform forward FFT on k-space twice to prepare cailbration data (Nkx x Nky x 1 x Nc)
%--------------------------------------------------------------------------
% Siemens: k-space ("ksp") <=> image-space
% BART:                        image-space <=> k-space ("ksp_calib")
%--------------------------------------------------------------------------
ksp_calib = ksp(:,:,:,:,1,1,1,1,1,1,1);
tstart = tic; fprintf('BART: performing FFT on calibration data... ');
ksp_calib = bart('fft -u 3', ksp_calib);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

tstart = tic; fprintf('BART: performing FFT on calibration data... ');
ksp_calib = bart('fft -u 3', ksp_calib);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%% ESPIRiT calibration
cal_size = [6 6 1]; % size of calibration region

%--------------------------------------------------------------------------
% ESPIRiT calibration (-I: with intensity correction)
%--------------------------------------------------------------------------
tstart = tic; fprintf('BART: ESPIRiT calibration... ');
maps = bart(sprintf('ecalib -t 0.02 -c 0 -r%d:%d:%d -m1', cal_size(1), cal_size(2), cal_size(3)), ksp_calib);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%% Perform FFT reconstruction (Nkx x Nky x 1 x Nc x ones(1,6) x Nt x nr_repetitions)
%--------------------------------------------------------------------------
% Siemens: k-space <=> image-space
%--------------------------------------------------------------------------
tstart = tic; fprintf('BART: performing FFT on "ksp"... ');
imc = bart('fft -u 3', ksp);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%% Perform coil combination
im = reshape(sum(bsxfun(@times, conj(maps), imc), 4), [N1 N2 Nt nr_repetitions]);

%% Calculate normalized phase images [rad]
im_phase_normalized = angle(im ./ im(:,:,1)); % [rad]

%% Calculate a mask
im_ = sum(sum(abs(im),3),4);
mask = (im_ > 10 * min(im_(:)));
mask_indices = find(mask > 0);

nan_mask = zeros(N1, N2, 'single');
nan_mask(mask) = 1;
nan_mask(~mask) = NaN;

%% Perform 2D phase unwrapping
im_phase_unwrapped = zeros(N1, N2, Nt, nr_repetitions, 'single');
for idx2 = 1:nr_repetitions
    for idx1 = 1:Nt
        tstart = tic; fprintf('(%d/%d): Performing 2D phase unwrapping (%d/%d)... ', idx2, nr_repetitions, idx1, Nt);
        im_phase_unwrapped(:,:,idx1,idx2) = unwrap_phase(nan_mask .* im_phase_normalized(:,:,idx1,idx2));
        fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
    end
end

%% Calculate the readout direction in the DCS
%--------------------------------------------------------------------------
% TRANSVERSAL, Phase enc. dir. = A >> P
% phase_sign = +1, read_sign = -1
%                    [ x ]   [   0.00000   -1.00000    0.00000][PE]
% R_gcs2dcs_ismrmrd: [ y ] = [  -1.00000   -0.00000    0.00000][RO]
%                    [ z ]   [   0.00000    0.00000   -1.00000][SL]
%--------------------------------------------------------------------------
[~,index] = max(abs(R_gcs2dcs_ismrmrd(:,2)));

if index == 1
    read_direction = 'x';
elseif index == 2
    read_direction = 'y';
elseif index == 3
    read_direction = 'z';
end

grad_sign = sign(R_gcs2dcs_ismrmrd(index,2));

%% Calculate the phase encoding direction in the DCS
[~,index] = max(abs(R_gcs2dcs_ismrmrd(:,1)));

if index == 1
    phase_direction = 'x';
elseif index == 2
    phase_direction = 'y';
elseif index == 3
    phase_direction = 'z';
end

%% Perform linear fitting of the phase
%--------------------------------------------------------------------------
% x in [m], kx in [rad/m]
% A * c = [x 1] [  kx  ] = rhs
%               [offset]
% OLS: A.' * A * c = A.' * rhs
%      c_hat = inv(A.' * A) * (A.' * rhs)
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Calculate A
%--------------------------------------------------------------------------
A = zeros(length(mask_indices), 2, 'single');
A(:,2) = 1;

if strcmp(read_direction, 'x')
    A(:,1) = x(mask_indices);
elseif strcmp(read_direction, 'y')
    A(:,1) = y(mask_indices);
elseif strcmp(read_direction, 'z')
    A(:,1) = z(mask_indices);
end

%--------------------------------------------------------------------------
% Perform linear fitting of the phase per time point
%--------------------------------------------------------------------------
k = zeros(Nt, nr_repetitions, 'single'); % [rad/m]
for idx2 = 1:nr_repetitions
    for idx1 = 1:Nt
        tstart = tic; fprintf('(%d/%d): Performing linear fitting of the phase (%d/%d)... ', idx2, nr_repetitions, idx1, Nt);
        rhs = reshape(im_phase_unwrapped(:,:,idx1,idx2), [N 1]);
        rhs = rhs(mask_indices); % [rad]
        c_hat = inv(A.' * A) * (A.' * rhs); % [rad/m]
        k(idx1,idx2) = c_hat(1) / (2 * pi); % [cycle/m]
        fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
    end
end

%% Calculate the measured k-space trajectories [cycle/m] (p: plus, m: minus)
if grad_sign == 1
    k_p_scope = k(:,1);
    k_n_scope = k(:,2);
else    
    k_p_scope = k(:,2);
    k_n_scope = k(:,1);
end

%% Caclulate the measured gradient waveforms [mT/m]
gamma = 42576000; % [Hz/T]

% 1 / ([Hz/T]) * [cycle/m] / [sec] * [1e3mT/T] => [mT/m]
g_p_scope = 1 / gamma * cat(1, zeros(1,1), diff(k_p_scope,1,1)) / real_dwell_time_scope * 1e3; % [mT/m]
g_n_scope = 1 / gamma * cat(1, zeros(1,1), diff(k_n_scope,1,1)) / real_dwell_time_scope * 1e3; % [mT/m]

%% Calculate the number of samples for ADC in bSTAR
adc_samples_bstar = nr_echoes * base_resolution * readout_os_factor;
adc_duration_bstar = adc_samples_bstar * real_dwell_time_bstar; % [sec]
fprintf('RUTime/TotalTime/ADC[0]/Resolution: %3.0f/%3.0f/%6.1f/%4.2f\n', ramp_time * 1e6, total_time * 1e6, adc_duration_bstar * 1e6, resolution * 1e3);

%% Calculate a time vector for scope [sec]
t_scope = ((0:Nt-1).' + 0.5) * real_dwell_time_scope; % [sec]

%% Calculate a time vector for bipolar k-space trajectories in bSTAR [sec]
% NOT WORKING?? A BUG?? Here, adc_delay must be a multiple of 1 us (=1000 ns) instead of 100 ns (ADC RASTER TIME). 
shift_adc_bstar = round((2 * total_time - adc_duration_bstar) / 2 / (adc_raster_time * 10)) * (adc_raster_time * 10); % [sec]
t_bstar1 = adc_grad_delay + shift_adc_bstar + ((0:adc_samples_bstar-1).' + 0.5) * real_dwell_time_bstar; % [sec]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5555
% Manual echo2 tuning!!!
t_bstar2 = t_bstar1 + real_dwell_time_bstar / 2; % [sec]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5555

%% Interpolate biploar k-space trajectories [cycle/m]
k_p_bstar_echo1 = interp1(t_scope, k_p_scope, t_bstar1, 'spline');
k_n_bstar_echo1 = interp1(t_scope, k_n_scope, t_bstar1, 'spline');

k_p_bstar_echo2 = interp1(t_scope, k_p_scope, t_bstar2, 'spline');
k_n_bstar_echo2 = interp1(t_scope, k_n_scope, t_bstar2, 'spline');

% Original code without tuning
% k_p_bstar = interp1(t_scope, k_p_scope, t_bstar, 'spline');
% k_n_bstar = interp1(t_scope, k_n_scope, t_bstar, 'spline');

%% Calculate scaled k-space trajectories per echo [cycle/m]
nr_dirs    = nr_repetitions;   % number of directions, e.g.) +X, -X, +Y, -Y, +Z and -Z physical axes
nr_samples = base_resolution * readout_os_factor; % number of samples per direction

scale_factor = amplitude_bstar / amplitude_scope;

echo1_range = (1:nr_samples).';
echo2_range = (1:nr_samples).' + nr_samples;

trj = zeros(nr_samples, 6 * nr_echoes, 'double');

if strcmp(read_direction, 'x')
    trj(:,1)  =  k_p_bstar_echo1(echo1_range) * scale_factor;
    trj(:,2)  = -k_n_bstar_echo1(echo1_range) * scale_factor;
    trj(:,7)  =  k_p_bstar_echo2(echo2_range) * scale_factor;
    trj(:,8)  = -k_n_bstar_echo2(echo2_range) * scale_factor;
elseif strcmp(read_direction, 'y')
    trj(:,3)  =  k_p_bstar_echo1(echo1_range) * scale_factor;
    trj(:,4)  = -k_n_bstar_echo1(echo1_range) * scale_factor;
    trj(:,9)  =  k_p_bstar_echo2(echo2_range) * scale_factor;
    trj(:,10) = -k_n_bstar_echo2(echo2_range) * scale_factor;
elseif strcmp(read_direction, 'z')
    trj(:,5)  =  k_p_bstar_echo1(echo1_range) * scale_factor;
    trj(:,6)  = -k_n_bstar_echo1(echo1_range) * scale_factor;
    trj(:,11) =  k_p_bstar_echo2(echo2_range) * scale_factor;
    trj(:,12) = -k_n_bstar_echo2(echo2_range) * scale_factor;
end

%% Create a .trj file
trj_header.id          = zeros(4,1);
trj_header.lNumEchoes  = nr_echoes;    % number of echoes
trj_header.lNumDirs    = nr_dirs;      % number of directions
trj_header.lNumSamples = nr_samples;   % number of samples per direction
trj_filename = sprintf('traj_b%4.0f_n%d_Zhao_2020_MRM_v2_fine_tuning.trj', bandwidth, nr_samples);
end