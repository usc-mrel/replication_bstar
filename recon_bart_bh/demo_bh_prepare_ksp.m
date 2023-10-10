% demo_bh_prepare_ksp.m
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
    ismrmrd_data_file  = strrep(json.ismrmrd_data_file, '/', '\');
    ismrmrd_noise_file = strrep(json.ismrmrd_noise_file, '/', '\');
    output_path        = strrep(json.output_path, '/', '\');
else
    ismrmrd_data_file  = json.ismrmrd_data_file;
    ismrmrd_noise_file = json.ismrmrd_noise_file;
    output_path        = json.output_path;
end

%--------------------------------------------------------------------------
% Sequence parameters
%--------------------------------------------------------------------------
flag_self      = json.seq.flag_self;
nr_contrasts   = json.seq.nr_contrasts;
nr_interleaves = json.seq.nr_interleaves;
nr_readouts    = json.seq.nr_readouts;

%--------------------------------------------------------------------------
% Coil selection
%--------------------------------------------------------------------------
if isfield(json, 'coil_selection')
    coil_selection = json.coil_selection;
else
    coil_selection = [];
end

%% Make an output directory
mkdir(output_path);

%% Read k-space data (ISMRMRD format)
tstart = tic; fprintf('%s: Reading an ISMRMRD file: %s... \n', datetime, ismrmrd_data_file);
if exist(ismrmrd_data_file, 'file')
    dset = ismrmrd.Dataset(ismrmrd_data_file, 'dataset');
    fprintf('%s: done! (%6.4f/%6.4f sec)\n', datetime, toc(tstart), toc(start_time));
else
    error('File %s does not exist.  Please generate it.' , ismrmrd_data_file);
end

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

%% Calculate the receiver noise matrix
[Psi,inv_L] = calculate_receiver_noise_matrix(ismrmrd_noise_file);

%% Get k-space imaging data (1 x Nk x NiNs x Nc)
%--------------------------------------------------------------------------
% ( 1)not used   x ( 2)readout dimension x ( 3)number of TRs  x ( 4)COIL_DIM
% ( 5)MAPS_DIM   x ( 6)TE_DIM            x ( 7)COEFF_DIM      x ( 8)COEFF2_DIM
% ( 9)ITER_DIM   x (10)CSHIFT_DIM        x (11)TIME_DIM       x (12)TIME2_DIM
% (13)LEVEL_DIM  x (14)SLICE_DIM         x (15)AVG_DIM        x (16)BATCH_DIM
%--------------------------------------------------------------------------
ksp = complex(zeros([1 Nk NiNs Nc], 'single'));

chunk_size = 10000;
nr_chunks = ceil(NiNp / chunk_size);

start_time = tic;
count = 1;
for chunk_nr = 1:nr_chunks
    tstart = tic; fprintf('%s: Reading a chunk of acquisitions (%2d/%2d)... \n', datetime, chunk_nr, nr_chunks);

    %----------------------------------------------------------------------
    % Calculate the index range of a current chunk
    %----------------------------------------------------------------------
    if chunk_nr < nr_chunks
        index_range = (1:chunk_size).' + (chunk_nr - 1) * chunk_size;
    else
        index_range = ((chunk_nr - 1) * chunk_size + 1 : NiNp).';
    end
    current_chunk_size = length(index_range);

    %----------------------------------------------------------------------
    % Read a chunk of acquisitions at a time
    %----------------------------------------------------------------------
    raw_data = dset.readAcquisition(index_range(1), index_range(end));
    scan_counter_list = raw_data.head.scan_counter.';

    for idx = 1:current_chunk_size
        tstart = tic; fprintf('%s:(%2d/%2d) Reading k-space data (%5d/%5d)... ', datetime, chunk_nr, nr_chunks, idx, current_chunk_size);

        %------------------------------------------------------------------
        % Determine the scan number
        %------------------------------------------------------------------
        scan_nr = find(imaging_data_indices == scan_counter_list(idx));

        if ~isempty(scan_nr)
            %--------------------------------------------------------------
            % Prewhiten k-space data
            %--------------------------------------------------------------
            profile = raw_data.data{idx}; % number_of_samples x nr_channels
            profile = (inv_L * profile.').';

            %--------------------------------------------------------------
            % Accumulate k-space
            %--------------------------------------------------------------
            ksp(1,:,count,:) = reshape(profile, [1 Nk 1 Nc]);
            count = count + 1;
        end
        fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
    end
end

%% Perform coil selection
if ~isempty(coil_selection)
    ksp = ksp(:,:,:,coil_selection);
    Nc = length(coil_selection);
end

%% Write a .cfl file
for echo = 1:nr_contrasts
    echo_range = (1:Nk/nr_contrasts) + (echo - 1) * Nk / nr_contrasts;

    %----------------------------------------------------------------------
    % ksp (1 x Nk/nr_contrasts x NiNs x Nc)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path, sprintf('ksp_echo%d', echo));
    tstart = tic; fprintf('%s:(e=%d/%d) Writing a .cfl file: %s... \n', datetime, echo, nr_contrasts, cfl_file);
    writecfl(cfl_file, ksp(:,echo_range,:,:));
    fprintf('%s: done! (%6.4f/%6.4f sec)\n', datetime, toc(tstart), toc(start_time));
end
