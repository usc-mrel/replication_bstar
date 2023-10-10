function save_images_as_dicom_files(twix, dicom_matrix_size, recon_resolution, recon_type, dicom_directory, cfl_file)

start_time = tic;
%% Calculate a reconstructed image for DICOM
%--------------------------------------------------------------------------
% Read in a reconstructed image
%--------------------------------------------------------------------------
img = readcfl(cfl_file);
recon_matrix_size = size(img);

%--------------------------------------------------------------------------
% Set the size of a dicom file
%--------------------------------------------------------------------------
%Nx = 2^(nextpow2(recon_matrix_size(1)));
%Ny = 2^(nextpow2(recon_matrix_size(2)));
%Nz = 2^(nextpow2(recon_matrix_size(2)));
Nx = dicom_matrix_size(1);
Ny = dicom_matrix_size(2);
Nz = dicom_matrix_size(3);

%--------------------------------------------------------------------------
% Zeropad in image-space
%--------------------------------------------------------------------------
N_max = max(recon_matrix_size);
img_zpad = complex(zeros(N_max, N_max, N_max, 'double'));
index2_range = (-floor(recon_matrix_size(2)/2):ceil(recon_matrix_size(2)/2)-1).' + floor(N_max/2) + 1;
img_zpad(:,index2_range,:) = img;

%--------------------------------------------------------------------------
% Crop an image
%--------------------------------------------------------------------------
index1_range = (-floor(Nx/2):ceil(Nx/2)-1).' + floor(N_max/2) + 1;
index2_range = (-floor(Ny/2):ceil(Ny/2)-1).' + floor(N_max/2) + 1;
index3_range = (-floor(Nz/2):ceil(Nz/2)-1).' + floor(N_max/2) + 1;
img_dicom = img_zpad(index1_range, index2_range, index3_range);

%--------------------------------------------------------------------------
% Convert an image from double to uint16
%--------------------------------------------------------------------------
img_dicom = abs(img_dicom) / max(abs(img_dicom(:)));
img_dicom = im2uint16(img_dicom);

%% Write a .dcm file
%--------------------------------------------------------------------------
% Make a dicom directory
%--------------------------------------------------------------------------
mkdir(dicom_directory);

for idx1 = 1:Nx
    %----------------------------------------------------------------------
    % Get study date
    %----------------------------------------------------------------------
    % ExamMemoryUID: '2_0_136761868_20230725_173617_925'
    under_loc = strfind(twix.hdr.Config.ExamMemoryUID, '_');
    StudyDate = sprintf('%s', twix.hdr.Config.ExamMemoryUID(under_loc(3)+1:under_loc(4)-1));

    %----------------------------------------------------------------------
    % Define a dicom filename
    %----------------------------------------------------------------------
    dicom_filename = sprintf('MR_%s_pulseq_bstar_%s_00000_000_000_000_%04d.dcm', lower(recon_type), StudyDate, idx1-1);
    dicom_file = fullfile(dicom_directory, dicom_filename);
    tstart = tic; fprintf('Creating a dicom file: %s... ', dicom_file);

    %----------------------------------------------------------------------
    % Write a dicom file
    %----------------------------------------------------------------------
    dicomwrite(reshape(img_dicom(:,idx1,:), [Nx Nz]), dicom_file);

    %----------------------------------------------------------------------
    % Get a dicom heder
    %----------------------------------------------------------------------
    dicom_info = dicominfo(dicom_file);

    if idx1 == 1
        StudyInstanceUID = dicom_info.StudyInstanceUID;
        SeriesInstanceUID = dicom_info.SeriesInstanceUID;
    end

    % FrameOfReference: '1.3.12.2.1107.5.2.18.41185.1.20230725173617988.0.0.0'
    %'1.3.12.2.1107.5.2.18.41185.   30000023071321585775300000086'
    %'1.3.12.2.1107.5.2.18.41185.   20230725 175246 91595886794.0.0.0'
    %'1.3.12.2.1107.5.2.18.41185.   20230725 174929 74104785507.0.0.0'
    dot_loc = strfind(twix.hdr.Config.FrameOfReference, '.');
    StudyInstanceUID = sprintf('%s.30000023071321585775300000086', twix.hdr.Config.FrameOfReference(1:dot_loc(9)-1));

    %----------------------------------------------------------------------
    % (0008,0020) Study Date
    %----------------------------------------------------------------------
    dicom_info.StudyDate = StudyDate;

    %----------------------------------------------------------------------
    % (0008,0021) Series Date
    %----------------------------------------------------------------------
    dicom_info.SeriesDate = StudyDate;

    %----------------------------------------------------------------------
    % (0008,0022) Acquisition Date
    %----------------------------------------------------------------------
    dicom_info.AcquisitionDate = StudyDate;

    %----------------------------------------------------------------------
    % (0008,0023) Content Date
    %----------------------------------------------------------------------
    dicom_info.ContentDate = StudyDate;

    % ExamMemoryUID: '2_0_104333414_20221031_171518_365'
    studyTime = sprintf('%s.%s', twix.hdr.Config.ExamMemoryUID(under_loc(4)+1:under_loc(5)-1), twix.hdr.Config.ExamMemoryUID(under_loc(end)+1:end));

    %----------------------------------------------------------------------
    % (0008,0030) Study Time
    %----------------------------------------------------------------------
    dicom_info.StudyTime = studyTime;

    %----------------------------------------------------------------------
    % (0008,0031) Series Time
    %----------------------------------------------------------------------
    dicom_info.SeriesTime = studyTime;

    %----------------------------------------------------------------------
    % (0008,0032) Acquisition Time
    %----------------------------------------------------------------------
    dicom_info.AcquisitionTime = studyTime;

    %----------------------------------------------------------------------
    % (0008,0033) Content Time
    %----------------------------------------------------------------------
    dicom_info.ContentTime = studyTime;

    %----------------------------------------------------------------------
    % (0008,0060) Modality
    %----------------------------------------------------------------------
    dicom_info.Modality = twix.hdr.Dicom.Modality;

    %------------------------------------------------------------------
    % (0008,103E) Series Description
    %------------------------------------------------------------------
    dicom_info.SeriesDescription = twix.hdr.Dicom.tProtocolName;

    %----------------------------------------------------------------------
    % (0018,5100) Patient Position
    %----------------------------------------------------------------------
    dicom_info.PatientPosition = twix.hdr.Dicom.tPatientPosition;

    %----------------------------------------------------------------------
    % (0008,0070) Manufacturer
    %----------------------------------------------------------------------
    dicom_info.Manufacturer = twix.hdr.Dicom.Manufacturer;

    %----------------------------------------------------------------------
    % (0008,0080) Institution Name
    %----------------------------------------------------------------------
    dicom_info.InstitutionName = twix.hdr.Dicom.InstitutionName;

    %----------------------------------------------------------------------
    % (0008,0081) Institution Address
    %----------------------------------------------------------------------
    dicom_info.InstitutionAddress = twix.hdr.Dicom.InstitutionAddress;

    %----------------------------------------------------------------------
    % (0008,1010) Station Name
    %----------------------------------------------------------------------
    dicom_info.StationName = sprintf('AWP%d', twix.hdr.Dicom.DeviceSerialNumber);

    %----------------------------------------------------------------------
    % (0008,1040) Institutional Department Name
    %----------------------------------------------------------------------
    dicom_info.InstitutionalDepartmentName = sprintf('MRI Research');

    %----------------------------------------------------------------------
    % (0008,1090) Manufacturer's Model Name
    %----------------------------------------------------------------------
    dicom_info.ManufacturerModelName = twix.hdr.Dicom.ManufacturersModelName;

    %----------------------------------------------------------------------
    % (0018,0084) Imaging Frequency
    %----------------------------------------------------------------------
    dicom_info.ImagingFrequency = twix.hdr.Dicom.lFrequency * 1e-6; % [Hz] * [1e-6MHz/Hz] => *1e-6 [MHz]

    %----------------------------------------------------------------------
    % (0018,0085) Imaged Nucleus
    %----------------------------------------------------------------------
    dicom_info.ImagedNucleus = '1H';

    %----------------------------------------------------------------------
    % (0018,0087) Magnetic Field Strength
    %----------------------------------------------------------------------
    dicom_info.MagneticFieldStrength = twix.hdr.Dicom.flMagneticFieldStrength;

    %----------------------------------------------------------------------
    % (0008,0080) Institution Name
    %----------------------------------------------------------------------
    dicom_info.InstitutionName = twix.hdr.Dicom.InstitutionName;

    %----------------------------------------------------------------------
    % (0018,1000) Device Serial Number
    %----------------------------------------------------------------------
    dicom_info.DeviceSerialNumber = sprintf('%d', twix.hdr.Dicom.DeviceSerialNumber);

    %----------------------------------------------------------------------
    % (0018,1314) Flip Angle
    %----------------------------------------------------------------------
    dicom_info.FlipAngle = twix.hdr.Dicom.adFlipAngleDegree;

    %----------------------------------------------------------------------
    % (0020,0011) Series Number
    %----------------------------------------------------------------------
    dicom_info.SeriesNumber = idx1; % fix? something wrong?

    %----------------------------------------------------------------------
    % (0020,0012) Acquisition Number
    %----------------------------------------------------------------------
    dicom_info.AcquisitionNumber = 1; % fix? something wrong?

    %----------------------------------------------------------------------
    % (0020,0013) Instance Number
    %----------------------------------------------------------------------
    dicom_info.InstanceNumber = idx1; % fix? something wrong?

    %----------------------------------------------------------------------
    % (0010,0020) Patient ID
    %----------------------------------------------------------------------
    dicom_info.PatientID = twix.hdr.Config.PatientID;

    %----------------------------------------------------------------------
    % (0010,0010) Patient's Name
    %----------------------------------------------------------------------
    dicom_info.PatientName.FamilyName = 'Anonymized';

    %------------------------------------------------------------------
    % (0018,1030) Protocol Name
    %------------------------------------------------------------------
    dicom_info.ProtocolName = twix.hdr.Dicom.tProtocolName;

    %------------------------------------------------------------------
    % (0010,0030) Patient's Birth Date
    %------------------------------------------------------------------
    dicom_info.PatientBirthDate = twix.hdr.Config.PatientBirthDay;

    %------------------------------------------------------------------
    % (0010,1030) Patient's Weight
    %------------------------------------------------------------------
    dicom_info.PatientWeight = twix.hdr.Dicom.flUsedPatientWeight;

    %------------------------------------------------------------------
    % (0010,0040) Patient's Sex
    %------------------------------------------------------------------
    if twix.hdr.Dicom.lPatientSex == 1
        PatientSex = 'F'; % female
    elseif twix.hdr.Dicom.lPatientSex == 2
        PatientSex = 'M'; % male
    elseif twix.hdr.Dicom.lPatientSex == 3
        PatientSex = 'O'; % other
    end
    dicom_info.PatientSex = PatientSex;

    %------------------------------------------------------------------
    % (0010,0010) Patient's Age
    %------------------------------------------------------------------
    dicom_info.PatientAge = sprintf('%dY', twix.hdr.Dicom.flPatientAge); % '27Y'

    %------------------------------------------------------------------
    % (0010,0010) Patient's Size
    %------------------------------------------------------------------
    dicom_info.PatientSize = twix.hdr.Meas.flPatientHeight * 1e-3;

    %------------------------------------------------------------------
    % (0010,0015) Body Part Examined
    %------------------------------------------------------------------
    dicom_info.BodyPartExamined = twix.hdr.Dicom.tBodyPartExamined;

    %------------------------------------------------------------------
    % (0020,000D) Study Instance UID
    %------------------------------------------------------------------
    dicom_info.StudyInstanceUID = StudyInstanceUID;

    %------------------------------------------------------------------
    % (0020,000D) Series Instance UID
    %------------------------------------------------------------------
    dicom_info.SeriesInstanceUID = SeriesInstanceUID;

    dicom_info.StudyID                 = '1';
    dicom_info.SeriesNumber            = 2000; % fix? something wrong?
    dicom_info.AcquisitionNumber       = 1;    % fix? something wrong?
    dicom_info.InstanceNumber          = idx1; % fix? something wrong?

    %----------------------------------------------------------------------
    % (0020,0050) Slice Thickness
    %----------------------------------------------------------------------
    dicom_info.SliceThickness = recon_resolution(2) * 1e3;

    %----------------------------------------------------------------------
    % (0020,0032) Image Position (Patient)
    % The x, y, and z coordinates of the upper left hand corner 
    % (center of the first voxel transmitted) of the image, in mm. 
    % See Section C.7.6.2.1.1 for further explanation.    
    %----------------------------------------------------------------------
    % (-floor(Ny/2):ceil(Ny/2)-1) + 0.5
    % Ny = 4 => (-2, -1, 0, 1) + 0.5 => (-1.5, -0.5, 0.5, 1.5)
    dicom_info.ImagePositionPatient = [ recon_resolution(1) * -floor(Nx/2); ...
                                        recon_resolution(2) * (-floor(Ny/2) + 0.5 + idx1 - 1 - 1); ...
                                       -recon_resolution(3) * -floor(Nx/2)] * 1e3;

    %----------------------------------------------------------------------
    % (0020,0037) Image Orientation (Patient)
    %----------------------------------------------------------------------
    dicom_info.ImageOrientationPatient = [1; 0; 0; 0; 0; -1]; % CORONAL

    %----------------------------------------------------------------------
    % (0020,1041) Slice Location
    %----------------------------------------------------------------------
    % (-floor(Ny/2):ceil(Ny/2)-1) + 0.5
    % Ny = 4 => (-2, -1, 0, 1) + 0.5 => (-1.5, -0.5, 0.5, 1.5)
    dicom_info.SliceLocation = recon_resolution(2) * (-floor(Ny/2) + 0.5 + idx1 - 1 - 1) * 1e3;

    %----------------------------------------------------------------------
    % (0028,0030) Pixel Spacing
    %----------------------------------------------------------------------
    dicom_info.PixelSpacing = recon_resolution([1,3]) * 1e3; % [mm]

    %------------------------------------------------------------------
    % Rewrite a dicom file
    %------------------------------------------------------------------
    dicomwrite(reshape(img_dicom(:,idx1,:), [Nx Nz]), dicom_file, dicom_info, 'CreateMode', 'copy');
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
end

end