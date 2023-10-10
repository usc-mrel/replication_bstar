function [header,trj] = read_trj_file(trj_fullpath)

%% Read a .trj file
%--------------------------------------------------------------------------
% traj*.trj file is the output of the calibration procedure and contains a header,
% and 6 vectors of floats (32 bit) with 416 samples each.
%
% The header format is:
% 
% struct sTrjHeader {
% 
%         char id[4];             // Identification bytes
%         int  lNumEchoes;        // Number of echoes
%         int  lNumDirs;          // Number of directions
%         int  lNumSamples;       // Number of samples per direction
% };
% 
% Followed by 6 float vectors with the information from +X, -X, +Y, -Y, +Z and -Z physical axes.
%--------------------------------------------------------------------------
% Nam's comments
% 1. It seems double instead of float.
% 2. Data storage order:
% +++++++++++  -----------  ...  +++++++++++  -----------  +++++++++++  -----------
% +X (echo 1), -X (echo 1), ..., +Z (echo 1), -Z (echo 1), +X (echo 2)  -X (echo 2), etc
%--------------------------------------------------------------------------
fid = fopen(trj_fullpath, 'r');

id          = fread(fid, 4, 'uchar');
lNumEchoes  = fread(fid, 1, 'int');
lNumDirs    = fread(fid, 1, 'int');
lNumSamples = fread(fid, 1, 'int');

fprintf('id[0/1/2/3] = %d/%d/%d/%d\n', id(1), id(2), id(3), id(4));
fprintf('lNumEchoes  = %3d // Number of echoes\n', lNumEchoes);
fprintf('lNumDirs    = %3d // Number of directions\n', lNumDirs);
fprintf('lNumSamples = %3d // Number of samples per direction\n', lNumSamples);

trj = fread(fid, [lNumSamples lNumDirs * lNumEchoes], '*double', 'ieee-le');

fclose(fid);

%% Create a header structure
header.id          = id;
header.lNumEchoes  = lNumEchoes;
header.lNumDirs    = lNumDirs;
header.lNumSamples = lNumSamples;

end