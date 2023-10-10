function create_trj_file(header, trj, trj_fullpath)
% Create a trj file
% Written by Nam Gyun Lee
% Email: namgyunl@usc.edu, ggang56@gmail.com (preferred)
% Started: 07/29/2022, Last modified: 07/29/2022


fid = fopen(trj_fullpath, 'w');

fwrite(fid, header.id(1), 'char');
fwrite(fid, header.id(2), 'char');
fwrite(fid, header.id(3), 'char');
fwrite(fid, header.id(4), 'char');

fwrite(fid, header.lNumEchoes , 'int');
fwrite(fid, header.lNumDirs   , 'int');
fwrite(fid, header.lNumSamples, 'int');

fwrite(fid, trj(:), 'double');

fclose(fid);

end