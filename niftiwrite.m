function [header] = niftiwrite(filename,data,sourceheader)
%Write a nifti file
%
%niftiwrite(filename,data) creates a nifti file with default parameters and
%writes the data.
%
%niftiwrite(filename,data,sourceheader) creates a nifti file with
%parameters inherited from sourceheader and writes the data.
%
%no x y swap!
%
%Written: 2009-08-10 by Ed Vessel

switch nargin
    case 1
        disp('Not enough input arguments!');
    case 2
        nfdout = niftifile(filename);
        [nfdout.nx nfdout.ny nfdout.nz nfdout.nt] = size(data);
    case 3
        nfdout = niftifile(filename,sourceheader);
        nfdout.nt = size(data,4);
end
nfdout.descrip = 'created in matlab by EV niftiwrite';
nfdout.nu = 1; nfdout.nv = 1; nfdout.nw = 1; %TEMP setting these so nvox is computed
nfdout.datatype = class(data);
nfdout = fopen(nfdout,'write'); % at this point the header is written

for i = 1:nfdout.nt
    wbuff =  reshape( squeeze(data(:,:,:,i)), 1, nfdout.nx*nfdout.ny*nfdout.nz);
    [nfdout] = fwrite(nfdout, wbuff, nfdout.nx*nfdout.ny*nfdout.nz);
end



nfdout = fclose(nfdout);
header = nfdout;
end