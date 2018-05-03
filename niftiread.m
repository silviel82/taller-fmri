function [data,header,data_vect] = niftiread(filename)
%[data,header,data_vect] = niftiread(filename)
%
%Reads in a nifti file specified by filename and returns the data, the
%header, and the data in a vectorized form
%
%NO x y swap!
%
%written 2009-08-10 by Ed Vessel

nfdin = niftifile(filename);
nfdin = fopen(nfdin,'read');
data = zeros(nfdin.nx, nfdin.ny, nfdin.nz, nfdin.nt);
for i = 1:nfdin.nt
    [nfdin, databuff] = fread(nfdin, nfdin.nx*nfdin.ny*nfdin.nz);
    databuff_reshape = reshape(databuff, [nfdin.nx nfdin.ny nfdin.nz]);
    data(:,:,:,i) = databuff_reshape;
    data_vect(i,:) = databuff_reshape(:);
end

nfdin = fclose(nfdin);
header = nfdin;
end