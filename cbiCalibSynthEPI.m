function cbiCalibSynthEPI(rhoName, boName, rsName, TE, kyDir, synthName, upFactor)
% Synthesize an EPI image at the desired TE from calibration data
% rhoName:
%    nifti file rho estimate from calibration scan
% boName:
%    nifti file B0 estimate from calibration scan
%    could be regularized
% rsName:
%    nifti file R2star estimate from calibration scan
%    could be regularized
% kyDir:
%     0 indicates setting B0=0, i.e. T2* only
%    +1 indicates conventional readout
%    -1 indicates reversed readout
% TE:
%    desired TE in units of echo spacing
% upFactor:
%    upsample factor (default=2)
%
% Note:
%   this program assumes that the EPI and the calibration scan have
%   the same acquisition protocol, i.e. FOV, resolution, echo-spacing, 
%   phase encode direction, etc.)
% v. 2.0: May 2013, PV.

% check upFactor:
if (nargin<7)
  warning('Unknown value for upFactor.  Using default (upFactor=2)')
  upFactor = 2;    % default
end


% Load the parameter estimates
nfd = niftifile(rhoName);
nfd = fopen(nfd,'read');
[nfd, rho] = fread(nfd,nfd.nx*nfd.ny*nfd.nz);
rho = reshape(rho, nfd.nx, nfd.ny, nfd.nz);
nfd = fclose(nfd);

nfd = niftifile(boName);
nfd = fopen(nfd,'read');
[nfd, bo] = fread(nfd,nfd.nx*nfd.ny*nfd.nz);
bo = reshape(bo, nfd.nx, nfd.ny, nfd.nz);
nfd = fclose(nfd);

nfd = niftifile(rsName);
nfd = fopen(nfd,'read');
[nfd, rs] = fread(nfd,nfd.nx*nfd.ny*nfd.nz);
rs = reshape(rs, nfd.nx, nfd.ny, nfd.nz);
nfd = fclose(nfd);

% switch the sign of B0 depending on the readout direction
switch kyDir
    case 0
        bo = 0.0*bo;
    case 1
        % do nothing
    case -1
        % reverse
        bo = -bo;
    otherwise
        % error
        error('kyDir should be 0, +1 or -1')
end

% Open the output file
nfd = niftifile(synthName,nfd);
nfd.descrip = 'NYU CBI synthetic EPI image';
nfd = fopen(nfd,'write');

% Forward model and inverse
% see recon code for sign convention, etc.
% we're going to apply the matrices on the right
% The fourier operator
Ny = double(nfd.ny);
% Ny = double(2*nfd.ny);
ky = (-Ny/2 : Ny/2-1)';
y  = (-upFactor*Ny/2 : upFactor*Ny/2-1)'/(upFactor*Ny);
zFy = zeros(upFactor*Ny,Ny);
for n = 1:Ny
   zFy(:,n) = exp(2*pi*1i*ky(n)*y(:));
end
% inverse is the adjoint
zPFy = zFy';
% forward is the integeral
zFy = 1/(upFactor*Ny) * zFy;

Nx = double(nfd.nx);

% build the 2D filter:
filt2D=my2Dfilt(Nx,Ny);
filt2D_up = zeros(upFactor*[Nx,Ny]);
filt2D_up( Nx*(upFactor-1)/2+(1:Nx), ...
           Ny*(upFactor-1)/2+(1:Ny)) = filt2D;

% tStart (first ky line, in units of Echo Spacing):
tStart = TE - double(nfd.ny)/2;

% Loop over the slices and process one at a time
epi = zeros(size(rho));
for nz = 1:nfd.nz
% for nz = 10
    rhoslice = upsample2D( rho(:,:,nz), filt2D_up, upFactor );
    if isreal(rho)
      rhoslice = abs(rhoslice);
    end
    boslice = real( upsample2D( bo(:,:,nz), filt2D_up, upFactor ));
    rsslice =  abs( upsample2D( rs(:,:,nz), filt2D_up, upFactor ));
    data = zeros(upFactor*Nx,Ny);   % this is the signal
    
    brfactor = exp(1i*boslice - rsslice);
    % phase evolution and decay at the beginning of the readout:
    bors = exp((tStart-1)*(1i*boslice-rsslice));    % just 1 ES before beginning of readout.
    
    for n = 1:Ny
        bors = bors.*brfactor;
        data(:,n) = (bors.*rhoslice)*zFy(:,n); % zFy was transposed!      
%         temp = fftshift(fft2(fftshift(bors.*rhoslice)));
%         temp2 = temp(nfd.nx/2+(1:nfd.nx),nfd.ny/2+(1:nfd.ny));
%         data(:,n) = temp2(:,n);%fftshift(fft(fftshift( (bors.*rhoslice)*zFy(:,n), 1),[],1),1); % zFy was transposed!
    end
    % the reconstructed image at double res would be abs(data * zPFy), but
    % we want to downsample and filter it:
%     epi(:,:,nz) = abs(fftshift(ifft2(fftshift(data.*filt2D))));
    hola = filt2D_up.*fftshift(fft2(fftshift(data*zPFy)))/numel(rhoslice);
    epi(:,:,nz) = abs( fftshift(ifft2(fftshift( hola(Nx*(upFactor-1)/2+(1:Nx),Ny*(upFactor-1)/2+(1:Ny)) ))) )*numel(epi(:,:,nz));
end

% Write to disk
nfd = fwrite(nfd, single(abs(epi)), nfd.nx*nfd.ny*nfd.nz);
fclose(nfd);
end  %of the main function

function filt2D=my2Dfilt(Nx,Ny)
  % PV: Filter in 2D:
  f1 = cbiSetKFilt(-Nx/2:Nx/2-1,(7/8)*Nx/2)';
  f2 = cbiSetKFilt(-Ny/2:Ny/2-1,(7/8)*Ny/2);
  filt2D = repmat(f1,1,Ny).*repmat(f2,Nx,1);
end

function myParam2x = upsample2D( myParam, filt2D_up, upFactor )
  fParam = fftshift(fft2(fftshift( myParam )))/numel(myParam);
  fParam2x = zeros(upFactor*size(myParam));
  fParam2x( size(myParam,1)*(upFactor-1)/2+(1:size(myParam,1)), ...
            size(myParam,2)*(upFactor-1)/2+(1:size(myParam,2))) = fParam;
  myParam2x = fftshift(ifft2(fftshift(fParam2x.*filt2D_up)))*numel(fParam2x);
end