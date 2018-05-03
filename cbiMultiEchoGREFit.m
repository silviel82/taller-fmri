function [status] = cbiMultiEchoGREFit(fmName, estName, to, fittype)
%
% status = cbiMultiEchoGREFit(fmName, estName)
%
%   Estimate the effective B0 and R2star from multi-echo GRE data
%   e.g. from images acquired using the cbi_seg_epi_fm sequence.
%
%   fmName
%     name of nifti file containing multiecho images
%   estName
%     basename of nifti file for storing the result
%     will generate files
%         baseName_rho     Spin-density (p)
%         baseName_bo      Magnetic Field (B0)
%         baseName_rs      Relaxation rate (R2* = 1/T2*)
%         baseName_oe      Odd-Even phase difference
%         baseName_mean    Mean image from calibration data
%         baseName_snr     An estimate of how much signal there is
%                          this goes between 0 and 1
%         baseName_err     Residual error in the fit
%         baseName_coils   Relative coil sensitivities
%
%  TO DO
%   Add error checking on the file IO
%   use the status for useful error checking
%

if nargin < 3
  to = 0;
end

if nargin < 4
    fittype = 'adjsimple';
end

status = 0;

% Open the FM data file
nfd = niftifile(fmName);
nfd = fopen(nfd,'read');

% FIXME:
% For now, we only handle the first "contrast"
nfd.nv = 1;


% Initialize the output arrays
rho  = zeros(nfd.nx, nfd.ny, nfd.nz);         % Spin density
bo   = zeros(nfd.nx, nfd.ny, nfd.nz);         % B0
rs   = zeros(nfd.nx, nfd.ny, nfd.nz);         % R2star
oe   = zeros(nfd.nx, nfd.ny, nfd.nz);         % Odd-even phase difference
snr  = zeros(nfd.nx, nfd.ny, nfd.nz);         % Effective SNR
cn   = zeros(nfd.nx, nfd.ny, nfd.nz, nfd.nt); % Relative coil sensitivity
meanimg  = zeros(nfd.nx, nfd.ny, nfd.nz);     % Mean image
synthimg  = zeros(nfd.nx, nfd.ny, nfd.nz);    % Synthetic image

% The number of points and the time in integer units
N = double(nfd.nu);
% The odd number case needs to be handled a bit carefully.
% In this version of the code, we only handle the even number of points.
N = 2*floor(N/2);

% The time vector for fitting
t = (0:N/2-1)';

% Do the calculation one slice at a time
% so that we don't use too much RAM
data = zeros(nfd.nx, nfd.ny, N, nfd.nt, 'single');
for zn = 1:nfd.nz
  % Read the data
  for tn   = 1:N
    for ncoil = 1:nfd.nt
      filepos =  (tn-1)*nfd.nt*nfd.nz*nfd.nx*nfd.ny + ...
                      (ncoil-1)*nfd.nz*nfd.nx*nfd.ny + ...
                                (zn-1)*nfd.nx*nfd.ny;
      nfd = fseek(nfd,filepos); %#ok<LTARG>
      [nfd,datbuff] = fread(nfd,nfd.nx*nfd.ny);
      data(:,:,tn,ncoil) = reshape(datbuff,nfd.nx, nfd.ny);
    end
  end

  % Loop over the pixels and calculate the parameters
  for yn = 1:nfd.ny
    for xn = 1:nfd.nx
      
      % Read the time series (for all coils)
      %   make 2D, (number of time points/2) * (2*number of coils)
      %   i.e. split the odd and even data by coil
      %   this way we have 2*ncoils "channels"
      s = zeros(N/2,2*nfd.nt);
      for ncoil = 1:nfd.nt
        s(:,2*(ncoil-1)+1) = double(squeeze(data(xn,yn,1:2:N,ncoil)));  % odd
        s(:,2*(ncoil-1)+2) = double(squeeze(data(xn,yn,2:2:N,ncoil)));  % even
      end

      % Compute the mean image across time
      mi = 1/N * sqrt(sum(sum(abs(squeeze(data(xn,yn,:,:))).^2,2)));
      
      % If the recon is correct then, the signal is in SNR units
      % and almost all of the background will look rayleigh and have
      % a mean that is less than 3.
      % if (norm(s,'fro')/N/double(nfd.nt) > 3)
      % Skip the blank voxels
      if (norm(s) > 0)
      
        % we have some signal

        if strcmp(fittype, 'hsvd')
            % HSVD fit.  This is expensive!
            
            % Form the Hankel matrix from the data matrix s
            L = floor(N/2/2);
            M = N/2 - L + 1;
            H = zeros(M,2*L*nfd.nt);
            for p = 1:size(s,2)
                for m = 1:L
                    H(:,m+(p-1)*L) = s(m:m+M-1,p);
                end
            end
            [U,S,V] = svd(H,0); %#ok<NASGU>
            
            % The following equivalent to the svd truncated hankel
            % Utop Ubottom, eigenvalue thing for a single component
            r = log((U(1:M-1,1)'*U(2:M,1))/(U(1:M-1,1)'*U(1:M-1,1)));
            
        elseif strcmp(fittype, 'simple')
            % Simple fit
            %   s(2:end,:) = beta*s(1:end-1,:);
            v = s(1:end-1,:); v = v(:);
            w = s(2:end  ,:); w = w(:);
            beta = v'*w / (v'*v);
            r = log(beta);

        else
            % Adapted simple fit            
            v = s(1:end-1,:); v = v(:);
            w = s(2:end  ,:); w = w(:);
            beta = v'*w / (v'*v);
            ri = imag(log(beta));
            sa = mean(sqrt(sum(abs(s(1:3,:)).^2,2)));
            sb = mean(sqrt(sum(abs(s(end-2:end,:)).^2,2)));
            ssum = sum(sqrt(sum(abs(s(2:end-1,:)).^2,2)));
            rr = (sb-sa)/ssum;
            r = rr + 1i*ri;

        end
        
        % Don't allow negative decay rates
        if (real(r)>0)
            r = 1i*imag(r);
        end

        % Estimate the amplitude by solving the linear system
        h = exp(r*t);
        a = h\s;

        % Compute the residual error, relative residual error,
        % and an estimate for the standard deviation on a
        shat = h*a;
        relres = norm(s,'fro')/norm(s-shat,'fro'); % |S| / |S+N|
        relres = relres - 1; % subtract 1 to approximate |S|/|N| 

        % Compute the fitted parameters
        % Field and R2*
        pbo = imag(r/2);  % time between the odd or even echoes
        prs = -real(r/2); % is twice the echo spacing

        % Odd-even fit: correct for bo evolution
        poe = angle(mean(conj(exp(1i*pbo)*a(1:2:end)).*(a(2:2:end))));

        % Individual coil images
        pc = 0.5*a(1:2:end) + 0.5*exp(-1i*(poe+pbo))*a(2:2:end);
        % Remove the phase due to bo since the excitation
        % Keep the amplitude that of the beginning of the readout
        % Remember that where there is no signal estimate of rs blows up
        pc = exp(-1i*pbo*to)*pc;

        % "Spin density" - weighted by the coil sensitivity
        prho = norm(pc);

        % "Synthetic image" - spin density weighted by R2*,
        % averaged across the readout
        if prs > 1e-6
          si = prho.*(1/2/N).*(1-exp(-2*prs*N))./(1-exp(-prs));
        else
          si = prho;
        end
 
      else
        % voxels with signal identically zero
        prho = 0;
        pbo  = 0;
        prs  = 0;
        poe  = 0;
        si = 0;
        relres = 0;
        pc = zeros(nfd.nt,1);
      end

      % Stuff into the output arrays    
      rho(xn,yn,zn)  = prho;
      bo(xn,yn,zn)   = pbo;
      rs(xn,yn,zn)   = prs;
      oe(xn,yn,zn)   = poe;
      meanimg(xn,yn,zn) = mi;
      synthimg(xn,yn,zn) = si;
      snr(xn,yn,zn)  = relres;
      cn(xn,yn,zn,:) = pc;
    
    end % xn loop
  end % yn loop
end % zn loop

% Close the data file
nfd = fclose(nfd);

% Create a new nifti file with the output name of the parameter
% estimated and the header of the original data

% Rho abs
nfdout = niftifile([estName '_rho.nii'],nfd);
nfdout.datatype='single';  % set to single precision
nfdout.nt = 1;             % get rid of the coil number.
nfdout.dt = 1;
nfdout.nu = 1;             % get rid of the number of time points
nfdout.du = 1;
nfdout = fopen(nfdout,'write');
nfdout = fwrite(nfdout, single(rho(:)), nfd.nx*nfd.ny*nfd.nz);
nfdout = fclose(nfdout); %#ok<NASGU>

% Bo
nfdout = niftifile([estName '_bo.nii'],nfd);
nfdout.datatype='single';
nfdout.nt = 1;
nfdout.dt = 1;
nfdout.nu = 1;
nfdout.du = 1;
nfdout = fopen(nfdout,'write');
nfdout = fwrite(nfdout, single(bo(:)), nfd.nx*nfd.ny*nfd.nz);
nfdout = fclose(nfdout); %#ok<NASGU>

% R2*
nfdout = niftifile([estName '_rs.nii'],nfd);
nfdout.datatype='single';
nfdout.nt = 1;
nfdout.dt = 1;
nfdout.nu = 1;
nfdout.du = 1;
nfdout = fopen(nfdout,'write');
nfdout = fwrite(nfdout, single(rs(:)), nfd.nx*nfd.ny*nfd.nz);
nfdout = fclose(nfdout); %#ok<NASGU> 

% Odd-Even correction
nfdout = niftifile([estName '_oe.nii'],nfd);
nfdout.datatype='single';
nfdout.nt = 1;
nfdout.dt = 1;
nfdout.nu = 1;
nfdout.du = 1;
nfdout = fopen(nfdout,'write');
nfdout = fwrite(nfdout, single(oe(:)), nfd.nx*nfd.ny*nfd.nz);
nfdout = fclose(nfdout); %#ok<NASGU>

% Mean
nfdout = niftifile([estName '_mean.nii'],nfd);
nfdout.datatype='single';
nfdout.nt = 1;
nfdout.dt = 1;
nfdout.nu = 1;
nfdout.du = 1;
nfdout = fopen(nfdout,'write');
nfdout = fwrite(nfdout, single(meanimg(:)), nfd.nx*nfd.ny*nfd.nz);
nfdout = fclose(nfdout); %#ok<NASGU>

% Synthetic image
nfdout = niftifile([estName '_synth.nii'],nfd);
nfdout.datatype='single';
nfdout.nt = 1;
nfdout.dt = 1;
nfdout.nu = 1;
nfdout.du = 1;
nfdout = fopen(nfdout,'write');
nfdout = fwrite(nfdout, single(synthimg(:)), nfd.nx*nfd.ny*nfd.nz);
nfdout = fclose(nfdout); %#ok<NASGU>

% Signal to Noise Ratio
nfdout = niftifile([estName '_snr.nii'],nfd);
nfdout.datatype='single';
nfdout.nt = 1;
nfdout.dt = 1;
nfdout.nu = 1;
nfdout.du = 1;
nfdout = fopen(nfdout,'write');
nfdout = fwrite(nfdout, single(snr(:)), nfd.nx*nfd.ny*nfd.nz);
nfdout = fclose(nfdout); %#ok<NASGU>

% Coil Sensitivities
nfdout = niftifile([estName '_coil.nii'],nfd);
nfdout.datatype='complex'; % set to single precision complex
nfdout.nt = nfd.nt;        % the number of coils
nfdout.dt = 1;
nfdout.nu = 1;             % 4D data
nfdout.du = 1;
nfdout = fopen(nfdout,'write');
nfdout = fwrite(nfdout, single(cn(:)), nfd.nx*nfd.ny*nfd.nz*nfd.nt);
nfdout = fclose(nfdout); %#ok<NASGU>

% Done
return


