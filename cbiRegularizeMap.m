function status = cbiRegularizeMap(fmParam, fmWeight, fmReg, fittype, varargin) 
%
% cbiRegularizeMap.m
%
% Purpose:
%   Regularize an estimate from a calibration sequence.
%   Smooth extension to regions where there is not sufficient signal
%
% Usage:
%   [status] = cbiRegularizeMap(fmParam, fmWeight, fmReg, fittype, bias)
%
%  fmParam:
%     nifti file containing parameter estimate
%     e.g. one of the output of cbiMultiEchoGREFit
%  fmWeight:
%     nifti file containing the trust weight
%     e.g. the truncated _relres image from cbiMultiEchoGREFit
%  fmReg:
%     Name of nifti file for the output
%  fittype:
%     poly  Low order polynomial
%     multi Multi-scale
%  scalemin:
%     minimum smoothness scale
%  bias:
%     if bias exists set the background to this value
%  status:
%    return value for error checking
%

%
% Author: Souheil Inati
%         New York University
%

% When run as a program, status = 0 indicates successful completion.
status = 0;

% Check the number of arguments
% if the scalemin was not provided, set it to 2
% if the bias was not provided, set it to Nan
if nargin == 4
  scalemin = 2;
  bias = NaN;
elseif nargin == 5
  scalemin = varargin{1};
  bias = NaN;
else
  scalemin = varargin{1};
  bias = varargin{2};
end


% Open the weight image and read
nfdin = niftifile(fmWeight);
nfdin = fopen(nfdin,'read');
[nfdin,datbuff] = fread(nfdin,nfdin.nx*nfdin.ny*nfdin.nz);
weight = reshape(datbuff,nfdin.nx, nfdin.ny, nfdin.nz);
fclose(nfdin);

% Open the map parameter estimate file for reading
nfdin = niftifile(fmParam);
nfdin = fopen(nfdin,'read');

% Create a new nifti file with the output name and the header of the
% original data.  Everything else about the header for the regularized
% parameter estiamtes is the same
nfdout = niftifile(fmReg,nfdin);

% Open the file
nfdout = fopen(nfdout,'write');

% Loop over the time points in the parameter file and regularize
% each volume

for n = 1:nfdin.nt

  % Read in this volume
  [nfdin,datbuff] = fread(nfdin, nfdin.nx*nfdin.ny*nfdin.nz);
  param = reshape(datbuff, nfdin.nx, nfdin.ny, nfdin.nz);

  % Regularize
  temp = fsmooth(param, weight, nfdin.dx, nfdin.dy, nfdin.dz, fittype, bias, scalemin);
  
  % Write out this regularized volume
  if strcmp(nfdout.datatype,'single') 
    nfdout = fwrite(nfdout, single(temp(:)), nfdin.nx*nfdin.ny*nfdin.nz);
  else
    nfdout = fwrite(nfdout, double(temp(:)), nfdin.nx*nfdin.ny*nfdin.nz);
  end
    
  
end

% Close the output file
fclose(nfdout);

% Close the data file
fclose(nfdin);

end % function cbiRegularizeMap

function rf = fsmooth(f,w,hx,hy,hz, fittype, bias, scalemin)
% Regularization
% f: function to regularize
% w: trust weight
% hx, hy, hz: resolution in x,y,z
%
% Multiscale (hierarchical) weighted smoothing
% at each scale compute:
% the residual from the previous scale (df)
% a blurred version of the weights (sw)
% a blurred version of the weighted residual (swdf)
% Update the estimate
% rf = rf + sw .* swdf ./ (sw.^2 + c^2)
% where c is a regularization parameter
% e.g. if w goes from 0 to 1, then c = .1 is reasonable
%
% Scale tree:
% initialize at infinity (i.e. mean)
% go from 1/4 FOV  (FOV = the smallest of nx*hx, ny*hy, nz*hz)
% down to h
% halving at every step

% Number of pixels
[Nx,Ny,Nz] = size(f);

switch fittype

  case 'poly'
    % done
    % use a low order polynomial
    
    % create a grid in x,y,z
    [x,y,z] = ndgrid( linspace(-1,1,Nx), ...
      linspace(-1,1,Ny), ...
      linspace(-1,1,Nz));

    % flatten to vector
    f = reshape(f,[Nx*Ny*Nz 1]);
    w = reshape(w,[Nx*Ny*Nz 1]);
    x = reshape(x,[Nx*Ny*Nz 1]);
    y = reshape(y,[Nx*Ny*Nz 1]);
    z = reshape(z,[Nx*Ny*Nz 1]);

    % P = 0
    c = sum(f.*w) / sum(w);
    rf = c*ones(size(x));

    % P = 1
    p = x;
    df = f - rf;
    c = sum(df.*w.*p) / sum(w.*(p.^2));
    rf = rf + c*p;

    p = y;
    df = f - rf;
    c = sum(df.*w.*p) / sum(w.*(p.^2));
    rf = rf + c*p;

    p = z;
    df = f - rf;
    c = sum(df.*w.*p) / sum(w.*(p.^2));
    rf = rf + c*p;

    % P = 2
    p = 3/2*x.^2 - 1/2;
    df = f - rf;
    c = sum(df.*w.*p) / sum(w.*(p.^2));
    rf = rf + c*p;

    p = 3/2*y.^2 - 1/2;
    df = f - rf;
    c = sum(df.*w.*p) / sum(w.*(p.^2));
    rf = rf + c*p;

    p = 3/2*z.^2 - 1/2;
    df = f - rf;
    c = sum(df.*w.*p) / sum(w.*(p.^2));
    rf = rf + c*p;

    p = x.*y;
    df = f - rf;
    c = sum(df.*w.*p) / sum(w.*(p.^2));
    rf = rf + c*p;

    p = x.*z;
    df = f - rf;
    c = sum(df.*w.*p) / sum(w.*(p.^2));
    rf = rf + c*p;

    p = y.*z;
    df = f - rf;
    c = sum(df.*w.*p) / sum(w.*(p.^2));
    rf = rf + c*p;

    % P = 3
    p = 5/2*x.^3 - 3/2*x;
    df = f - rf;
    c = sum(df.*w.*p) / sum(w.*(p.^2));
    rf = rf + c*p;

    p = 5/2*y.^3 - 3/2*y;
    df = f - rf;
    c = sum(df.*w.*p) / sum(w.*(p.^2));
    rf = rf + c*p;

    p = 5/2*z.^3 - 3/2*z;
    df = f - rf;
    c = sum(df.*w.*p) / sum(w.*(p.^2));
    rf = rf + c*p;

    p = (3/2*x.^2 - 1/2).*y;
    df = f - rf;
    c = sum(df.*w.*p) / sum(w.*(p.^2));
    rf = rf + c*p;

    p = (3/2*x.^2 - 1/2).*z;
    df = f - rf;
    c = sum(df.*w.*p) / sum(w.*(p.^2));
    rf = rf + c*p;

    p = (3/2*y.^2 - 1/2).*x;
    df = f - rf;
    c = sum(df.*w.*p) / sum(w.*(p.^2));
    rf = rf + c*p;

    p = (3/2*y.^2 - 1/2).*z;
    df = f - rf;
    c = sum(df.*w.*p) / sum(w.*(p.^2));
    rf = rf + c*p;

    p = (3/2*z.^2 - 1/2).*x;
    df = f - rf;
    c = sum(df.*w.*p) / sum(w.*(p.^2));
    rf = rf + c*p;

    p = (3/2*z.^2 - 1/2).*y;
    df = f - rf;
    c = sum(df.*w.*p) / sum(w.*(p.^2));
    rf = rf + c*p;

    % Put back into the grid
    rf = reshape(rf,[Nx Ny Nz]);

  case 'multi'
    % hierarchical fit

    % initial guess - weighted mean or zero
    if isfinite(bias)
      rf = bias*ones(size(f));
    else
      rfinit = mean(f(:).*w(:))/mean(w(:));
      rf = rfinit*ones(size(f));
    end
    
    % The hierarchy of smoothnesses
    minres = min([hx hy hz]);
    scalevec = [10:-2:scalemin];%[10 8 6 4];%2.^([4:-1:1]);

    for scale = 1:length(scalevec)

      % set the width
      width = scalevec(scale)*minres;

      % Calculate the speparable 3D convolution kernel compute the kernel
      % (in pixels) each direction.
      % note: the resolution in mm is hx, hy, hz
      Wh = ceil(width/hx);
      if mod(Wh,2) == 0
        Wh = Wh+1;
      end
      % kernel in x
      % kx = exp(-linspace(-2,2,Wh).^2);
      kx = 1-abs([-(Wh+1)/2:(Wh+1)/2]/Wh);
      kx = kx/sum(kx);                 % normalize
      Wh  = ceil(width/hy);
      if mod(Wh,2) == 0
        Wh = Wh+1;
      end
      % kernel in y
      % ky = exp(-linspace(-2,2,Wh).^2);
      ky = 1-abs([-(Wh+1)/2:(Wh+1)/2]/Wh);
      ky = ky/sum(ky);                 % normalize
      Wh  = ceil(width/hz);
      if mod(Wh,2) == 0
        Wh = Wh+1;
      end
      % kernel in z
      % kz = exp(-linspace(-2,2,Wh).^2);
      kz = 1-abs([-(Wh+1)/2:(Wh+1)/2]/Wh);
      kz = kz/sum(kz);                 % normalize

      % the residual
      drf = f - rf;
      % the initial error
      errprev = norm(w(:).*drf(:))/norm(w(:).*f(:));
      fprintf('  Scale %d, initial error %f\n',scale, errprev);
      
      if errprev < .001
        % if error is small enough skip this scale
        break
      else

        % otherwise, iterate
        for m = 1:10
          % weight the residual
          wdrf = w.*drf;
          % blur the weighted residual
          swdrf = sepconv(wdrf);
          % project
          alpha = swdrf(:)'*wdrf(:) / (swdrf(:)'*swdrf(:));
          % step in that direction
          x = rf + alpha * swdrf;
          % new residual
          drf = f - x;
          % check for convergence
          err = norm(w(:).*drf(:))/norm(w(:).*f(:));
          fprintf('    iter %d, err = %f, delta = %f\n',m, err, err-errprev)
          if  (err-errprev) > -.001
            break
          else
            rf = x;
            errprev = err;
          end
        end

      end
        
    end

end

  function rg = sepconv(g, varargin)
    % Separable convolution

    % Symmetric
    % don't need transp_flag

    % keep track of the input size
    insize = size(g);

    % Reshape onto the grid
    rg  = reshape(g, [Nx Ny Nz]);

    xoff = (length(kx)-1)/2;
    yoff = (length(ky)-1)/2;
    zoff = (length(kz)-1)/2;
    for my = 1:Ny
      for mz = 1:Nz
        temp = conv(squeeze(rg(:,my,mz)),kx);
        rg(:,my,mz) = temp(xoff+(1:Nx));
      end
    end
    for mx = 1:Nx
      for mz = 1:Nz
        temp = conv(squeeze(rg(mx,:,mz)),ky);
        rg(mx,:,mz) = temp(yoff+(1:Ny));
      end
    end
    for mx = 1:Nx
      for my = 1:Ny
        temp = conv(squeeze(rg(mx,my,:)),kz);
        rg(mx,my,:) = temp(zoff+(1:Nz));
      end
    end

    % Reshape back to the input size
    rg  = reshape(rg, insize);

  end

end
