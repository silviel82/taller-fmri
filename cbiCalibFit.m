
function cbiCalibFit(calName, varargin)
%
% function cbiCalFit(calName, varargin)
%
% Fit calibration data and produced regularized
% B0, R2*, Coil Sensitivity and EPI correction maps
%
% calName
%   Name of the calibration data file
%
% Optional parameters
%   fitName
%      default: 'cal'
%   regFitName
%      default: fitName_reg
%   SNRCutoff - this is the point at which the weight is 0.5
%      default: 5
%      should work even for very high res sequences
%   regOnlyFlag
%      default: 0
%      if 1, only do the regularization, useful for debugging

% Use the inputParser thingy to deal with the optional inputs
p = inputParser;
p.FunctionName = 'cbiCalibFit';
p.addParamValue('fitName',     'cal',   @ischar);
p.addParamValue('regOnlyFlag',     0,   @(x) x==0 || x==1);
p.parse(varargin{:});

% Set the  intermediate 
fitName     = p.Results.fitName;
regFitName  = [p.Results.fitName '_reg'];
regOnlyFlag = p.Results.regOnlyFlag;

% If we haven't fit the multi-echo data we do it here
if regOnlyFlag == 0
  % Get the delay from the RF to the start of the GRE train
  [pathstr, filestr, ext] = fileparts(calName);
  paramFileName = fullfile(pathstr,[filestr '_params.mat']);
  load(paramFileName)
  to = sp.r.tStart;
  % Fit the data to get the maps
  fprintf('Estimating the maps...\n');
  cbiMultiEchoGREFit(calName,fitName,to);
end

% Compute a mask by truncating the spin density image
nfd = niftifile([fitName '_rho.nii']);
nfd = fopen(nfd,'read');
[nfd, rho] = fread(nfd,nfd.nx*nfd.ny*nfd.nz);
nfd = fclose(nfd);
mask = single((rho>3));
nfd = niftifile([fitName '_mask.nii'],nfd);
nfd = fopen(nfd,'write');
nfd = fwrite(nfd,mask,nfd.nx*nfd.ny*nfd.nz);
nfd = fclose(nfd);


% % Regularize the maps
% Spin density
fprintf('Regularizing the spin density...\n');
cbiRegularizeMap([fitName '_rho.nii'], ...
  [fitName '_mask.nii'], [regFitName '_rho.nii'], 'multi', 4, 0); 
 
% Coil sensitivity
fprintf('Regularizing the coil sensitivity...\n');
cbiRegularizeMap([fitName '_coil.nii'], ...
  [fitName '_mask.nii'], [regFitName '_coil.nii'], 'multi', 8, 0); 
 
% Odd-even ghost
fprintf('Regularizing the odd-even phase difference...\n');
cbiRegularizeMap([fitName '_oe.nii'], ...
  [fitName '_mask.nii'], [regFitName '_oe.nii'], 'multi', 8);

% Bo
fprintf('Regularizing the magnetic field (B0)...\n');
cbiRegularizeMap([fitName '_bo.nii'], ...
  [fitName '_mask.nii'], [regFitName '_bo.nii'], 'multi', 4);      

% R2*
fprintf('Regularizing the decay rate (R2*)...\n');
cbiRegularizeMap([fitName '_rs.nii'], ...
  [fitName '_mask.nii'], [regFitName '_rs.nii'], 'multi', 4);

fprintf('\nDone.\n');

return