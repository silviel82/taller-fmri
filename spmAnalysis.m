%%%% spmAnalysis.m
%%%% Silvia.11.2.15
% spmAnalysis  is a function that takes in a single fMRI run, a
% prespecified design matrix and allows for specification of HIRF
% parameters, and returns a volume for parameter estimates from the model
% fit, a volume for t statistics and a volume for p values per voxel.
% INPUTS:   - designMatrix: is a design matrix 
%           - variable arguments: tau and delta for gamma function
% 
% OUTPUT:   - betas:    a 3-D array with the amplitude of neural activity
%                       parameter obtained from the model fit for each voxel
%           - tMaps:    a 3-D array with the t statistic for each voxel
%           - pMaps:    a 3-D array with the p values for each voxel 

function [betas, tMaps, pMaps] = spmAnalysis(designMatrix,epi,varargin)

% parameters for function inputs
p = inputParser; % handle for function input parser
defaultTau = 2; % set default tau
defaultDelta = 0.3; % set default delta

addRequired(p,'designMatrix',@isnumeric);
addRequired(p,'epi',@isnumeric);
addParameter(p,'tau',defaultTau);
addParameter(p,'delta',defaultDelta);

parse(p,designMatrix,epi,varargin{:});

% HIRF computation according to specified or default parameters:
time = [0:2:30];
tshift = max(time-p.Results.delta,0);
HIRF = (tshift/p.Results.tau).^2 .* exp(-tshift/p.Results.tau) / (2*p.Results.tau);
HIRF = HIRF/sum(HIRF);

% Convolve specified designMatrix with computed HIRF for first component of
% model
modelActivity = conv(designMatrix,HIRF);
modelActivity = modelActivity(1:length(designMatrix));
nTime = length(designMatrix);
% Model drift:
modelDrift = [1:nTime];
% Model constant:
modelConstant = ones(1,nTime);
% Model
model = [modelActivity(:) modelDrift(:) modelConstant(:)];

% Extract timeseries from fMRI run
numFrames = size(epi,4); % number of measurements
justVoxels = epi(:,:,:,1); % 3-D volume
epiSize = length(find(justVoxels)); %number of voxels
% create empty matrix to populate
tSeries = zeros(numFrames,epiSize);
% obtain x,y,z coordinates of each voxel in the volume
[x, y, z] = ind2sub(size(justVoxels),find(justVoxels));
% now we create the matrix by finding the timeseries corresponding to each
% of the identified voxels according to their 3-D coordinates
for voxel = 1:epiSize
    tSeries(:,voxel) = squeeze(epi(x(voxel),y(voxel),z(voxel),:));
end
% dump the first cycle (first 5 TRs = 10 seconds)
tSeries = tSeries(6:end,:);

% estimate the beta weights (b). The equation is:
%    y = X b where we want to solve for b. So we compute:
%    b = pinv(X) * y
modelInv = pinv(model);
b = modelInv * tSeries;
% calculate model predictions by multiplying model by weights, then compute
% the residuals (difference between fMRI data and model prediction)
modelPredictions = model * b;
residuals = (tSeries - modelPredictions);
% Compute the SD of the residuals:
residualSD = std(residuals);
residualVar = residualSD.*residualSD;

% Compute the SD for the amplitude of neural activity parameter only (pick
% this one (while ignoring the others) we define a 'contrast' vector):
c = [1 0 0]';
% Compute the parameter SD for each voxel:
bSD = zeros(size(residualSD));
for voxel=1:epiSize
    bSD(voxel) = sqrt(c' * (modelInv * modelInv') * c * residualVar(voxel));
end
% Compute tStat and pValue using the parameter SD
tStat = b(1,:)./bSD;
pValue = 1-tcdf(tStat,117);

% FUNCTION OUTPUT

% 3-D Array of Betas
% Take only the amplitude of neural activity parameter
b1 = b(1,:)';
% populate empty 3-D volume to place each voxel's amplitude parameter
betas = zeros(size(justVoxels));
for voxel= 1:epiSize
    betas(x(voxel),y(voxel),z(voxel)) = b1(voxel);
end

% 3-D Array of tStats
% populate empty 3-D volume to place each voxel's tStat
tStat = tStat';
tMaps = zeros(size(justVoxels));
for voxel= 1:epiSize
    tMaps(x(voxel),y(voxel),z(voxel)) = tStat(voxel);
end


% 3-D Array of p-values
% populate empty 3-D volume to place each voxel's p value
pValue = pValue';
pMaps = zeros(size(justVoxels));
for voxel= 1:epiSize
    pMaps(x(voxel),y(voxel),z(voxel)) = pValue(voxel);
end


end
