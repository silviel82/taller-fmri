%%%% simulfMRI.m
%%%% Silvia.10.13.15
% simulfMRI  is a function that simulates fMRI data based on specified
% parameters and runs a regression against a defined model.
% INPUTS:   - noiseSD: a value for the noise standard deviation in the fMRI data
%           - neuralAmp: a value that scales the amplitude of the fMRI data
%           - drift: drift in fMRI signal
%               'linear': linear drift (default)
%               'nonlinear': second order polynomial drift
%           - HIRF: Hemodynamic response function, by default is set to
%               single 'gamma' function.
%               'gamma': single gamma function of the type (t/tau)^2 * e^(t/tau)
%               'gaussian': gaussian function of the type
%                           (1/2pi*sigma) * e^-(x-mu)^2/(2sigma^2)
%           - variable arguments: tau and delta for gamma function
% OUTPUT:   all output are figures
%           - figure(1): HIRF
%           - figure(2): Time series for 1000 active and 1000 non active
%               voxels
%           - figure(3): Parameter estimates for neural activity, drift and
%               baseline
%           - figure(4): Parameter estimates with error bars for neural
%               activity, drift and baseline
%           - figure(5): Histogram of parameter estimates for neural
%               activity of active versus non active voxels
%           - figure(6): Residuals of model fit
%           - figure(7): Distribution of residuals
%           - figure(8): T statistics and p values for each voxels

function simulfMRI(noiseSD,neuralAmp,drift,HIRF,varargin)

% parameters for function inputs
p = inputParser; % handle for function input parser
defaultHIRF = 'gamma'; % set the default HIRF as a gamma
defaultTau = 2; % set default tau
defaultDelta = 2; % set default delta
defaultDrift = 'linear'; % set the default HIRF as a gamma

addRequired(p,'noiseSD',@isnumeric);
addRequired(p,'neuralAmp',@isnumeric);
addOptional(p,'HIRF',defaultHIRF,@ischar);
addOptional(p,'drift',defaultDrift,@ischar);
addParameter(p,'tau',defaultTau);
addParameter(p,'delta',defaultDelta);

parse(p,noiseSD,neuralAmp,drift,HIRF,varargin{:});

% Parameters for block alternation simulation of neural activity:
% block-alternation experiment in which the visual stimulus contrast (and
% hence the neural activity) changes from low to high every 6 seconds:

t = [1:120]; % time of the experiment in seconds
blockDuration = 6;  % duration in seconds on 1 block
% simulate neural activity that goes from 0 to ceiling with every block:
neuralActivity = neuralAmp*mod(ceil(t/blockDuration),2);
% We will assume there is a baseline neural actvity:
baseline = 100;

% parameters for HIRF computation
time = [0:1:30];
tshift = max(time-p.Results.delta,0);
ampscale = 1.9;
sd = 2.1;
mu = 6;
x = linspace(0,30,1000);

% Define HIRF depending to function input
if strcmp(HIRF,'gamma');
    HIRF = (tshift/p.Results.tau).^2 .* exp(-tshift/p.Results.tau) / (2*p.Results.tau);
    figure(1); clf;
    plot(time,HIRF);
    title('Hemodynamic Impulse Response Function')
    ylabel('Hemodynamic response')
    xlabel('Time (sec)')
elseif strcmp(HIRF,'gaussian')
    HIRF = ampscale*(1/(2*pi*sd)*exp(-(x-mu).^2/(2*sd^2)));
    figure(1); clf;
    plot(x,HIRF);
    title('Hemodynamic Impulse Response Function')
    ylabel('Hemodynamic response')
    xlabel('Time (sec)')
end


% Now we will convolve the neural activity with the selected HIRF:
fmriSignal = baseline + conv(neuralActivity,HIRF);
% We will dump the first 24 seconds:
fmriSignal = fmriSignal(25:length(neuralActivity));

% Now we will generate a volume of time series for 2000 voxels, half of
% which will be activated and the other half not activated:
nTime = length(fmriSignal);
% Fill nonactive voxels with baseline image intensity
nonactiveVoxels = baseline * ones(nTime,1000);
% Fill active voxels, each with a copy of fmriSignal
activeVoxels = repmat(fmriSignal(:),[1,1000]);
% put the two together
data = [activeVoxels nonactiveVoxels];

% add noise
noise = noiseSD * randn(size(data));
data = data + noise;
 
% Define and add drift depending to function input
if strcmp(drift,'linear');
    driftRate = 0.01;
    for td=1:nTime
        data(td,:) = data(td,:) + td*driftRate;
    end
elseif strcmp(drift,'nonlinear')
    driftRate = 0.001;
    for td=1:nTime
        data(td,:) = data(td,:) + (td*driftRate + td.^2*driftRate);
    end
end

% OUTPUT PART 1

% Plot the HIRF
figure(1)
plot(time,HIRF);
title('Hemodynamic Impulse Response Function')
ylabel('Hemodynamic response')
xlabel('Time (sec)')

% Plot time series
figure(2); clf;
subplot(2,1,1)
plot(data(:,1:1000));
title('Active Voxels')
ylabel('fMRI response (raw image intensity units)')
xlabel('Time (sec)')
subplot(2,1,2)
plot(data(:,1001:2000));
title('Nonactive Voxels')
ylabel('fMRI response (raw image intensity units)')
xlabel('Time (sec)')

% Now regress the simulated fMRI activity against a generalized linear
% model where our design matrix has a constant and the baseline activity,
% the neural activity, and the drift:

% For the model we will always assume a gamma function, with the following
% parameters:
time = [0:1:30];
tau2 = 2; 
delta2 = 2; 
tshift2 = max(time-delta2,0);
HIRF2 = (tshift2/tau2).^2 .* exp(-tshift2/tau2) / (2*tau2);

% Column 1: neural activity.
modelActivity = conv(neuralActivity,HIRF2);
modelActivity = modelActivity(25:length(neuralActivity));
% Column 2: linear drift.
modelDrift = [1:nTime];
% Column 3: constant, baseline image intensity.
modelConstant = ones(1,nTime);

model = [modelActivity(:) modelDrift(:) modelConstant(:)];

% Now we will estimate the beta weights (b). The equation is:
%    y = X b where we want to solve for b. So we compute:
%    b = pinv(X) * y

nVoxels = size(data,2);
modelInv = pinv(model);
b = modelInv * data;

% OUTPUT PART 2

% Plot the parameter estimates
figure(3); clf;
subplot(3,1,1)
plot(b(1,:))
title('Estimates of Neural Activity')
ylabel('Amplitude of Neural Activity (arb units)')
xlabel('Position (voxel #)')
subplot(3,1,2)
plot(b(2,:))
title('Estimates of Drift')
ylabel('Drift rate (delta image intensity/sec)')
xlabel('Position (voxel #)')
subplot(3,1,3)
plot(b(3,:))
title('Estimates of Baseline Intensity')
ylabel('Mean image intensity (raw image intensity units)')
xlabel('Position (voxel #)')

% Get a confidence interval for the parameter estimates b.
b = zeros(3,nVoxels);
bmin = zeros(3,nVoxels);
bmax = zeros(3,nVoxels);
for voxel=1:nVoxels
    [btmp,bint,r,rint,stats] = regress(data(:,voxel),model,0.05);
    b(:,voxel) = btmp;
    bmin(:,voxel) = bint(:,1);
    bmax(:,voxel) = bint(:,2);
end

% Plot the parameter estimates (from every 20th voxel) with error bars
subVox = [1:20:nVoxels];
figure(4); clf;
subplot(3,1,1)
errorbar(subVox,b(1,subVox),b(1,subVox)-bmin(1,subVox),bmax(1,subVox)-b(1,subVox))
set(gca,'xlim',[-20 2000]);
set(gca,'ylim',[-0.3,1.3]);
title('Estimates of Neural Activity')
ylabel('Amplitude of Neural Activity (arb units)')
xlabel('Position (voxel #)')
subplot(3,1,2)
errorbar(subVox,b(2,subVox),b(2,subVox)-bmin(2,subVox),bmax(2,subVox)-b(2,subVox))
set(gca,'xlim',[-20 2000]);
title('Estimates of Drift')
ylabel('Drift rate (image intensity/sec)')
xlabel('Position (voxel #)')
subplot(3,1,3)
errorbar(subVox,b(3,subVox),b(3,subVox)-bmin(3,subVox),bmax(3,subVox)-b(3,subVox))
set(gca,'xlim',[-20 2000]);
title('Estimates of Baseline Intensity')
ylabel('Mean image intensity (raw image intensity units)')
xlabel('Position (voxel #)')

% Plot the histogram of the active and non active voxel parameter estimates
figure(5); clf;
rgb1 = [82 117 181]./255; 
rgb3 = [248 90 78]./255; 
a = histogram(b(1,1:1000));
hold on
na = histogram(b(1,1001:2000));
a.FaceColor = rgb3;
na.FaceColor = rgb1;
title('Active and Nonactive Voxels')
ylabel('Number of voxels')
xlabel('Estimated neural activity')
legend('active voxels','nonactive voxels');

% Get the residuals of the model (difference between data and predictions
% made by the model, part of the data not explained):
modelPredictions = model * b;
residuals = (data - modelPredictions);

% Plot the SD of the residuals for each voxel:
figure(6); clf;
plot(std(residuals))
title('Residuals of the model fit')
ylabel('SD of residuals (image intensity units)')
xlabel('Position (voxel #)')

% scale the residuals by the noise standard deviation that was added in the
% first place:
residualsZscore = residuals(:)/noiseSD;

% Then plot a histogram of the residuals superimposed with the normal pdf:
nSamples = length(residualsZscore);
step = 0.1;
x = [-5:step:5];
histResiduals = hist(residualsZscore,x)/(step*nSamples);
normalPDF = normpdf(x,0,1);

figure(7); clf;
bar(x,histResiduals);
hold on;
plot(x,normalPDF,'r');
hold off;
title('Residuals Compared with Normal Distribution')
ylabel('Probability/Frequency')
xlabel('Z score')

% We compute the SD of the residuals, separately for each voxel:
residualSD = std(residuals);
residualVar = residualSD.*residualSD;
c = [1 0 0]';
bSD = zeros(size(residualSD));
for voxel=1:nVoxels
    bSD(voxel) = sqrt(c' * modelInv * modelInv' * c * residualVar(voxel));
end

% Compute T statistics using SD of residuals for each voxel, and p values
% for each voxel T stat
tStat = b(1,:)./bSD;
pValue = 1-tcdf(tStat,117);

% Plot the t statistics and p values for each voxel:
figure(8); clf;
subplot(2,1,1)
plot(tStat);
title('T statistic')
ylabel('T value')
xlabel('Position (voxel #)')
subplot(2,1,2)
plot(pValue);
title('P value')
ylabel('P value')
xlabel('Position (voxel #)')

sprintf('the number of active voxels that are statiscally significant is: %d\n',sum(pValue(1:1000) < 0.05))
sprintf('the number of nonactive voxels that are statiscally significant is: %d\n',sum(pValue(1001:2000) < 0.05))
end


