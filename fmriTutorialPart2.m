% fmriTutorialPart2.m
% DJH 1/26/2004

%%%%%%%%%%%%%%%%%
% Preliminaries %
%%%%%%%%%%%%%%%%%

% This tutorial demonstrates how to do basic fMRI data analysis in matlab.
% In the process it provides lots of examples of matlab code that you can
% use to build your own more sophisticated analysis tools. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generating simulated fMRI data set %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Simulating an impulse of neural activity

% Consider a vision experiment in which the stimulus display is blank (mean
% gray) for the first 10 sec, then there is a brief (1 sec) presentation of a
% high contrast visual stimulus, followed once again by a blank screen for
% the remaining 39 sec. If we acquire data with TR=1 sec, then we could
% characterize the stimulus protocol with the following vector:

neuralActivity = zeros(60,1);
neuralActivity(11) = 1;
% Plot it
figure(1); clf;
plot(neuralActivity);
title('Impulse of Neural Activity')
ylabel('Relative neural activity')
xlabel('Time (sec)')

%% The hemodynamic impulse response

% We can think of this as the time-course of the underlying neural
% activity. Because of the sluggish hemodynamic response, we would expect
% the fMRI response to be delayed and spread out over time. Formally, we
% use convolution to characterize the relationship between the fMRI
% response and the underlying neural activity. To simulate the effect of
% the hemodynamics, we use a model for the hemodynamic impulse response.

% The gamma function is one of several popular approximations for the 
% hemodynamic impulse response. There's nothing particularly special about
% the gamma function. There are other functions that do a better job of
% approximating the hemodynamics. But at least it is simple with only 2
% free parameters:
%     tau: time constant (sec)
%     delta: pure delay after stimulus onset (sec)
% Choose some values for these parameters:
tau = 2;
delta = 2;

% Plot the HIRF with these parameter values:
t = [0:1:30];
tshift = max(t-delta,0);
HIRF = (tshift/tau).^2 .* exp(-tshift/tau) / (2*tau);
figure(1); clf;
% Plot it
plot(t,HIRF);
title('Hemodynamic Impulse Response Function')
ylabel('Hemodynamic response')
xlabel('Time (sec)')
% where HIRF stands for "hemodynamic impulse response function".

% Change the values of the parameters and replot the function to see what
% effect its parameter has on the shape of the function.

% For every choice of parameters, the volume of the function is 1 (or very
% close to 1 - it would be exactly one if we sampled more finely than once
% per second):
sum(HIRF)

%% Computing the fMRI response from the neural activity and the
%  hemodynamic impulse response

% Now we use convolution to transform the neuralActivity into an
% fmri signal:
fmriSignal = conv(neuralActivity,HIRF);
fmriSignal = fmriSignal(1:length(neuralActivity));
% Plot it
figure(1); clf;
plot(fmriSignal);
title('fMRI Signal for an Impulse of Neural Activity')
ylabel('fMRI signal (arb units)')
xlabel('Time (sec)')
% Note that the second line is needed because of the way in which matlab's
% 'conv' function works. if you type 'help conv' at your matlab prompt:
%
% >> help conv
%
%  CONV Convolution and polynomial multiplication.
%    C = CONV(A, B) convolves vectors A and B.  The resulting
%    vector is length LENGTH(A)+LENGTH(B)-1.
%  
% Because the output is longer then the input, we throw away the part of
% the fmriSignal that is padded at the end.

%% Baseline
% In real fMRI data, the change in image intensity over time is riding on
% top of a baseline image (see above, Loading and displaying MRI images).
% We typically plot the percent signal change in the image intensity. So
% let's add a value of 100 for the baseline or mean image intensity. And
% then divide by the mean and subtract 1 to convert to percent signal
% change. Note on terminology: I am using 'fmriSignal' to refer to the
% raw image intensity values and 'fmriResponse' to refer to time series
% after they have been converted to units of percent signal change.
baseline = 100;

fmriSignal = baseline + conv(neuralActivity,HIRF);
fmriSignal = fmriSignal(1:length(neuralActivity));
fmriResponse = 100 * ((fmriSignal/(mean(fmriSignal)) - 1));
% Plot it
figure(1); clf;
plot(fmriResponse);
title('fMRI Response for an Impulse of Neural Activity')
ylabel('fMRI response (% change in image intensity)')
xlabel('Time (sec)')
% Note that we are now back pretty close to where we were before; because
% we used a value of 100 for the baseline intensity the original fMRI
% signal was already in units of percent signal change. But note that the
% fMRI response now takes on values slightly less than zero. Why did this
% happen? Is this a problem or is it OK?

%% Sustained response
% See what happens to the fMRI response if we make the stimulus longer so
% that the neural activity is sustained for a longer period of time:
neuralActivity = zeros(80,1);
neuralActivity(21:50) = 1;
fmriSignal = baseline + conv(neuralActivity,HIRF);
fmriSignal = fmriSignal(1:length(neuralActivity));
fmriResponse = 100 * ((fmriSignal/(mean(fmriSignal)) - 1));
% Plot it
figure(1); clf;
subplot(2,1,1)
plot(neuralActivity)
ylim([0,1.1])
title('Sustained Neural Activity')
ylabel('Relative neural activity')
xlabel('Time (sec)')
subplot(2,1,2)
plot(fmriResponse);
title('fMRI Response to Sustained Neural Activity')
ylabel('fMRI response (% change in image intensity)')
xlabel('Time (sec)')
% The peak-to-peak amplitude of the fMRI response was about 1% change in
% image intensity:
max(fmriResponse)-min(fmriResponse)

%% Block alternation

% Now let's simulate the fMRI response for a block-alternation experiment
% in which the visual stimulus contrast (and hence the neural activity)
% changes from low to high every 6 seconds:
t = [1:120];
blockDuration = 6;
neuralActivity = mod(ceil(t/blockDuration),2);
fmriSignal = baseline + conv(neuralActivity,HIRF);
fmriSignal = fmriSignal(1:length(neuralActivity));
fmriResponse = 100 * ((fmriSignal/(mean(fmriSignal)) - 1));
% Plot it
figure(1); clf;
subplot(2,1,1)
plot(neuralActivity)
ylim([0,1.1])
title('Block Alternation in Neural Activity')
ylabel('Relative neural activity')
xlabel('Time (sec)')
subplot(2,1,2)
plot(fmriResponse);
title('fMRI Response to Block Alternation in Neural Activity')
ylabel('fMRI response (% change in image intensity)')
xlabel('Time (sec)')

% There are 10 cycles of neural activity and 10 cycles in the response. But
% that the fMRI response is delayed in time relative to the neural
% activity.

%% The First Cycle
% The first cycle of the response looks different. This is not an artifact
% of the simulation. It is really what one would expect to happen in a real
% fMRI data set because the brain is starting from an extended period with
% no stimulation before the first cycle, whereas there is a history of
% stimulation before each of the other cycles. Because the hemodynamics
% have such a long memory (20-30 sec or so), you have to deal with this
% transient that will occur at the beginning of each experiment. One way to
% deal with this transient in the hemodynamics is to try to explicitly take
% it into account in the analysis (which most analysis programs try to do).
% But, in fact, you have no idea what the subject was doing/thinking before
% the experiment began so there is no way to model it accurately. That's
% why in my lab (DJH), we typically toss the first chunk (20 sec or so) of
% the data to avoid the transient entirely.

% By tossing the first 24 seconds, we get a nice periodic sinusoidal
% response:
neuralActivity = neuralActivity(25:120);
fmriSignal = fmriSignal(25:120);
fmriResponse = 100 * ((fmriSignal/(mean(fmriSignal)) - 1));
% Plot it
figure(1); clf;
subplot(2,1,1)
plot(neuralActivity)
ylim([0,1.1])
title('Block Alternation in Neural Activity')
ylabel('Relative neural activity')
xlabel('Time (sec)')
subplot(2,1,2)
plot(fmriResponse);
title('fMRI Response to Block Alternation in Neural Activity')
ylabel('fMRI response (% change in image intensity)')
xlabel('Time (sec)')

% The amplitude of the fMRI response was about plus or minus 0.2% change in
% image intensity:
max(fmriResponse)
min(fmriResponse)

%% Temporal Fourier transform of the fMRI responses
% (skip this part unless you know what a Fourier transform is)

% Another way to look at the fMRI response is to plot its Fourier
% transform:
nTime = length(neuralActivity);
neuralActivityFFT = abs(fft(neuralActivity)) / (nTime/2);
fmriResponseFFT = abs(fft(fmriResponse)) / (nTime/2);
neuralActivityFFT = neuralActivityFFT(1:nTime/2);
fmriResponseFFT = fmriResponseFFT(1:nTime/2);
frequencies = [0:(nTime/2)-1]/nTime;
% Plot it
figure(2); clf;
subplot(2,1,1)
plot(frequencies,neuralActivityFFT)
title('Fourier Transform of Neural Activity')
ylabel('Relative amplitude of neural activity')
xlabel('Frequency (Hz)')
subplot(2,1,2)
plot(frequencies,fmriResponseFFT);
title('Fourier Transform of fMRI Response')
ylabel('fMRI response amplitude (% change in image intensity)')
xlabel('Frequency (Hz)')

% The (nTime/2) is needed to scale the FFT properly. Note that the response
% amplitude plotted in figure 2 is about 0.2% change in image intensity,
% just what we got before.

% Note that there are lots of frequency components in the block alternation
% of the stimulus and neural activity. But there is really only one
% component in the fMRI response. Why?

% Try varying the 'blockDuration' above to see what that does to the
% response amplitudes, to the transient in the fMRI responses, and to the
% frequency decomposition of the fMRI respones. What happens and why? It
% can be helpful to use blockDurations that divide evenly into the overall
% duration (e.g., 6, 12, 24). Why?

%% Noise

% Go back and regenerate the block-alternation data with blockDuration=6
% (as you did originally before varying the blockDuration).

% Real data, of course, is noisy. We can simulate that too. The
% 'randn' function in matlab produces normally distributed noise with
% mean=0 and variance=1. Make some noise and add it to the responses:
noiseSD = 0.1;
noise = noiseSD * randn(size(fmriSignal));
noisyFmriSignal = fmriSignal + noise;
fmriResponse = 100 * ((noisyFmriSignal/(mean(noisyFmriSignal)) - 1));
% Plot it
figure(1); clf;
subplot(2,1,1)
plot(noise)
title('Noise Added to fMRI Signal')
ylabel('Noise (raw image intensity units)')
xlabel('Time (sec)')
subplot(2,1,2)
plot(fmriResponse);
title('fMRI Response with Noise Added')
ylabel('fMRI response (% change in image intensity)')
xlabel('Time (sec)')
% Note that each time you evaluate the above code, you get a somewhat
% different result because of the noise.

%% Drift

% Real data also exhibits drift, that is largely caused by slow changes in
% the functioning of the MRI scanner. The drift is often modeled as being a
% linear drift, although that is not particularly accurate because
% sometimes it drifts up, sometimes down, and sometimes up and down. But
% we'll stick with the linear drift model for now.
noiseSD = 0.1;
driftRate = 0.01;
noise = noiseSD * randn(size(fmriSignal));
drift = driftRate * [1:length(fmriSignal)];
noisyFmriSignal = fmriSignal + noise + drift;
fmriResponse = 100 * ((noisyFmriSignal/(mean(noisyFmriSignal)) - 1));
% Plot it
figure(1); clf;
subplot(2,1,1)
plot(drift)
title('Drift in fMRI Time Series')
ylabel('Drift (raw image intensity)')
xlabel('Time (sec)')
subplot(2,1,2)
plot(fmriResponse);
title('fMRI Response with Noise and Drift Added')
ylabel('fMRI response (% change in image intensity)')
xlabel('Time (sec)')

% Go back up above and compute the FFT of the new fMRI response (with noise
% and drift added). What does it look like now and why?

%% Generating a whole pile of simulated data

% To explore the statistics of fMRI data analysis, it is useful to generate
% a large number of simulated fMRI time series, one for each of a large
% number of simulated voxels. To keep things simple, the noise will be
% statistically independent and identically distributed. Half of the voxels
% will correspond to tissue that is activated by the block-alternation
% experimental protocol. The other half will not be active. 

% Recompute the block alternation
t = [1:120];
blockDuration = 6;
baseline = 100;
neuralActivity = mod(ceil(t/blockDuration),2);
fmriSignal = baseline + conv(neuralActivity,HIRF);
fmriSignal = fmriSignal(1:length(neuralActivity));

% Make a time series of "images", each with 2000 voxels, half of which will be
% activated and the other half not activated:
nTime = length(fmriSignal);
% Fill nonactive voxels with baseline image intensity
nonactiveVoxels = baseline * ones(nTime,1000);
% Fill active voxels, each with a copy of fmriSignal
activeVoxels = repmat(fmriSignal(:),[1,1000]);

% put the two together, one above the other
data = [activeVoxels nonactiveVoxels];
% The result is a 2d array: nTimePoints x nVoxels
size(data)

% Add noise and drift
noiseSD = 0.1;
driftRate = 0.01;
% add noise
noise = noiseSD * randn(size(data));
data = data + noise;
% add drift
for t=1:nTime
    data(t,:) = data(t,:) + t*driftRate;
end

% Plot the time series for all of the active voxels and all of the inactive
% voxels:
figure(1); clf;
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Introduction to fMRI data analysis %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Regression

% Using linear regression (also called a general linear model or GLM) to
%  estimate model parameters from the data.

% We have a model for how the bold image intensities should change over
% time. The model, in fact, was used to generate the simulated data above.
% The model consists of: 1) the time course of the underlying neural
% activity, 2) the hemodynamic impulse response, 3) linear drift, 4)
% independent and identically distributed (Gaussian) noise. In a real data
% set, there are several unknown parameters in this model. The drift rate
% is unknown and may vary from one voxel to the next. This is not much of a
% problem (as we will see below) as long as it is actually a linear drift.
% But in reality it isn't a linear drift which is more of a problem. The
% parameters of the hemodynamic impulse response are unknown and they may
% vary from subject to subject, session to session, or even from one brain
% area to another. Most data analysis algorithms simply assume a fixed
% hemodynamic impulse response (which is almost certainly incorrect). The
% noise is even more of a problem. The noise, in reality, is neither
% independent (there are both spatial and temporal correlations in the
% noise) nor identically distributed (the noise variance can vary from one
% location to another). Most of the spiffy statistical techniques in fMRI
% data analysis have to do with the fact that the noise is not IID. But we
% will ignore these complications because the noise in our simulated data
% is IID and of known variance. Likewise, our hemodynamic impulse response
% function is also known. So we're in luck.

% Linear regression (GLM) can always be written as a matrix multiplication:
%     y = X b
% where y is a vector of measurements, X is a matrix that consists of the
% model, and b are the unknown parameters of the model. The unknown
% parameters in b are sometimes referred to as "explanatory variables"
% because they are used to explain the variability in the data y. The
% values of b are sometimes referred to as "beta weights".

% The basic idea is to regress the measured time series at a voxel against
% the model. This gives you estimates for the model's unknown parameters
% and it gives you a residual time series (the left over stuff that the
% model failed to account for). If the parameter estimate is large and the
% residuals are small relative to the noise variance, then your voxel is
% behaving as an "active" voxel should. If the both the parameter estimates
% and the residuals are small, then you are looking at a voxel that is not
% exhibiting activity. If the residuals are large, then you probably have a
% lousy model. 

% In our case y is the measured time series at a single voxel, a vector
% with 120 time points. The model in X (often referred to as the design
% matrix because it has to do with the experimental design) is a matrix
% with 3 columns. The first column is the part of the model that accounts
% for the neural activity convolved with the HIRF. The second column is the
% part of the model that accounts for the linear drift, and the third
% column is the part of the model that accounts for the mean baseline image
% intensity. The vector b contains the 3 parameters (beta weights) that we
% want to solve for. Parameter one estimates the amplitude of the
% underlying neural activity. Parameter two estimates the drift rate.
% Parameter 3 estimates the mean image intensity.

% So all we need to do is build the 3 columns of X and use regression to
% solve for b at each voxel.

% Column 1: model of activity. Here again we will compute our model for the
% fMRI signal, computed by convolving our model for the underlying neural
% activity with our model for the hemodynamic impulse response function:
modelActivity = conv(neuralActivity,HIRF);
modelActivity = modelActivity(1:length(neuralActivity));
% This is the same as what was done to produce the simulated data above,
% except that we didn't bother adding the baseline image intensity nor the
% drift (as above) because we will deal with both of those things in the
% other two columns of X.

% Column 2: linear drift.
modelDrift = [1:nTime];

% Column 3: constant, baseline image intensity.
modelConstant = ones(1,nTime);

% Build the design matrix by putting the 3 columns together:
model = [modelActivity(:) modelDrift(:) modelConstant(:)];
% Look at the model
model

% Now estimate b (the beta weights). The equation is:
%    y = X b
% where we want to solve for b. So we compute:
%    b = pinv(X) * y
% where the 'pinv' function in matlab computes the pseudo-inverse of a
% matrix.
nVoxels = size(data,2);
b = zeros(3,nVoxels);
modelInv = pinv(model);
for voxel=1:nVoxels
    b(:,voxel) = modelInv * data(:,voxel);
end

% Note that you can also do this without the loop:
%   b = modelInv * data;

% Plot the parameter estimates
figure(2); clf;
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

% Note that the estimates of neural activity are around 1 for the active
% voxels and around 0 for the nonactive voxels. The estimated drift rates
% are all around 0.01. The mean image intensity estimates are all around
% 100. These are the values that were used above to generate the simulated
% data so we have succeeded in extracting reasonably good estimates for the
% parameters. 

%% Putting error bars on the parameter estimates

% Of course, it's not enough just to estimate the amplitudes of the
% underlying neural activity. Because the data are noisy, we need to put
% error bars on our measurements. Fortunately, matlab provides a function
% called 'regress' that allows you to do this in one step. Typing 'help
% regress' at your matlab prompt gives the following:
%
% >> help regress
% 
%  REGRESS Multiple linear regression using least squares.
%     b = REGRESS(y,X) returns the vector of regression coefficients, b,
%     in the linear model  y = Xb, (X is an nxp matrix, y is the nx1
%     vector of observations).
%  
%     [B,BINT,R,RINT,STATS] = REGRESS(y,X,alpha) uses the input, ALPHA
%     to calculate 100(1 - ALPHA) confidence intervals for B and the
%     residual vector, R, in BINT and RINT respectively.  The vector
%     STATS contains the R-square statistic along with the F and p
%     values for the regression.
%  
%     The X matrix should include a column of ones so that the model
%     contains a constant term.  The F and p values are computed under
%     the assumption that the model contains a constant term, and they
%     are not correct for models without a constant.  The R-square
%     value is the ratio of the regression sum of squares to the
%     total sum of squares.
%
% bint provides a confidence interval for the parameter estimates b. If we
% set alpha=0.05 then this will be a 95% confidence interval.

% Let's try it (but be patient because this function is a bit slow):
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
figure(2); clf;
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
title('Estimates of Drift')
ylabel('Mean image intensity (raw image intensity units)')
xlabel('Position (voxel #)')

% For some (perhaps all) applications of fMRI, this is all that is really
% needed -  estimates for the amplitude of the underlying neural activity,
% ideally with error bars on the estimates. But for some reason, the
% neuroimaging field is obsessed with so-called statistical parameter
% maps...

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Statistical parameter maps %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Histograms of parameter estimates

% Because error bars are difficult to visualize on a two-dimensional image,
% neuroimaging data analysis programs often create statistical-parameter
% maps that represent in one way or another whether or not the activity in
% each voxel is statistically significantly different between two or more
% conditions. This is, of course, equivalent to finding whether or not the
% error bars overlap.

% The key statistical observation is that according to the null hypothesis
% (no activity), the parameter estimates are themselves normally
% distributed because they were computed using linear regression from data
% with IID normally distributed noise. 

% Plot histograms of first parameter estimate, the one that models neural
% activity:
figure(3); clf;
subplot(2,1,1)
hist(b(1,1:1000))
title('Active Voxels')
ylabel('Number of voxels')
xlabel('Estimated neural activity')
subplot(2,1,2)
hist(b(1,1001:2000))
title('Non Active Voxels')
ylabel('Number of voxels')
xlabel('Estimated neural activity')
% The two separate histograms, of course, correspond to the active voxels
% and the non-active voxels. The means of these two distributions are
% different from one another but we only know that because we knew which
% voxels to put in each plot separately. For a real data set, we won't know
% a priori which voxels are active or not. To determine if the parameter
% estimates are different from one another, or different from zero, we need
% to know their standard deviation. And that has to be estimated from the
% data.

%%% Residuals and noise estimation

% If our model (neural activiy, hemodynamics, etc.) is accurate, then it
% should do a pretty good job of fitting the data. To find out how good the
% fit is, we can look at the residuals, that is, the difference between the
% observed data and the model fit at each voxel.
modelPredictions = model * b;
residuals = (data - modelPredictions);
% Plot the SD of the residuals for each voxel:
figure(3); clf;
plot(std(residuals))
title('Residuals of the model fit')
ylabel('SD of residuals (image intensity units)')
xlabel('Position (voxel #)')

%% Distribution of residuals

% The residuals ought to be normally distributed. To check this, we scale
% the residuals by the noise standard deviation that was added in the first
% place:
residualsZscore = residuals(:)/noiseSD;
% Then plot a histogram superimposed with the normal pdf:
nSamples = length(residualsZscore);
step = 0.1;
x = [-5:step:5];
histResiduals = hist(residualsZscore,x)/(step*nSamples);
normalPDF = normpdf(x,0,1);
figure(3); clf;
bar(x,histResiduals);
hold on;
plot(x,normalPDF,'r');
hold off;
title('Residuals Compared with Normal Distribution')
ylabel('Probability/Frequency')
xlabel('Z score')
% And, as you can see, there is a very nice match between the normal pdf
% and the histogram of the residuals. Conclusion: the standard deviation of the
% residuals is a good estimate of the standard deviation of the noise.

%% Transforming the noise estimates

% The residuals tell us how good the fit is and inform us about the
% standard deviation of the noise in the data, but this is not what we
% really need to know. Rather, we need to know how good the parameter
% estimates are. Fortunately, this is pretty simple to characterize. The
% parameter estimates were computed as a weighted sum of the data. See
% above where we computed: b = modelInv * data. Because the data are
% normally distributed and the parameter estimates are a weighted sum of
% the data, then the parameter estimates are also normally distributed.
% Phew. Now all we need is to figure out the standard deviation of the
% parameter estimates. The standard deviation of the parameter estimates
% can be computed from the standard deviations of the residuals by
% subjecting it to the same modelInv transformation.

% Compute the SD of the residuals, separately for each voxel:
residualSD = std(residuals);
residualVar = residualSD.*residualSD;

% Next we will transform the residualSD to compute the SD of the parameter
% estimates. We really care only about the 1st parameter that represents
% the amplitude of neural activity so that's the only one we'll compute. To
% pick this one (while ignoring the others) we define a 'contrast' vector:
c = [1 0 0]';

% Compute the parameter SD for each voxel:
bSD = zeros(size(residualSD));
for voxel=1:nVoxels
    bSD(voxel) = sqrt(c' * modelInv * modelInv' * c * residualVar(voxel));
end
figure(3); clf;
plot(bSD)
title('SD of Estimated Neural Activity')
ylabel('SD of neural activity (arb units)')
xlabel('Position (voxel #)')

% These estimates for the SDs are close to the sample standard deviations,
% computed separately for the active and non-active voxels:
mean(bSD)
std(b(1,1:1000))
std(b(1,1001:2000))
% The 3 numbers that get printed should be pretty close to one another. But
% note that we can compute the sample standard deviations from b only if we
% already know which voxels are active versus inactive, which is why we
% needed to do all the work computing bSD from the residuals.

% Now we have everything we need and we can do a t-test for each
% voxel. The T statistic is the ratio of the parameter estimate to the
% estimated SD. The p value is computed from the T statistic using the
% cumulative distribution function.
tStat = b(1,:)./bSD;
pValue = 1-tcdf(tStat,117);

% Plot 'em
figure(3); clf;
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

% Count the number of "false alarms" among the pixels that are not
% active. There should be about 50 out of 1000 with pvalues < 0.05 and
% there should be about 10 with pvalues < 0.01.
sum(pValue(1001:2000) < 0.05)
sum(pValue(1001:2000) < 0.01)

