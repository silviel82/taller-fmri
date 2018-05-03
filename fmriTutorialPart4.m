% fmriTutorialPart4.m
% DJH 1/31/2005

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulate some experiments with different trial sequences %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% TR = 1 sec
% nFrames = 480 sec (8 mins)
TR = 1;
nFrames = 480; 

% trialSequence1
% Stim on for 1s (1 TR) followed by stim off for 5 sec at regular
% intervals. Start times of each trial indicated by 1's.
trialSequence1 = zeros(nFrames,1);
trialSequence1(1:6:nFrames) = 1;
plot(trialSequence1)

% trialSequence2
% Pseudo-randomized version of trialSequence1
trialSequence2 = trialSequence1(randperm(nFrames));
plot(trialSequence2);

% trialSequence3
% Three different trial types in pseudo-random order. Start times of the
% two trial types indicated by 1's, 2's, and 3's.
%    Trial type 1: 2 sec stimulus presentations
%    Trial type 2: 4 sec stimulus presentations
%    Trial type 3: No stimulus
% All trials are separated by a 2 sec inter-trial interval.
%
% So the duration of each trial type is
%    Trial type 1: 4 sec (2 sec stimulus followed by 2 sec ITI)
%    Trial type 2: 6 sec (4 sec stimulus followed by 2 sec ITI)
%    Trial type 3: 2 sec (0 sec stimulus followed by 2 sec ITI)
%
% We want an equal number of ecah trial type. Presenting each trial once
% takes 12 sec, so if the duration of the run is 240 sec, then we have
% time for 20 trials of each type (20*12=240).
%
% To build the trial sequence we start by making a pseudo-random trial order. Then we
% pad with the appropriate number of zeros.
trialOrder = [ones(nFrames/12,1); 2*ones(nFrames/12,1); 3*ones(nFrames/12,1)];
nTrials = length(trialOrder);
trialOrder = trialOrder(randperm(nTrials));
% Note that at this point trialOrder is the same as what is in the stimvec
% file output by the stimulus presentation program.
trialSequence3 = zeros(nFrames,1);
frame = 1;
for trial = 1:nTrials
    trialType = trialOrder(trial);
    trialSequence3(frame) = trialType;
    if trialType == 1
        % duration of trial type 1 is 4 TRs
        frame = frame + 4;
    elseif trialType ==2
        % duration of trial type 2 is 6 TRs
        frame = frame + 6;
    else
        % duration of trial type 3 is 2 TRs
        frame = frame + 2;
    end
end
plot(trialSequence3)

% Compute the underlying neural responses for each of these three
% trialSequences. For trial sequences 1 and 2, this is trivially just a
% scaled copy of the trialSequence. 
neuralResponseAmplitude = 2;
neuralResponse1 = neuralResponseAmplitude * trialSequence1;
neuralResponse2 = neuralResponseAmplitude * trialSequence2;

% For trial type 3, the idea is to incorporate some dynamics (adaptation)
% into the responses. For each trial, there is an intial transient response
% followed by a decay such that for the short trials we have a neural
% response that looks like: 
%    2 1.4 0 0 
% And for the long trials we have neural response that looks like:
%    2 1.4 1.1 1 0 0
% We can express the neural response in matrix notation which is worthwhile
% because we will use the same notation to express the design matrix
% below:
%      (neural response col vector) = (response matrix) * amps
% where response matrix is nFrames X 4 and contains appropriately placed
% 1's and 0's and amps is a vector of 4 numbers: 
amps = [2 1.37 1.14 1.05]';
responseMatrix = zeros(nFrames,4);
frame = 1;
for trial = 1:nTrials
    trialType = trialOrder(trial);
    if trialType == 1
        responseMatrix(frame,1) = 1;
        responseMatrix(frame+1,2) = 1;
        frame = frame + 4;
    elseif trialType ==2
        responseMatrix(frame,1) = 1;
        responseMatrix(frame+1,2) = 1;
        responseMatrix(frame+2,3) = 1;
        responseMatrix(frame+3,4) = 1;
        frame = frame + 6;
    else
        frame = frame + 2;
    end
end
% Compute the simulated neural responese
neuralResponse3 = responseMatrix * amps;
% Plot it
plot(neuralResponse3);
% Print the first part of it so that you can see that it did the right
% thing
neuralResponse3(1:20)


% Compute the fMRI responses from the underlying neural responses by
% convolving with the hrf. For the purpose of this tutorial, we won't
% bother adding baseline or drift (because you already know how to do that
% and what to do about it). 

% HIRF parameters
tau = 2;
delta = 2;

% Noise Parameters
noiseSD = 0.1;

% The three fmri responses
fmriResponse1 = hrfconv(neuralResponse1,tau,delta,TR) + noiseSD*randn(nFrames,1);
fmriResponse2 = hrfconv(neuralResponse2,tau,delta,TR) + noiseSD*randn(nFrames,1);
fmriResponse3 = hrfconv(neuralResponse3,tau,delta,TR) + noiseSD*randn(nFrames,1);
% Plot each of them
plot(fmriResponse1)
plot(fmriResponse2)
plot(fmriResponse3)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Basic regression analysis %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Build the corresponding design matrices. This is actually pretty trivial
% now that we've done all of the work above. The design matrices for
% experiements 1 and 2 are simply hrf convolved versions of the trial
% sequences (note no need for drift and mean in the design matrices);
designMatrix1 = hrfconv(trialSequence1,tau,delta,TR);
designMatrix2 = hrfconv(trialSequence2,tau,delta,TR);

% The design matrix for experiment 3 is an hrf convolved version of
% the above responseMatrix.
designMatrix3 = zeros(nFrames,4);
designMatrix3(:,1) = hrfconv(responseMatrix(:,1),tau,delta,TR);
designMatrix3(:,2) = hrfconv(responseMatrix(:,2),tau,delta,TR);
designMatrix3(:,3) = hrfconv(responseMatrix(:,3),tau,delta,TR);
designMatrix3(:,4) = hrfconv(responseMatrix(:,4),tau,delta,TR);
imagesc(designMatrix3)
colormap gray

% Estimate (using linear regression) the underlying neuralResponse
% amplitudes. For experiments 1 and 2, you should find beta1 = beta2 = 2
% (or at least something close to 2). Why isn't your answer exactly equal
% to 2? For experiment 3, you should find that beta3 is a vector with
% values close to the 4 numbers specified above: [2 1.37 1.14 1.05].
beta1 =  designMatrix1 \ fmriResponse1
beta2 =  designMatrix2 \ fmriResponse2
beta3 =  designMatrix3 \ fmriResponse3


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run a bunch of repeated simulations to characterize the %
% variability in the parameter estimates.                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

numRuns = 1000;
% betas1 = zeros(1,numRuns);
% betas2 = zeros(1,numRuns);
betas3 = zeros(4,numRuns);
for run = 1:numRuns
%     fmriResponse1 = hrfconv(neuralResponse1,tau,delta,TR) + noiseSD*randn(nFrames,1);
%     fmriResponse2 = hrfconv(neuralResponse2,tau,delta,TR) + noiseSD*randn(nFrames,1);
    fmriResponse3 = hrfconv(neuralResponse3,tau,delta,TR) + noiseSD*randn(nFrames,1);
%     betas1(run) =  designMatrix1 \ fmriResponse1;
%     betas2(run) =  designMatrix2 \ fmriResponse2;
    betas3(:,run) =  designMatrix3 \ fmriResponse3;
end

% Evaluate these one at a time to plot histograms, labeled with the mean
% and std of the parameter estimates. Are the mean values what they should
% be or is there bias in the estimates? Is there a difference between the
% reliability of the parameter estimates in the different experiments? If
% so, why?
hist(betas1,50); 
title(['mean = ',num2str(mean(betas1)),'  SD  = ',num2str(std(betas1))]);

hist(betas2,50); 
title(['mean = ',num2str(mean(betas2)),'  SD  = ',num2str(std(betas2))]);

hist(betas3(1,:),50); 
title(['mean = ',num2str(mean(betas3(1,:))),'  SD  = ',num2str(std(betas3(1,:)))]);

hist(betas3(2,:),50); 
title(['mean = ',num2str(mean(betas3(2,:))),'  SD  = ',num2str(std(betas3(2,:)))]);

hist(betas3(3,:),50); 
title(['mean = ',num2str(mean(betas3(3,:))),'  SD  = ',num2str(std(betas3(3,:)))]);

hist(betas3(4,:),50); 
title(['mean = ',num2str(mean(betas3(4,:))),'  SD  = ',num2str(std(betas3(4,:)))]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimating the HIRF in rapid event-related experiments:          %
% Why you have to randomize the trial types, or jitter the timing. %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Experiment 1 was done with regular-spaced trials, all the same as one
% another. What impact does this have on the average response, averaged
% across repeated trials.

% Clip out 24 sec epochs following the start of each trial, ignoring the 1st
% 24 sec to allow the hemodynamics to reach steady state and stopping 24
% sec before the end of the run just to simplify things so that we don't
% need to deal with the partial epochs.
epochDuration = 24;
numTrials1 = sum(trialSequence1(epochDuration:nFrames-epochDuration));
responses1 = zeros(numTrials1,epochDuration);
trialNum = 1;
for frame = epochDuration:nFrames-epochDuration
    if trialSequence1(frame) == 1
        epoch = fmriResponse1(frame:frame+epochDuration-1)';
        responses1(trialNum,:) = epoch;
        trialNum = trialNum + 1;
    end
end
meanResponse1 = mean(responses1);
semResponse1 = std(responses1)/sqrt(numTrials1);
% Plot the mean response with error bars
errorbar([1:epochDuration],meanResponse1,semResponse1);

% Compare this to what happens when the trials are pseudo-randomized in
% experiment 2.
epochDuration = 24;
numTrials2 = sum(trialSequence2(epochDuration:nFrames-epochDuration));
responses2 = zeros(numTrials2,epochDuration);
trialNum = 1;
for frame = epochDuration:nFrames-epochDuration
    if trialSequence2(frame) == 1
        epoch = fmriResponse2(frame:frame+epochDuration-1)';
        responses2(trialNum,:) = epoch;
        trialNum = trialNum + 1;
    end
end
meanResponse2 = mean(responses2);
semResponse2 = std(responses2)/sqrt(numTrials2);
% Plot the mean response with error bars
errorbar([1:epochDuration],meanResponse2,semResponse2);

% Why is it important to randomize? What do you have to randomize? Would it
% be ok to randomize EITHER the trial sequence or the time-intervals
% between trials?


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Use nonlinear least-squares to estimate parameters of HIRF %
% (amplitude, tau and delta) along with estimates of the     %
% parameters characterizing the underlying neural activity.  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% To use the lsqnonlin routine (part of the optimization toolbox)
% we need to define a few things
% X=LSQNONLIN(FUN,X0,LB,UB,OPTIONS,P1,P2,..)
x0(1) = 1; % The first element is the initial guess at the amplitude
x0(2) = 3; % The second element is the initial guess at tau
x0(3) = 3; % The third element is the initial guess at delta

lb = []; % Lower bound
ub = []; % Upper bound
options = optimset('lsqnonlin'); % The default options for the lsqnonlin function

% Parameter estimates for experiment 1.
% You should get an answer close to [2,2,2]. The first number is the
% estimate of neural response amplitude. Second is tau estimate. Third is
% delta estimate.
x1 = lsqnonlin(@fitfunc,x0,lb,ub,options,fmriResponse1,trialSequence1,TR)

% Parameter estimates for experiment 2.
% Again, answer should be close to [2,2,2]
x2 = lsqnonlin(@fitfunc,x0,lb,ub,options,fmriResponse2,trialSequence2,TR)

% For experiment 3 we have to work a bit harder.
% The first four are initial guesses of the neural response parameters.
x0(1) = 2;
x0(2) = 1.4;
x0(3) = 1.1;
x0(4) = 1;
% The next two are initial guesses for tau and delta
x0(5) = 2;
x0(6) = 2;

% Now we need to make sure that all of the parameter estimates are positive
lb = [0 0 0 0 0 0]; % Lower bound
ub = [5 5 5 5 5 10]; % Upper bound

% Get the estimates.
x3 = lsqnonlin(@fitfunc2,x0,lb,ub,options,fmriResponse3,responseMatrix,TR)
% Notice that this might fail (depending on your particular randomization
% of the trial sequence). The fit might hit the maximum number of iterations and
% fail to come up with a reasonable estimate even though we started very
% close to the correct solution. If lsqnonlin does come up with estimates,
% the estimated neural responses might not be monotically decreasing. Why
% should it fail like this?

% Let's try a different function which is constrained to be monotonically
% decreasing a+b*exp(-c*t)
clear x0,x3;
% The first two are initial guesses of the exponential paremeters (a, b,
% and c in the above equation).
x0(1) = 1;
x0(2) = 0.1;
x0(3) = 0.1;
% The next two are initial guesses of tau and delta
x0(4) = 3;
x0(5) = 1;

% Again, we need to make sure that parameter estimates are are positive
lb = [0 0 0 0 0]; % Lower bound
ub = [5 5 5 5 10]; % Upper bound

x3 = lsqnonlin(@fitfunc3,x0,lb,ub,options,fmriResponse3,responseMatrix,TR)

% Check to see of the estimated exponential parameters produce something
% close to the expected fall-off in neural responses: [2 1.37 1.14 1.05]
t = [0 1 2 3];
a = x3(1);
b = x3(2)
c = x3(3)
a + b*exp(-c*t)
