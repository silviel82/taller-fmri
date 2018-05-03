%%%%%%%%%%%%%%%%%%%%%%%
% Simulate experiment %
%%%%%%%%%%%%%%%%%%%%%%%

% TR = 1 sec
% nFrames = 480 sec (8 mins)
TR = 1;
nFrames = 2000; 

% trialSequence: randomized sequence of 0,1,2,3
trialSequence = zeros(nFrames,1);
trialSequence(1:nFrames/8) = 1;
trialSequence(nFrames/8+1:nFrames/4) = 2;
trialSequence(nFrames/4+1:3*nFrames/8) = 3;
trialSequence = trialSequence(randperm(nFrames));

% neuralResponse is same as trialSequence
neuralResponse = trialSequence;

% Compute the fMRI responses from the underlying neural responses by
% convolving with the hrf and add noise.
tau = 2;
delta = 2;
lowNoiseSD = 0.1;
highNoiseSD = 1;
fmriResponseLowNoise = hrfconv(neuralResponse,tau,delta,TR) + lowNoiseSD*randn(nFrames,1);
fmriResponse = hrfconv(neuralResponse,tau,delta,TR) + highNoiseSD*randn(nFrames,1);

% Regression analysis
designMatrix = zeros(nFrames,3);
for trialType = 1:3
    designMatrix(:,trialType) = (trialSequence == trialType);
end
for trialType = 1:3
    designMatrix(:,trialType) = hrfconv(designMatrix(:,trialType),tau,delta,TR);
end
% Parameter estimates 
betaLowNoise =  designMatrix \ fmriResponseLowNoise;
beta =  designMatrix \ fmriResponse;
% Model response using betas
modelResponseLowNoise = designMatrix * betaLowNoise;
modelResponse = designMatrix * beta;

% Plot trial sequence and fmriResponse
figure(1); clf;
set(gcf,'DefaultLineLineWidth',2)
time = [100:200];
subplot(3,1,1)
plot(time,trialSequence(time));
ylabel('Neural activity')
subplot(3,1,2)
plot(time,[fmriResponseLowNoise(time) modelResponseLowNoise(time)]);
ylabel('fMRI response (low noise)')
xlabel('Time (sec)')
subplot(3,1,3)
plot(time,[fmriResponse(time) modelResponse(time)]);
ylabel('fMRI response (high noise)')
xlabel('Time (sec)')

% Trial-triggered averages
epochLength = 20;
startTrial = 30;
endTrial = length(trialSequence) - epochLength;
responses = cell(1,3);
modelResponses = cell(1,3);
for trialType = 1:3
    trialNumbers = startTrial + find(trialSequence(startTrial:endTrial) == trialType);
    numTrials = length(trialNumbers);
    epochs = zeros(epochLength,numTrials);
    modelEpochs = zeros(epochLength,numTrials);
    for n=1:numTrials
        trialNum = trialNumbers(n);
        epochs(:,n) = fmriResponse(trialNum:trialNum+epochLength-1);
        modelEpochs(:,n) = modelResponse(trialNum:trialNum+epochLength-1);
    end
    responses{trialType} = epochs;
    modelResponses{trialType} = modelEpochs;
end
meanResponses = zeros(epochLength,3);
semResponses = zeros(epochLength,3);
meanModelResponses = zeros(epochLength,3);
for trialType = 1:3
    meanResponses(:,trialType) = mean(responses{trialType},2);
    semResponses(:,trialType) = std(responses{trialType},1,2) / sqrt(size(responses{trialType},2));
    meanModelResponses(:,trialType) = mean(modelResponses{trialType},2);
end
% Plot it
figure(2); clf; hold on;
set(gcf,'DefaultLineLineWidth',4)
time = [1:epochLength];
plot(time,meanModelResponses(:,1),'k');
%plot(time,meanResponses(:,1),'r');
dsErrorsurface(time,meanResponses(:,1),semResponses(:,1),[1 0 0],0.3);
plot(time,meanModelResponses(:,2),'k');
%plot(time,meanResponses(:,2),'b');
dsErrorsurface(time,meanResponses(:,2),semResponses(:,2),[0 1 0],0.3);
plot(time,meanModelResponses(:,3),'k');
%plot(time,meanResponses(:,3),'b');
dsErrorsurface(time,meanResponses(:,3),semResponses(:,3),[0 0 1],0.3);
hold off
xlabel('Time (sec)');
ylabel('fMRI response');

% Randomization test
[beta,betas,pvalues] = randomizationTest(trialSequence,fmriResponse,tau,delta,TR,1000);
beta
pvalues
figure(3); clf;
subplot(3,1,1)
hist(betas(1,:),20)
subplot(3,1,2)
ylabel('Number of occurrences')
hist(betas(2,:),20)
subplot(3,1,3)
hist(betas(3,:),20)
xlabel('Parameter estimate')

% Covariance and linear xforms on random variables
N = 1000;
data = randn(2,N);
figure(1); clf;
plot(data(1,:),data(2,:),'.')
xlim([-4 4])
ylim([-4 4])

xform = [sqrt(2)/2 sqrt(2)/2; 1 0];
dataXform = xform * data;
figure(2); clf;
plot(dataXform(1,:),dataXform(2,:),'.')
xlim([-4 4])
ylim([-4 4])

covData = (data * data')/N
covDataXform = (dataXform * dataXform')/N

estCov = xform * eye(2) * xform'

