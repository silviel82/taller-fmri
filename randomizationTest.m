function [beta,betas,pvalues] = randomizationTest(trialSequence,fmriResponse,tau,delta,TR,nReps)

nTrialTypes = max(trialSequence);
nFrames = length(trialSequence);
seq = trialSequence;
betas = zeros(nTrialTypes+1,nReps+1);
fmriResponse = fmriResponse - mean(fmriResponse);

for rep = 1:nReps+1
    
    % Build design matrix
    designMatrix = ones(nFrames,nTrialTypes+1);
    for trialType = 1:nTrialTypes
        designMatrix(:,trialType) = (seq == trialType);
    end
    for trialType = 1:nTrialTypes
        designMatrix(:,trialType) = hrfconv(designMatrix(:,trialType),tau,delta,TR);
    end
    
    % Parameter estimates
    beta =  designMatrix \ fmriResponse;
    betas(:,rep) = beta;
    
    % Randomize trial sequence
    seq = trialSequence(randperm(nFrames));

end

beta = betas(:,1);
betas = betas(:,2:end);

pvalues = zeros(nTrialTypes+1,1);
for trialType = 1:nTrialTypes
    pvalues(trialType) = sum(beta(trialType,1) < betas(trialType,:))
end
pvalues = pvalues / nReps;
