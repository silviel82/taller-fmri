%%%%%%%% Taller fMRI UTP 05.03.18 %%%%%%%%%

%%% creado por David Heeger
%%% actualizado silvia 01.19.16

% Event-Related Design

%%%% PART 1: ROI ANALYSIS %%%%

% Define ROIs
%% ------- Trial-tiggered average --------

%% 1. Load the ROIs from the localizer scans and load the data from the
% event-related scans

% First, define the paths for the program to run, must change according to
% these file locations on each machine
% fmriToolsPath = '/Users/silvia/Dropbox\ \(Personal\)/MATLAB/fmriTools'; 
% dataPath = '"/Users/silvia/Dropbox (Personal)/MATLAB/fmriTools/CSBEventRelated"';
% addpath('/Users/silvia/Dropbox (Personal)/MATLAB/fmriTools/CSBEventRelated/resp');

% Load the event-related fMRI data:

[run1, ~] = niftiread(fullfile(dataPath,'09+Experimentalrun','151201153312.nii'));
[run2, ~] = niftiread(fullfile(dataPath,'10+Experimentalrun','151201153844.nii'));
[run3, ~] = niftiread(fullfile(dataPath,'11+Experimentalrun','151201154512.nii'));
[run4, ~] = niftiread(fullfile(dataPath,'12+Experimentalrun','151201155114.nii'));
[run5, ~] = niftiread(fullfile(dataPath,'13+Experimentalrun','151201155702.nii'));
[run6, ~] = niftiread(fullfile(dataPath,'14+Experimentalrun','151201160235.nii'));
[run7, ~] = niftiread(fullfile(dataPath,'15+Experimentalrun','151201160828.nii'));
[run8, ~] = niftiread(fullfile(dataPath,'16+Experimentalrun','151201161352.nii'));

% Load the masks: VC = visual cortex (Calcarine sulcus, L/R: left/right,
% L/H: low/high threshold
[roiVCLL, ~] = niftiread(fullfile(dataPath,'08+LocalizerRun_2TR','VC_LvR_T8mx.img'));
[roiVCLH, ~] = niftiread(fullfile(dataPath,'08+LocalizerRun_2TR','VC_LvR_T10mx.img'));
[roiVCRL, ~] = niftiread(fullfile(dataPath,'08+LocalizerRun_2TR','VC_RvL_T8mx.img'));
[roiVCRH, ~] = niftiread(fullfile(dataPath,'08+LocalizerRun_2TR','VC_RvL_T10mx.img'));

[roiMCLL, ~] = niftiread(fullfile(dataPath,'08+LocalizerRun_2TR','MC_LvR_T8mx.img'));
[roiMCLH, ~] = niftiread(fullfile(dataPath,'08+LocalizerRun_2TR','MC_LvR_T10mx.img'));
[roiMCRL, ~] = niftiread(fullfile(dataPath,'08+LocalizerRun_2TR','MC_RvL_T8mx.img'));
[roiMCRH, ~] = niftiread(fullfile(dataPath,'08+LocalizerRun_2TR','MC_RvL_T10mx.img'));

% Load the response output (2 columns: first trial types, second responses)

resp_run1 = importdata('respRun1.txt');
resp_run2 = importdata('respRun2.txt');
resp_run3 = importdata('respRun3.txt');
resp_run4 = importdata('respRun4.txt');
resp_run5 = importdata('respRun5.txt');
resp_run6 = importdata('respRun6.txt');
resp_run7 = importdata('respRun7.txt');
resp_run8 = importdata('respRun8.txt');

%% 2. Divide by the mean and convert the time series to percent signal change

% We will use a short function that extracts the run time series from a
% specified ROI and divides by the mean and converts them to percent signal
% change and then computes the mean across the voxels in the ROI.

help roiAnalysis3

% Combine all runs and rois in single arrays
allRuns = {run1;run2;run3;run4;run5;run6;run7;run8};
allRois = {roiVCLL;roiVCLH;roiVCRL;roiVCRH;roiMCLL;roiMCLH;roiMCRL;roiMCRH};

meanRoiTS = cell(length(allRuns),length(allRois));
% array with Runs on the rows and ROIs on the columns, each cell has the
% corresponding timeseries

for j = 1:length(allRois);
    for i = 1:length(allRuns);
        meanRoiTS{i,j} = roiAnalysis3(allRuns{i},allRois{j});
    end
end


%% 3. Remove the drift/trend and get the epochs

% We split the time series into epochs corresponding to each trial and
% subtracting the 1st time point of each epoch from the all of the time
% points in that epoch (i.e., so that each epoch of the time series starts
% at 0).

epochDuration1 = 12; % 24 second epochs
% epochDuration2 = 6; % to try a different epoch duration

ROI_epochs = cell(length(allRuns),length(allRois));
% Data array with Runs on the rows and ROIs on the columns, each cells
% contains epoch data is on rows. each row is a different epoch

% for each ROI and run time series split them into epochs corresponding to
% each trial and subtract the 1st time point
for j = 1:length(allRois);
    for i = 1:length(allRuns);
        sequence = meanRoiTS{i,j};
        nFrames = length(sequence);
        allEpochs = zeros(length(epochDuration1+1:nFrames-epochDuration1),epochDuration1);
        for frame = epochDuration1+1:nFrames-epochDuration1
            allEpochs(frame-epochDuration1,:) = sequence(frame:frame+epochDuration1-1)-sequence(frame);
        end
        
        ROI_epochs{i,j} = allEpochs;
        
    end
end

%% 4. Now split the epochs according to each trial type

% Combine the trial sequence in one array and the responses for each run
% (we'll use them later)

% There are 5 different trial types, numbered 0 through 4, corresponding
% to:
% '0': no stim (just fixation)
% '1': right visual field, left  orientation
% '2': right    "     "  , right     "
% '3': left     "     "  , left      "
% '4': left     "     "  , right     "

allTTypes = {resp_run1.data(:,1);resp_run2.data(:,1);resp_run3.data(:,1);...
    resp_run4.data(:,1);resp_run5.data(:,1);resp_run6.data(:,1);resp_run7.data(:,1)...
    ;resp_run8.data(:,1)};
allTResp = {resp_run1.data(:,2);resp_run2.data(:,2);resp_run3.data(:,2);...
    resp_run4.data(:,2);resp_run5.data(:,2);resp_run6.data(:,2);resp_run7.data(:,2)...
    ;resp_run8.data(:,2)};

trialClass = [0 1 2 3 4];
% Intermediate output array with trial class and all runs 
TT_inRun = cell(length(trialClass),length(allRuns));
% Final output array containing the epochs from each ROI per trial class
TT_Epochs = cell(length(trialClass),length(allRois)); 


for j = 1:length(allRois);
    for k = 1:length(trialClass)
        class = trialClass(k);
        for i = 1:length(allRuns);
            tType = allTTypes{i};
            tType = tType(epochDuration1+1:nFrames-epochDuration1,1);
            currentEpoch = ROI_epochs{i,j};
            
            tIndex = find(tType==class);    
            TT_inRun{k,i} = currentEpoch(tIndex,:); 
             
        end
      TT_Epochs{k,j} = cat(1,TT_inRun{k,:});
    end
    
end

%% 5. For each ROI, compute the mean time series

% 3d array contains the time points of each mean epoch for each ROI per
% trial class.
mean_epochs = zeros(5,12,8);
% mean_epochs = zeros(5,6,8); % for different epoch duration

for j = 1:length(allRois);
    for k = 1:length(trialClass)
        ROI_trialEp = TT_Epochs{k,j};
        mean_epochs(k,:,j) = mean(ROI_trialEp,1);
    end
end

%% 6. For each ROI, compute the SEM time series

% 3d array contains the SEM around each time point of mean epochs for each
% ROI per trial class.
SEM_epochs = zeros(5,12,8);
% SEM_epochs = zeros(5,6,8); % for different epoch duration 

for j = 1:length(allRois);
    for k = 1:length(trialClass)
        ROI_trialEp = TT_Epochs{k,j};
        SEM_epochs(k,:,j) = std(ROI_trialEp,0,1)/sqrt(size(ROI_trialEp,1));
    end
end

%% 7. Plot graphs of the mean time series with error bars


% Plot all the trial triggered averages per ROI
ROIS = {'Left Vis LT','Left Vis HT','Right Vis LT','Right Vis HT',...
    'Left Mot LT','Left Mot HT','Right Mot LT','Right Mot HT'};
x = [1:epochDuration1];

figure('Name','Trial triggered averages for each ROI','NumberTitle','off'); clf;
for r=1:length(allRois);
    subplot(2,4,r) 
    plot(x,mean_epochs(1,:,r),'-k','LineWidth', 1.6);
    hold on 
    plot(x,mean_epochs(2,:,r),'-r','LineWidth', 1.6);
    plot(x,mean_epochs(3,:,r),'-m','LineWidth', 1.6);
    plot(x,mean_epochs(4,:,r),'-b','LineWidth', 1.6);
    plot(x,mean_epochs(5,:,r),'-c','LineWidth', 1.6);
    leg = legend('Null','R-field/L-Ori','R-field/R-Ori','L-field/L-Ori','L-field/R-Ori');
    plotErrorBars(x,mean_epochs(1,:,r),SEM_epochs(1,:,r),'k','k',1);
    plotErrorBars(x,mean_epochs(2,:,r),SEM_epochs(2,:,r),'r','r',1);
    plotErrorBars(x,mean_epochs(3,:,r),SEM_epochs(3,:,r),'m','m',1);
    plotErrorBars(x,mean_epochs(4,:,r),SEM_epochs(4,:,r),'b','b',1);
    plotErrorBars(x,mean_epochs(5,:,r),SEM_epochs(5,:,r),'c','c',1);
    title(sprintf('%s ROI',ROIS{r}),'FontSize',14)
    xlabel ('Time (TRs)','FontSize',12)
    ylabel ('Percent signal change (%)','FontSize',12)
    set (leg,'location', 'northeast','LineWidth', 2)
    legend('boxoff')
    box off
    set(gcf,'color','white')
    set(gca,'tickdir','out','Layer','top', 'LineWidth', 1.6)
end

%% ------- Deconvolution --------

%% 1. Divide by the mean and covert to percent change

% We have done this already by using the function roiAnalysis3 that
% extracts the run time series from a specified ROI and divides by the mean
% and converts them to percent signal change and then computes the mean
% across the voxels in the ROI.

% This array has Runs on the rows and ROIs on the columns, each cell has the
% corresponding timeseries
meanRoiTS;

%% 2. Remove the drift from the data by linear trend removal

% create an empty cell array for detrending
detrendData = cell(length(allRuns),length(allRois));

% for each run and roi detrend using linear regression

for j = 1:length(allRois);
    for i = 1:length(allRuns);
            tSeries = meanRoiTS{i,j}; 
            Drift = [1:nFrames]';    
            driftPar = Drift \ tSeries;
            detrend_tSeries = tSeries - driftPar*Drift; 
            detrendData{i,j} = detrend_tSeries;  
    end
end

%% 3. Linear regression for deconvolution: build design matrix

% First, I'll create a Nx4L deconvolution design matrix for each run
% depending on the trial sequence. N = number of frames (TRs) of each run.
% L = number of time points in the HRF (epoch duration, in this case 12
% TRs).

% We will drop the first 5 frames of each run

% create an empty matrix for all runs
DeconvDMatrix = zeros((150-12)*8,12*4);

% stacked trial sequences (droping 5 first trials of each run):
allTSequence = cat(1,allTTypes{1,1}(13:end,:),allTTypes{2,1}(13:end,:),...
    allTTypes{3,1}(13:end,:),allTTypes{4,1}(13:end,:),allTTypes{5,1}(13:end,:),...
    allTTypes{6,1}(13:end,:),allTTypes{7,1}(13:end,:),allTTypes{8,1}(13:end,:));

% trial sequence matrix where for each column (trial type 1 through 4)
% sequence is denoted by a list of ones.
ttSeq = zeros((150-12)*8,4);
for trialType =1:4;
    ttSeq(:,trialType) = (allTSequence == trialType);
end

% Create a toeplix matrix for each trial type that compose the big
% deconvolution design matrix.
for n = [0:3];
    T = toeplitz(ttSeq(:,n+1));
    m = 1+12*n;
    DeconvDMatrix(:,m:(m+11)) = T(:,1:12);
end
    
%% 4. Concatenate time series and deconvolution matrices

% Concatenate all 8 runs for each ROI and drop first 5 frames of each run,
% combine in a 1160x8 matrix, 1160 frames from all runs, each column
% corresponds to an ROI.

concTSeries = zeros((150-12)*8,8);

for r=1:length(allRois);
    concTSeries(:,r) = cat(1,detrendData{1,r}(13:end,:),detrendData{2,r}(13:end,:),...
    detrendData{3,r}(13:end,:),detrendData{4,r}(13:end,:),detrendData{5,r}(13:end,:),...
    detrendData{6,r}(13:end,:),detrendData{7,r}(13:end,:),detrendData{8,r}(13:end,:));
end

% Regress data against deconvolution matrix and get parameters and
% confidence intervals for that estimation. I then retrieve the standard
% error from those confidence intervals and plot a figure with the
% parameters and their error for comparison with the trial tiggered
% average.

betas = zeros(48,8);
deconvResponses = zeros(4,12,8);
SEMdeconv = zeros(4,12,8);

x = [1:epochDuration1]; 

figure('Name','Deconvolution Results','NumberTitle','off'); clf;
for j = 1:length(allRois);
    m = epochDuration1;
    fmriResponse = concTSeries(:,j);
    [betas(:,j),bint] = regress(fmriResponse,DeconvDMatrix,0.05);
    bmin = bint(:,1); % upper limit of 95% confidence interval from regression
    bmax = bint(:,2); % lower limit of 95% confidence interval from regression
    
    % split parameters into epochs
    deconvResponses(1,:,j) = betas(1:m,j);
    deconvResponses(2,:,j) = betas(m+1:2*m,j);
    deconvResponses(3,:,j) = betas(2*m+1:3*m,j);
    deconvResponses(4,:,j) = betas(3*m+1:end,j);
    
    % Compute standard error from 95% confidence intervals
    SEMs = (bmax - bmin) / 3.92;
    SEMdeconv(1,:,j) = SEMs(1:m,1);
    SEMdeconv(2,:,j) = SEMs(m+1:2*m,1);
    SEMdeconv(3,:,j) = SEMs(2*m+1:3*m,1);
    SEMdeconv(4,:,j) = SEMs(3*m+1:end,1);
  
    subplot(2,4,j) 
    plot(x,deconvResponses(1,:,j)','-r','LineWidth', 1.6);
    hold on
    plot(x,deconvResponses(2,:,j)','-m','LineWidth', 1.6);
    plot(x,deconvResponses(3,:,j)','-b','LineWidth', 1.6);
    plot(x,deconvResponses(4,:,j)','-c','LineWidth', 1.6);
    leg = legend('R-field/L-Ori','R-field/R-Ori','L-field/L-Ori','L-field/R-Ori');
    plotErrorBars(x,deconvResponses(1,:,j),SEMdeconv(1,:,j),'r','r',1);
    plotErrorBars(x,deconvResponses(2,:,j),SEMdeconv(2,:,j),'m','m',1);
    plotErrorBars(x,deconvResponses(3,:,j),SEMdeconv(3,:,j),'b','b',1);
    plotErrorBars(x,deconvResponses(4,:,j),SEMdeconv(4,:,j),'c','c',1);
    title(sprintf('%s ROI',ROIS{j}),'FontSize',14)
    xlabel ('Time (TRs)','FontSize',12)
    ylabel ('beta estimates','FontSize',12)
    set (leg,'location', 'northeast','LineWidth', 2)
    legend('boxoff')
    box off
    set(gcf,'color','white')
    set(gca,'tickdir','out','Layer','top', 'LineWidth', 1.6)
end

%% 5. Compare results of this analysis to trial-triggered average
% see assignment

%% ------- Fit the Data --------

%% 1. build the design matrix, convolve it by a model of the hemodynamics

% define HIRF parameters:
TR =2;
tau=1;
delta=1.;

% Build a stacked design matrix, composed of 138 x 4 matrices,
% nFrames x trial types for each run.
stackedDesignMatrix = zeros((nFrames-12)*8,4);


for j=1:length(allRois);
    for i=1:length(allRuns);
        
        tType = allTTypes{i,1}(13:end,:);
        
        designMatrix = zeros(nFrames-12,4);
        for trialType = 1:4
            designMatrix(:,trialType) = (tType(:,1) == trialType);
        end
        for trialType = 1:4
            designMatrix(:,trialType) = hrfconv(designMatrix(:,trialType),tau,delta,TR);
        end
        stackedDesignMatrix(((nFrames-12)*i-(nFrames-12))+1:(nFrames-12)*i,:) = designMatrix;
    end
end
%% 2. Stacking the time series and corresponding design matrices ...use pseudo-inverse to compute the parameter estimates (beta weights)

% concatenate vertically the design matrix such that all runs are together:
% this results in a matrix that is 1200 x 4, nFrames (all runs) x trial
% types.

stackedTSallROI = cell(1,8);
for j=1:length(allRois);
    for i=1:length(allRuns);
        fmriSignal = detrendData{i,j}(13:end,:);
        stackedTSallROI{1,j} = [stackedTSallROI{1,j}; fmriSignal];
    end
end
% add two columns to design matrix and one for de-meaning (column of ones).
stackedDesignMatrix = [stackedDesignMatrix , ones((nFrames-12)*8,1)];

% compute parameter estimates by stacking the time series for each ROI and
% regressing it against the stacked design matrix.
betasReg = cell(1,8);

for j = 1:length(allRois);
        betasReg{1,j} = stackedDesignMatrix \ stackedTSallROI{1,j};
end

%% 3. Use the parameter estimates to compute a model for the fMRI time-series from each run

% compute the model fmri responses by multiplying the estimated beta
% weights by the design matrix

modelfMRI = cell(1,8);

for j=1:length(allRois);
    modelfMRI{1,j} = stackedDesignMatrix * betasReg{1,j};
end

%% 4. Compute the trial-triggered average on these model response time course

% We will start by de-stacking the model response time courses:
modelbyRun = cell(length(allRuns),length(allRois));

for j=1:length(allRois)
    fmriSignal = modelfMRI{1,j};
    for i = 1:length(allRuns);
        modelbyRun{i,j} = fmriSignal((nFrames-12)*(i-1)+1:i*(nFrames-12));
    end
end

% Break into epochs as we had done before and as we did in
% fmriTutorialPart4: Clip out 24 sec epochs following the start of each
% trial, ignoring the 1st 24 sec to allow the hemodynamics to reach steady
% state and stopping 24 sec before the end of the run just to simplify
% things so that we don't need to deal with the partial epochs.

modelEpochs = cell(length(allRuns),length(allRois));

for j=1:length(allRois)
    for i = 1:length(allRuns);
        sequence = modelbyRun{i,j};
        nFrames = length(sequence); 
       
        mEpochs = zeros(length(epochDuration1+1:nFrames-epochDuration1),epochDuration1);
        for frame = 1:nFrames-epochDuration1
            mEpochs(frame,:) = sequence(frame:frame+epochDuration1-1)'; 
        end
        
        modelEpochs{i,j} = mEpochs;
    end
end

% Now split the epochs into trial type as we had done before:

trialClass = [1 2 3 4];
% Intermediate output array with trial class and all runs 
modTT_inRun = cell(length(trialClass),length(allRuns));
% Final output array containing the epochs from each ROI per trial class
modTT_Epochs = cell(length(trialClass),length(allRois)); 


for j = 1:length(allRois);
    for k = 1:length(trialClass)
        class = trialClass(k);
        for i = 1:length(modelbyRun);
            tType = allTTypes{i,1}(13:end,:);
            tType = tType(epochDuration1+1:nFrames-epochDuration1,1);
            currentEpoch = modelEpochs{i,j};
            
            tIndex = find(tType==class);    
            modTT_inRun{k,i} = currentEpoch(tIndex,:); 
             
        end
      modTT_Epochs{k,j} = cat(1,modTT_inRun{k,:});
    end
    
end

% For each ROI, compute the mean model time series for each epoch per trial
% type and the standard error for each mean time point.

% 3d arrays each contain the time points of each mean epoch and sem for
% each ROI per trial class.
mean_modEpochs = zeros(4,12,8);
SEM_modEpochs = zeros(4,12,8);

for j = 1:length(allRois);
    for k = 1:length(trialClass)
        ROI_tmodlEp = modTT_Epochs{k,j};
        mean_modEpochs(k,:,j) = mean(ROI_tmodlEp,1);
        SEM_modEpochs(k,:,j) = std(ROI_tmodlEp,0,1)/sqrt(size(ROI_tmodlEp,1));
    end
end

% Plot all the model trial triggered averages per ROI
ROIS = {'Left Vis LT','Left Vis HT','Right Vis LT','Right Vis HT',...
    'Left Mot LT','Left Mot HT','Right Mot LT','Right Mot HT'};
x = [1:epochDuration1];

figure('Name','Model response trial triggered averages for each ROI','NumberTitle','off'); clf;
for j=1:length(allRois);
    subplot(2,4,j)
    plot(x,mean_modEpochs(1,:,j),'-r','LineWidth', 1.6);
    hold on 
    plot(x,mean_modEpochs(2,:,j),'-m','LineWidth', 1.6);
    plot(x,mean_modEpochs(3,:,j),'-b','LineWidth', 1.6);
    plot(x,mean_modEpochs(4,:,j),'-c','LineWidth', 1.6);
    leg = legend('R-field/L-Ori','R-field/R-Ori','L-field/L-Ori','L-field/R-Ori');
    plotErrorBars(x,mean_modEpochs(1,:,j),SEM_modEpochs(1,:,j),'r','r',1);
    plotErrorBars(x,mean_modEpochs(2,:,j),SEM_modEpochs(2,:,j),'m','m',1);
    plotErrorBars(x,mean_modEpochs(3,:,j),SEM_modEpochs(3,:,j),'b','b',1);
    plotErrorBars(x,mean_modEpochs(4,:,j),SEM_modEpochs(4,:,j),'c','c',1);
    title(sprintf('%s ROI',ROIS{j}),'FontSize',14)
    xlabel ('Time (TRs)','FontSize',12)
    ylabel ('Percent signal change (%)','FontSize',12)
    set (leg,'location', 'northeast','LineWidth', 2)
    legend('boxoff')
    box off
    set(gcf,'color','white')
    set(gca,'tickdir','out','Layer','top', 'LineWidth', 1.6)
end

%% Compare (by overlaying them on the same graph) the actual and model trial-triggered average plots

% Drop the null trial from the data epochs array for the comparison
mean_Epochs = mean_epochs(2:5,:,:);
SEM_Epochs = SEM_epochs(2:5,:,:);

Trials = {'R-field/L-Ori','R-field/R-Ori','L-field/L-Ori','L-field/R-Ori'};

figure('Name','Comparison data and model trial triggered averages for each ROI','NumberTitle','off'); clf;
% suptitle('data and model comparison of trial triggered averages')
for j=1:length(allRois);
    for trialT=1:4
    subplot(8,4,(j-1)*4+trialT)
    axis([0 14 -1 1.2])
    plot(x,mean_Epochs(trialT,:,j),'-b','LineWidth', 1.6);
     hold on
    plot(x,mean_modEpochs(trialT,:,j),'-m','LineWidth', 1.6);
    plotErrorBars(x,mean_Epochs(trialT,:,j),SEM_Epochs(trialT,:,j),'b','b',1);
    plotErrorBars(x,mean_modEpochs(trialT,:,j),SEM_modEpochs(trialT,:,j),'m','m',1);
    title(sprintf('Roi %s Trial %s',ROIS{j},Trials{trialT}),'FontSize',10)
    box off
    set(gcf,'color','white')
    set(gca,'tickdir','out','Layer','top', 'LineWidth', 1.6)
    end
end 
leg = legend('Data','Model');
set (leg,'Location', 'NorthEastOutside','LineWidth', 2)
legend('boxoff')

%% Separate comparison

figure('Name','Comparison data and model trial triggered averages for each ROI','NumberTitle','off'); clf;
for j=1:length(allRois);
    subplot(4,2,j)
        plot(x,mean_Epochs(1,:,j),'-r','LineWidth', 1.6);
         hold on
        plot(x,mean_Epochs(2,:,j),'-m','LineWidth', 1.6);
        plot(x,mean_Epochs(3,:,j),'-b','LineWidth', 1.6);
        plot(x,mean_Epochs(4,:,j),'-c','LineWidth', 1.6);
        leg = legend('R-field/L-Ori','R-field/R-Ori','L-field/L-Ori','L-field/R-Ori');
        plotErrorBars(x,mean_Epochs(1,:,j),SEM_Epochs(1,:,j),'r','r',1);
        plotErrorBars(x,mean_Epochs(2,:,j),SEM_Epochs(2,:,j),'m','m',1);
        plotErrorBars(x,mean_Epochs(3,:,j),SEM_Epochs(3,:,j),'b','b',1);
        plotErrorBars(x,mean_Epochs(4,:,j),SEM_Epochs(4,:,j),'c','c',1);
        hold on 
        plot(x,mean_modEpochs(1,:,j),'--dr','LineWidth', 3);
        plot(x,mean_modEpochs(2,:,j),'--dm','LineWidth', 3);
        plot(x,mean_modEpochs(3,:,j),'--db','LineWidth', 3);
        plot(x,mean_modEpochs(4,:,j),'--dc','LineWidth', 3);
        title(sprintf('Roi %s',ROIS{j}),'FontSize',14);
        xlabel ('Time (TRs)','FontSize',12)
        ylabel ('Percent signal change (%)','FontSize',12)
        set (leg,'Location', 'NorthEast','LineWidth', 2)
        legend('boxoff')
        box off
        set(gcf,'color','white')
        set(gca,'tickdir','out','Layer','top', 'LineWidth', 1.6)
    
end 

%% 5. Perform a deconvolution of the model response time course. 

modbetas = zeros(48,8);
model_deconvResponses = zeros(4,12,8);
SEMmoddeconv = zeros(4,12,8);

figure('Name','Model Response Deconvolution Results','NumberTitle','off'); clf;
for j=1:length(allRois);
    fmriSignal = modelfMRI{1,j};
    m = epochDuration1;

    [modbetas(:,j),bint] = regress(fmriSignal,DeconvDMatrix,0.05);
    bmin = bint(:,1);
    bmax = bint(:,2);
    
    % split the parameters into epochs
    model_deconvResponses(1,:,j) = modbetas(1:m,j);
    model_deconvResponses(2,:,j) = modbetas(m+1:2*m,j);
    model_deconvResponses(3,:,j) = modbetas(2*m+1:3*m,j);
    model_deconvResponses(4,:,j) = modbetas(3*m+1:end,j);
    
    % Compute standard error from 95% confidence intervals
    SEMs = (bmax - bmin) / 3.92;
    SEMmoddeconv(1,:,j) = SEMs(1:m,1);
    SEMmoddeconv(2,:,j) = SEMs(m+1:2*m,1);
    SEMmoddeconv(3,:,j) = SEMs(2*m+1:3*m,1);
    SEMmoddeconv(4,:,j) = SEMs(3*m+1:end,1);
    
    % Plot it
    subplot(2,4,j) 
    plot(x,model_deconvResponses(1,:,j)','-r','LineWidth', 1.6);
    hold on
    plot(x,model_deconvResponses(2,:,j)','-m','LineWidth', 1.6);
    plot(x,model_deconvResponses(3,:,j)','-b','LineWidth', 1.6);
    plot(x,model_deconvResponses(4,:,j)','-c','LineWidth', 1.6);
    leg = legend('R-field/L-Ori','R-field/R-Ori','L-field/L-Ori','L-field/R-Ori');
    plotErrorBars(x,model_deconvResponses(1,:,j),SEMmoddeconv(1,:,j),'r','r',1);
    plotErrorBars(x,model_deconvResponses(2,:,j),SEMmoddeconv(2,:,j),'m','m',1);
    plotErrorBars(x,model_deconvResponses(3,:,j),SEMmoddeconv(3,:,j),'b','b',1);
    plotErrorBars(x,model_deconvResponses(4,:,j),SEMmoddeconv(4,:,j),'c','c',1);
    title(sprintf('%s ROI',ROIS{j}),'FontSize',14)
    xlabel ('Time (TRs)','FontSize',12)
    ylabel ('beta estimates','FontSize',12)
    set (leg,'location', 'northeast','LineWidth', 2)
    legend('boxoff')
    box off
    set(gcf,'color','white')
    set(gca,'tickdir','out','Layer','top', 'LineWidth', 1.6)
end


%% Separate comparison

figure('Name','Comparison data and model deconvolutions for each ROI','NumberTitle','off'); clf;

for j=1:length(allRois);
    subplot(4,2,j) 
        plot(x,deconvResponses(1,:,j),'-r','LineWidth', 1.6);
         hold on
        plot(x,deconvResponses(2,:,j),'-m','LineWidth', 1.6);
        plot(x,deconvResponses(3,:,j),'-b','LineWidth', 1.6);
        plot(x,deconvResponses(4,:,j),'-c','LineWidth', 1.6);
        leg = legend('R-field/L-Ori','R-field/R-Ori','L-field/L-Ori','L-field/R-Ori');
        set(leg,'Location', 'NorthEast','LineWidth', 2);
        legend('boxoff')
        plotErrorBars(x,deconvResponses(1,:,j),SEMdeconv(1,:,j),'r','r',1);
        plotErrorBars(x,deconvResponses(2,:,j),SEMdeconv(2,:,j),'m','m',1);
        plotErrorBars(x,deconvResponses(3,:,j),SEMdeconv(3,:,j),'b','b',1);
        plotErrorBars(x,deconvResponses(4,:,j),SEMdeconv(4,:,j),'c','c',1);
        hold on 
        plot(x,model_deconvResponses(1,:,j),'--dr','LineWidth', 3);
        plot(x,model_deconvResponses(2,:,j),'--dm','LineWidth', 3);
        plot(x,model_deconvResponses(3,:,j),'--db','LineWidth', 3);
        plot(x,model_deconvResponses(4,:,j),'--dc','LineWidth', 3);
        title(sprintf('Roi %s',ROIS{j}),'FontSize',14);
        set(gca,'tickdir','out','Layer','top', 'LineWidth', 1.6)
        xlabel ('Time (TRs)','FontSize',12)
        ylabel ('Percent signal change (%)','FontSize',12)
        box off
        set(gcf,'color','white')
        
end

%% Compare (by overlaying them on the same graph) the actual and model deconvolutions

Trials = {'R-field/L-Ori','R-field/R-Ori','L-field/L-Ori','L-field/R-Ori'};

figure('Name','Comparison data and model deconvolutions for each ROI','NumberTitle','off'); clf;
% suptitle('data and model comparison of trial triggered averages')
for j=1:length(allRois);
    for trialT=1:4
    subplot(8,4,(j-1)*4+trialT)
%     axis([0 14 -1 1.2])
    plot(x,deconvResponses(trialT,:,j),'-b','LineWidth', 1.6);
     hold on
    plot(x,model_deconvResponses(trialT,:,j),'-m','LineWidth', 1.6);
    plotErrorBars(x,deconvResponses(trialT,:,j),SEMdeconv(trialT,:,j),'b','b',1);
    plotErrorBars(x,model_deconvResponses(trialT,:,j),SEMmoddeconv(trialT,:,j),'m','m',1);
    title(sprintf('Roi %s Trial %s',ROIS{j},Trials{trialT}),'FontSize',10)
    box off
    set(gcf,'color','white')
    set(gca,'tickdir','out','Layer','top', 'LineWidth', 1.6)
    end
end 
leg = legend('Data','Model');
set (leg,'Location', 'NorthEastOutside','LineWidth', 2)
legend('boxoff')


%%
%%%% PART 2: STATISTICAL PARAMETER MAPPING %%%%

% Use the same stacked design matrix we used for regression but only the
% trial sequence columns.
stackedSPMMatrix = stackedDesignMatrix(:,1:4);

% We will stack the time series but will convert them to percent signal
% change and detrend by regression first:
stackedRuns = zeros(size(allRuns{1,1},1),size(allRuns{1,1},2),...
    size(allRuns{1,1},3),length(allRuns{1,1}(:,:,:,13:end))*8);

Drift = [1:length(allRuns{1,1}(:,:,:,13:end))]';
Model = [ones(length(allRuns{1,1}(:,:,:,13:end)),1) Drift];


% this for loop is slow!
tic;
for r=1:length(allRuns)
    data = allRuns{r,1}(:,:,:,13:end);
    for x=1:size(data,1)
        for y=1:size(data,2)
            for z=1:size(data,3)
                tSeries = squeeze(data(x,y,z,:)); 
                baseline = mean(tSeries);
                percentTseries = 100 * (tSeries/baseline - 1);
                betas = Model \ percentTseries; 
                Tseries = percentTseries - Model*betas;
                stackedRuns(x,y,z,(r-1)*size(data,4)+1:r*size(data,4)) = tSeries;
                
            end
        end
    end
    
end
toc;

% For regression, estimate the beta weights using the pseudo inverse
% Calculate pseudo inverse of model
modelInv = pinv(stackedSPMMatrix);
% contrast matrix for computing the SD of each parameter
c1 = [1 0 0 0]';
c2 = [0 1 0 0]';
c3 = [0 0 1 0]';
c4 = [0 0 0 1]';

% create empty matrices for t Stat and p value maps for each trial type:
tMapTT1 = zeros(size(stackedRuns(:,:,:,1)));
tMapTT2 = zeros(size(tMapTT1));
tMapTT3 = zeros(size(tMapTT1));
tMapTT4 = zeros(size(tMapTT1));

pMapTT1 = zeros(size(tMapTT1));
pMapTT2 = zeros(size(tMapTT1));
pMapTT3 = zeros(size(tMapTT1));
pMapTT4 = zeros(size(tMapTT1));

% WARNING: This for loop is really really slow, takes about 2.5 minutes!
tic;
for x = 1:size(stackedRuns,1)
    for y = 1:size(stackedRuns,2)
        for z = 1:size(stackedRuns,3)
            tSeries = squeeze(stackedRuns(x,y,z,:));
            % regression b = pinv(X) * y
            bet = stackedSPMMatrix \ tSeries;
            % calculate model predictions by multiplying model by weights, then compute
            % the residuals (difference between fMRI data and model prediction)
            modelPrediction = stackedSPMMatrix * bet;
            residuals = tSeries - modelPrediction;
            residualSD = std(residuals);
            residualsVar = residualSD*residualSD;
            
            b1SD = sqrt(c1' * (modelInv * modelInv') * c1 * residualsVar);
            b2SD = sqrt(c2' * (modelInv * modelInv') * c2 * residualsVar);
            b3SD = sqrt(c3' * (modelInv * modelInv') * c3 * residualsVar);
            b4SD = sqrt(c4' * (modelInv * modelInv') * c4 * residualsVar);
            
            tMapTT1(x,y,z) = bet(1)./b1SD;
            tMapTT2(x,y,z) = bet(2)./b2SD;
            tMapTT3(x,y,z) = bet(3)./b3SD;
            tMapTT4(x,y,z) = bet(4)./b4SD;
     
            pMapTT1(x,y,z) = 1-tcdf(tMapTT1(x,y,z),nFrames-4); % number of degrees of freedom = n - 4 parameters
            pMapTT2(x,y,z) = 1-tcdf(tMapTT2(x,y,z),nFrames-4);
            pMapTT3(x,y,z) = 1-tcdf(tMapTT3(x,y,z),nFrames-4);
            pMapTT4(x,y,z) = 1-tcdf(tMapTT4(x,y,z),nFrames-4);
           
        end
    end
end
toc;

niftiwrite('tMapTrial1.nii',tMapTT1);
niftiwrite('tMapTrial2.nii',tMapTT2);
niftiwrite('tMapTrial3.nii',tMapTT3);
niftiwrite('tMapTrial4.nii',tMapTT4);

niftiwrite('pMapTrial1.nii',pMapTT1);
niftiwrite('pMapTrial2.nii',pMapTT2);
niftiwrite('pMapTrial3.nii',pMapTT3);
niftiwrite('pMapTrial4.nii',pMapTT4);

%%
%%%% PART 3: ADVANCED TOPICS IN EVENT-RELATED DATA ANALYSIS AND SIMULATION %%%%

%% ------- Simultaneous estimation of HIRF and response amplitudes --------
% We will start by initializing the parameters with a guess at their
% values. Parameters will be in vector beta the first four will correspond
% to the fmri response to each trial and the last two will be tau and delta
% parameters for the HIRF.

% Initialize neural response parameters for each trial
x0(1) = 1;
x0(2) = 1;
x0(3) = 1;
x0(4) = 1;
% Initialize tau and delta parameters for HIRF
x0(5) = 1;
x0(6) = 1;

% Make sure that all of the parameter estimates are positive
lb = [0 0 0 0 0 0]; % Lower bound
ub = [5 5 5 5 5 10]; % Upper bound
options = optimset('lsqnonlin'); % The default options for the lsqnonlin function
nFrames = 150;

nlDesignMatrix = zeros((nFrames)*8,4);
for i=1:length(allRuns);
        
    tType = allTTypes{i,1};
        
    designMatrix = zeros(nFrames,4);
    for trialType = 1:4
         designMatrix(:,trialType) = (tType(:,1) == trialType);
    end
        
    nlDesignMatrix(((nFrames)*i-(nFrames))+1:(nFrames)*i,:) = designMatrix;
end

nlstackedTS = cell(1,8);
for j=1:length(allRois);
    for i=1:length(allRuns);
        fmriSignal = detrendData{i,j};
        nlstackedTS{1,j} = [nlstackedTS{1,j}; fmriSignal];
    end
end

roi1 = lsqnonlin(@fitfunc2,x0,lb,ub,options,nlstackedTS{1,1},nlDesignMatrix,TR);
roi2 = lsqnonlin(@fitfunc2,x0,lb,ub,options,nlstackedTS{1,2},nlDesignMatrix,TR);
roi3 = lsqnonlin(@fitfunc2,x0,lb,ub,options,nlstackedTS{1,3},nlDesignMatrix,TR);
roi4 = lsqnonlin(@fitfunc2,x0,lb,ub,options,nlstackedTS{1,4},nlDesignMatrix,TR);
roi5 = lsqnonlin(@fitfunc2,x0,lb,ub,options,nlstackedTS{1,5},nlDesignMatrix,TR);
roi6 = lsqnonlin(@fitfunc2,x0,lb,ub,options,nlstackedTS{1,6},nlDesignMatrix,TR);
roi7 = lsqnonlin(@fitfunc2,x0,lb,ub,options,nlstackedTS{1,7},nlDesignMatrix,TR);
roi8 = lsqnonlin(@fitfunc2,x0,lb,ub,options,nlstackedTS{1,8},nlDesignMatrix,TR);

nlResults = [roi1;roi2;roi3;roi4;roi5;roi6;roi7;roi8];

betasnlReg = cell(1,8);
for j = 1:length(allRois);
        betasnlReg{1,j} = nlResults(j,1:4);
end

Tau = [roi1(5) roi2(5) roi3(5) roi4(5) roi5(5) roi6(5) roi7(5) roi8(5)];
Delta = [roi1(6) roi2(6) roi3(6) roi4(6) roi5(6) roi6(6) roi7(6) roi8(6)];

avgTau = mean([roi1(5) roi2(5) roi3(5) roi4(5) roi5(5) roi6(5) roi7(5) roi8(5)]);
avgDelta = mean([roi1(6) roi2(6) roi3(6) roi4(6) roi5(6) roi6(6) roi7(6) roi8(6)]);

% Build a stacked design matrix, composed of 138 x 4 matrices,
% nFrames x trial types for each run.

stackednlDesignMatrix = cell(1,8);
internlDesignMatrix = zeros((nFrames)*8,4);

for j=1:length(allRois);
    for i=1:length(allRuns);
        
        tType = allTTypes{i,1};
        
        designMatrix = zeros(nFrames,4);
        for trialType = 1:4
            designMatrix(:,trialType) = (tType(:,1) == trialType);
        end
        for trialType = 1:4
            designMatrix(:,trialType) = hrfconv(designMatrix(:,trialType),Tau(j),Delta(j),TR);
%             designMatrix(:,trialType) = hrfconv(designMatrix(:,trialType),avgTau,avgDelta,TR);
        end
        internlDesignMatrix(((nFrames)*i-(nFrames))+1:(nFrames)*i,:) = designMatrix;
    end
    
    stackednlDesignMatrix{1,j} = internlDesignMatrix;
end
% compute the model fmri responses by multiplying the estimated beta
% weights by the design matrix

modelnlfMRI = cell(1,8);

for j=1:length(allRois);
    modelnlfMRI{1,j} = stackednlDesignMatrix{1,j} * betasnlReg{1,j}';
end


% We will start by de-stacking the model response time courses:
modelnlbyRun = cell(length(allRuns),length(allRois));

for j=1:length(allRois)
    fmriSignal = modelnlfMRI{1,j};
    for i = 1:length(allRuns);
        modelnlbyRun{i,j} = fmriSignal((nFrames)*(i-1)+1:i*(nFrames));
    end
end

% Break into epochs as we had done before and as we did in
% fmriTutorialPart4: Clip out 24 sec epochs following the start of each
% trial, ignoring the 1st 24 sec to allow the hemodynamics to reach steady
% state and stopping 24 sec before the end of the run just to simplify
% things so that we don't need to deal with the partial epochs.

modelnlEpochs = cell(length(allRuns),length(allRois));

for j=1:length(allRois)
    for i = 1:length(allRuns);
        sequence = modelnlbyRun{i,j};
        nFrames = length(sequence); 
       
        nlmEpochs = zeros(length(epochDuration1+1:nFrames-epochDuration1),epochDuration1);
        for frame = epochDuration1+1:nFrames-epochDuration1
            nlmEpochs(frame-epochDuration1,:) = sequence(frame:frame+epochDuration1-1)-sequence(frame); 
        end
        
        modelnlEpochs{i,j} = nlmEpochs;
    end
end

% Now split the epochs into trial type as we had done before:

trialClass = [1 2 3 4];
% Intermediate output array with trial class and all runs 
nlmodTT_inRun = cell(length(trialClass),length(allRuns));
% Final output array containing the epochs from each ROI per trial class
nlmodTT_Epochs = cell(length(trialClass),length(allRois)); 


for j = 1:length(allRois);
    for k = 1:length(trialClass)
        class = trialClass(k);
        for i = 1:length(modelnlbyRun);
            tType = allTTypes{i,1};
            tType = tType(epochDuration1+1:nFrames-epochDuration1,1);
            currentEpoch = modelnlEpochs{i,j};
            
            tIndex = find(tType==class);    
            nlmodTT_inRun{k,i} = currentEpoch(tIndex,:); 
             
        end
      nlmodTT_Epochs{k,j} = cat(1,nlmodTT_inRun{k,:});
    end
    
end

% For each ROI, compute the mean model time series for each epoch per trial
% type and the standard error for each mean time point.

% 3d arrays each contain the time points of each mean epoch and sem for
% each ROI per trial class.
mean_modnlEpochs = zeros(4,12,8);
SEM_modnlEpochs = zeros(4,12,8);

for j = 1:length(allRois);
    for k = 1:length(trialClass)
        ROI_tmodlEp = nlmodTT_Epochs{k,j};
        mean_modnlEpochs(k,:,j) = mean(ROI_tmodlEp,1);
        SEM_modnlEpochs(k,:,j) = std(ROI_tmodlEp,0,1)/sqrt(size(ROI_tmodlEp,1));
    end
end

% Plot all the model trial triggered averages per ROI
ROIS = {'Left Vis LT','Left Vis HT','Right Vis LT','Right Vis HT',...
    'Left Mot LT','Left Mot HT','Right Mot LT','Right Mot HT'};
x = [1:epochDuration1];

figure('Name','NL Model response trial triggered averages for each ROI','NumberTitle','off'); clf;
for j=1:length(allRois);
    subplot(2,4,j)
    plot(x,mean_modnlEpochs(1,:,j),'-r','LineWidth', 1.6);
    hold on 
    plot(x,mean_modnlEpochs(2,:,j),'-m','LineWidth', 1.6);
    plot(x,mean_modnlEpochs(3,:,j),'-b','LineWidth', 1.6);
    plot(x,mean_modnlEpochs(4,:,j),'-c','LineWidth', 1.6);
    leg = legend('R-field/L-Ori','R-field/R-Ori','L-field/L-Ori','L-field/R-Ori');
    plotErrorBars(x,mean_modnlEpochs(1,:,j),SEM_modnlEpochs(1,:,j),'r','r',1);
    plotErrorBars(x,mean_modnlEpochs(2,:,j),SEM_modnlEpochs(2,:,j),'m','m',1);
    plotErrorBars(x,mean_modnlEpochs(3,:,j),SEM_modnlEpochs(3,:,j),'b','b',1);
    plotErrorBars(x,mean_modnlEpochs(4,:,j),SEM_modnlEpochs(4,:,j),'c','c',1);
    title(sprintf('%s ROI',ROIS{j}),'FontSize',14)
    xlabel ('Time (TRs)','FontSize',12)
    ylabel ('Percent signal change (%)','FontSize',12)
    set (leg,'location', 'NorthEast','LineWidth', 2)
    legend('boxoff')
    box off
    set(gcf,'color','white')
    set(gca,'tickdir','out','Layer','top', 'LineWidth', 1.6)
end

Trials = {'R-field/L-Ori','R-field/R-Ori','L-field/L-Ori','L-field/R-Ori'};

figure('Name','Comparison data and NL model trial triggered averages for each ROI','NumberTitle','off'); clf;
% suptitle('data and model comparison of trial triggered averages')
for j=1:length(allRois);
    for trialT=1:4
    subplot(8,4,(j-1)*4+trialT)
    axis([0 14 -1 1.2])
    plot(x,mean_Epochs(trialT,:,j),'-b','LineWidth', 1.6);
    hold on
    plot(x,mean_modnlEpochs(trialT,:,j),'-m','LineWidth', 1.6);
    plotErrorBars(x,mean_Epochs(trialT,:,j),SEM_Epochs(trialT,:,j),'b','b',1);
    hold on
    plotErrorBars(x,mean_modnlEpochs(trialT,:,j),SEM_modnlEpochs(trialT,:,j),'m','m',1);
    title(sprintf('Roi %s Trial %s',ROIS{j},Trials{trialT}),'FontSize',10)
    box off
    set(gcf,'color','white')
    set(gca,'tickdir','out','Layer','top', 'LineWidth', 1.6)
    end
end 
leg = legend('Data','Model');
set (leg,'Location', 'NorthEastOutside','LineWidth', 2)
legend('boxoff')

%Separate comparison

figure('Name','Comparison data and NL model trial triggered averages for each ROI','NumberTitle','off'); clf;
for j=1:length(allRois);
    subplot(4,2,j)
        plot(x,mean_Epochs(1,:,j),'-r','LineWidth', 1.6);
         hold on
        plot(x,mean_Epochs(2,:,j),'-m','LineWidth', 1.6);
        plot(x,mean_Epochs(3,:,j),'-b','LineWidth', 1.6);
        plot(x,mean_Epochs(4,:,j),'-c','LineWidth', 1.6);
        leg = legend('R-field/L-Ori','R-field/R-Ori','L-field/L-Ori','L-field/R-Ori');
        plotErrorBars(x,mean_Epochs(1,:,j),SEM_Epochs(1,:,j),'r','r',1);
        plotErrorBars(x,mean_Epochs(2,:,j),SEM_Epochs(2,:,j),'m','m',1);
        plotErrorBars(x,mean_Epochs(3,:,j),SEM_Epochs(3,:,j),'b','b',1);
        plotErrorBars(x,mean_Epochs(4,:,j),SEM_Epochs(4,:,j),'c','c',1);
        hold on 
        plot(x,mean_modnlEpochs(1,:,j),'--dr','LineWidth', 3);
        plot(x,mean_modnlEpochs(2,:,j),'--dm','LineWidth', 3);
        plot(x,mean_modnlEpochs(3,:,j),'--db','LineWidth', 3);
        plot(x,mean_modnlEpochs(4,:,j),'--dc','LineWidth', 3);
        title(sprintf('Roi %s',ROIS{j}),'FontSize',14);
        xlabel ('Time (TRs)','FontSize',12)
        ylabel ('Percent signal change (%)','FontSize',12)
        set (leg,'Location', 'NorthEast','LineWidth', 2)
        legend('boxoff')
        box off
        set(gcf,'color','white')
        set(gca,'tickdir','out','Layer','top', 'LineWidth', 1.6)
    
end 


%% ------- Orientation classification --------

% PART 1: Decoding hemifield (this is just to check our code does what it's
% supposed to).
% We will use 2 of our previously selected ROIs from the visual cortex (in
% proximity to the Calcarine Sulcus): one for the right visual cortex
% (VCRL), and one for the left visual cortex (VCLL) - both selected with a
% low threshold. 

VisROI{1} = roiVCRL;
VisROI{2} = roiVCLL;

% For each ROI, there will be two trial types (TYPE 1 = RIGHT HEMIFIELD,
% TYPE 2 + LEFT HEMIFIELD). 

% Create 4 empty matrices, 8(#runs) x # of voxels per ROI, for each
% roi/trial type combination:

r1B1 = nan(8,length(find(VisROI{1}))); % ROI 1 Trial Type 1
r1B2 = nan(8,length(find(VisROI{1}))); % ROI 1 Trial Type 2
r2B1 = nan(8,length(find(VisROI{2}))); % ROI 2 Trial Type 1
r2B2 = nan(8,length(find(VisROI{2}))); % ROI 2 Trial Type 2

nFrames = size(run1,4);
Drift = [1:nFrames]'; 

% We will begin by classifying the side of the stimulus presentation:

percentTseries_cell = cell(8,2); % number of functional runs and number of ROIs
detrend_percentTseries_cell = cell(8,2);

% Use regression to estimate the response amplitudes to each trial type,
% separately for each run and separately for each voxel in each ROI.

for j = 1:2;
    roi = VisROI{j};
    roiSize = length(find(roi));
    
    for i=1:length(allRuns)
        data = allRuns{i};
        numFrames = length(data);
        trialClass = allTTypes{i};
        tType = trialClass(:,1);
        
        tSeries = zeros(numFrames,roiSize);
        percentTseries = zeros(size(tSeries));
        [x,y,z] = ind2sub(size(roi),find(roi));
        
        for voxel = 1:roiSize
            tSeries(:,voxel) = squeeze(data(x(voxel),y(voxel),z(voxel),:));
            baseline = mean(tSeries(:,voxel));
            percentTseries(:,voxel) = 100 * (tSeries(:,voxel)/baseline - 1);
        end
 
        percentTseries_cell{i,j} = percentTseries;
        
        % first detrend using regression
        dBeta = Drift \ percentTseries; 
        detrend_percentTseries = percentTseries - repmat(dBeta,[nFrames,1]).* repmat(Drift,[1,roiSize]);
        detrend_percentTseries_cell{i,j} = detrend_percentTseries;
        
        % then run the analysis regression
        designMatrix = zeros(nFrames,2);
        
        designMatrix(:,1) = or(tType == 1, tType==2); %right hemifield
        designMatrix(:,2) = or(tType == 3, tType==4); %left hemifield
        designMatrix(:,1) = hrfconv(designMatrix(:,1),tau,delta,TR);
        designMatrix(:,2) = hrfconv(designMatrix(:,2),tau,delta,TR);
        
        modelInv = pinv(designMatrix);
        
        betas = modelInv * detrend_percentTseries;
        if j==1
            r1B1(i,:) = betas(1,:); % betas for trial type 1 (right hemifield)
            r1B2(i,:) = betas(2,:); % betas for trial type 2 (left hemifield)
        elseif j==2
            r2B1(i,:) = betas(1,:); % betas for trial type 1 (right hemifield)
            r2B2(i,:) = betas(2,:); % betas for trial type 2 (left hemifield)
        end
   
    end
end

% Use the response amplitudes from all but one run to train the classifier,
% and test the classifier using the remaining run. Repeat for each run in
% turn (leaving that run out of training).

group = cell(14,1);
group(1:7) = {'Right Hemifield'};
group(8:end) = {'Left Hemifield'};

GROUNDTRUTH = cell(2,1);
GROUNDTRUTH(1,1) = {'Right Hemifield'};
GROUNDTRUTH(2,1) = {'Left Hemifield'};
cp = classperf(GROUNDTRUTH);
classHemifield = num2cell(nan(2,8));

for lorun = 1:8 % left out sample run
    sample = [r1B1(lorun,:) r2B1(lorun,:) ; r1B2(lorun,:) r2B2(lorun,:)];
          
    % take out the sample run out of the training set:
    
    % Betas trial type 1
    trainR1B1 = r1B1;
    trainR1B1(lorun,:) = [];
    trainR2B1 = r2B1;
    trainR2B1(lorun,:) = [];
    % Betas trial type 2
    trainR1B2 = r1B2;
    trainR1B2(lorun,:) = [];
    trainR2B2 = r2B2;
    trainR2B2(lorun,:) = [];
    
    training = [trainR1B1 trainR2B1 ; trainR1B2 trainR2B2];

    class = classify(sample,training,group,'diaglinear');
    
    classHemifield{1,lorun} = class{1};
    classHemifield{2,lorun} = class{2};
    classperf(cp,class);
end
cp
fprintf('The classifier performance for hemifield is %3.1f%%\n',cp.CorrectRate*100);

%%
% PART 2: Decoding Orientation   
% For each ROI, there will be two trial types (TYPE 1 = RIGHT TILT,
% TYPE 2 + LEFT TILT). 

% Create 4 empty matrices, 8(#runs) x # of voxels per ROI, for each
% roi/trial type combination:

r1B1 = nan(8,length(find(VisROI{1}))); % ROI 1 Trial Type 1
r1B2 = nan(8,length(find(VisROI{1}))); % ROI 1 Trial Type 2
r2B1 = nan(8,length(find(VisROI{2}))); % ROI 2 Trial Type 1
r2B2 = nan(8,length(find(VisROI{2}))); % ROI 2 Trial Type 2


nFrames = size(run1,4);

for j = 1:2;
    roi = VisROI{j};
    roiSize = length(find(roi));
    
    for i=1:length(allRuns)
        data = allRuns{i};
        numFrames = length(data);
        trialClass = allTTypes{i};
        tType = trialClass(:,1);
        
        detrend_percentTseries = detrend_percentTseries_cell{i,j};
        
        % then run the analysis regression
        designMatrix = zeros(nFrames,2);
        
        designMatrix(:,1) = or(tType == 1, tType==3); % left orientation
        designMatrix(:,2) = or(tType == 2, tType==4); % right orientation
        designMatrix(:,1) = hrfconv(designMatrix(:,1),tau,delta,TR);
        designMatrix(:,2) = hrfconv(designMatrix(:,2),tau,delta,TR);
        
        modelInv = pinv(designMatrix);
        
        betas = modelInv * detrend_percentTseries;
        if j==1
            r1B1(i,:) = betas(1,:); % betas for left tilt
            r1B2(i,:) = betas(2,:); % betas for right tilt
        elseif j==2
            r2B1(i,:) = betas(1,:); % betas for left tilt
            r2B2(i,:) = betas(2,:); % betas for right tilt
        end
   
    end
end

group = cell(14,1);
group(1:7) = {'Right Tilt'};
group(8:end) = {'Left Tilt'};

GROUNDTRUTH = cell(2,1);
GROUNDTRUTH(1,1) = {'Right Tilt'};
GROUNDTRUTH(2,1) = {'Left Tilt'};
cp = classperf(GROUNDTRUTH);
classOrientation = num2cell(nan(2,8));

for lorun = 1:8 % left out sample run
    sample = [r1B1(lorun,:) r2B1(lorun,:) ; r1B2(lorun,:) r2B2(lorun,:)];
          
    % take out the sample run out of the training set:
    
    % Betas trial type 1
    trainR1B1 = r1B1;
    trainR1B1(lorun,:) = [];
    trainR2B1 = r2B1;
    trainR2B1(lorun,:) = [];
    % Betas trial type 2
    trainR1B2 = r1B2;
    trainR1B2(lorun,:) = [];
    trainR2B2 = r2B2;
    trainR2B2(lorun,:) = [];
    
    training = [trainR1B1 trainR2B1 ; trainR1B2 trainR2B2];

    class = classify(sample,training,group,'diaglinear');
    
    classOrientation{1,lorun} = class{1};
    classOrientation{2,lorun} = class{2};
    classperf(cp,class);
end
cp
fprintf('The classifier performance for orientation is %3.1f%%\n',cp.CorrectRate*100);

