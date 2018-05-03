%%%% roiAnalysis2.m
%%%% Silvia.11.29.15
% roiAnalysis2 is a function modification for function roiAnalysis.m that
% extracts the fMRI data from a single run based on a prespecified Region
% Of Interest (ROI). This new version drops first cycle as specified and
% computes the mean and sem of the entire time series across all voxels of
% the ROI then computes the average and sem across all trials
% (trial-triggered average).
% 
% INPUTS:   - epi:              the fMRI data from the chosen .nii file
%           - roi:              the mask from the selected .img file
%           - cycleDuration:    time duration of 1 cycle (in units of TR)
% 
% OUTPUT:   - newData:          mean of time series for each trial 
%  
function [newData] = roiAnalysis2(epi,roi,cycleDuration)
% Select the time series for the voxels defined within the ROI. This is a
% matrix with m rows of measurements (numFrames) and n columns of active
% voxels (roiSize). 
roiSize = length(find(roi));
numFrames = size(epi,4);
% create empty matrix to populate
tSeries = zeros(numFrames,roiSize);
% returns 3-D arrays containing the equivalent 3-D array subscripts
% equivalent to coordinates for the roi array of specified size
[x, y, z] = ind2sub(size(roi),find(roi));
% now we create the matrix by finding the timeseries corresponding to each
% of the identified voxels according to their 3-D coordinates
    for voxel = 1:roiSize
        tSeries(:,voxel) = squeeze(epi(x(voxel),y(voxel),z(voxel),:));
    end

% We will dump the first cycle, the entire run has 9 cycles and each cycle
% has 30 TRs (first 30 TRs = 22.5 seconds)
tSeries = tSeries(cycleDuration+1:end,:);

% We will reexpress the timeseries as percent signal change to correct for
% different baselines (mean levels of activity).
% Populate a new matrix:
percentTseries = zeros(size(tSeries));
% calculate the mean activity for each voxel and divide by the mean,
% express as percent change:
    for voxel = 1:roiSize
        baseline = mean(tSeries(:,voxel));
        percentTseries(:,voxel) = 100 * (tSeries(:,voxel)/baseline - 1);
    end

% Calculate the mean activity (as percent signal change) accross all voxels
data = mean(percentTseries,2);

% Calculate the average activity across all trials (trial-triggered average) 
% Compute number of cycles (trials) for this run.
% The ouput will be an mxn matrix with m trials and n time points

numCycles = length(data)/cycleDuration;
newData = zeros(numCycles,cycleDuration);
    for i = 1:numCycles;
        s = cycleDuration*(i-1)+1;
        newData(i,:) = data(s:(s+cycleDuration-1));
    end

end

