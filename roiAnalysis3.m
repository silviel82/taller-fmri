%%%% roiAnalysis3.m
%%%% Silvia.12.06.15
% roiAnalysis3 is a function modification for function roiAnalysis2.m that
% extracts the fMRI data from a single run based on a prespecified Region
% Of Interest (ROI). 
% 
% INPUTS:   - epi:              the fMRI data from the chosen .nii file
%           - roi:              the mask from the selected .img file
%           - runNumber:        number of fMRI run for figure
%           - onDuration:       stimulus on duration (in units of full TRs)
%           - offDuration:      stimulus off duration (in units of full TRs)
%           - numberOfCycles:   number of cycles per run
% 
% OUTPUT:   - figure(1):    mean time series from data with shaded error
%                           bar
%           - estMean:      the estimated response amplitude (% signal
%                           change)averaged across all voxels in the ROI.

%  
function meanTS = roiAnalysis3(epi,roi)
% Select the time series for the voxels defined within the ROI. This is a
% matrix with m rows of measurements (numFrames) and n columns of active
% voxels (roiSize). 
roiSize = length(find(roi));
numFrames = size(epi,4);
% create empty matrix to populate
tSeries = zeros(numFrames,roiSize);
% returns 3-D arrays containing the equivalent 3-D array subscripts
% equivalent to coordinates for the roi array of size 80 64 35
[x, y, z] = ind2sub(size(roi),find(roi));
% now we create the matrix by finding the timeseries corresponding to each
% of the identified voxels according to their 3-D coordinates
for voxel = 1:roiSize
    tSeries(:,voxel) = squeeze(epi(x(voxel),y(voxel),z(voxel),:));
end

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
meanTS = mean(percentTseries,2);
% Calculate the sem around the mean activity accross all voxels
error = std(percentTseries,1,2)/sqrt(length(percentTseries(1,:)));

% % plot mean data with shaded error bar
% figure('Name',sprintf('Mean Response and Error Run:%g\n',runNumber),'NumberTitle','off'); clf;
% plot(data)
% hold on
% shadedErrorBar([1:240],data,error,'-b',1)
% ylabel('fMRI response (% change in image intensity');
% xlabel('Time Course (TRs)');
% title('Mean activity across voxels in the ROI');

end