% fmriTutorialPart3.m
% DJH 2/11/2004
% EV  8/5/2013 modified to use niftiread

%%%%%%%%%%%%%%%%%
% Preliminaries %
%%%%%%%%%%%%%%%%%

% To begin, you need to download the Canned data set from the course
% website.

% Once you have downloaded the fmriTools and CannedData, you must define
% two pathnames to keep track of where they are. On a pc running windows,
% the pathnames will look something like this:
%
%     fmriToolsPath = 'c:\users\david\Matlab\fmriTools';
%     dataPath = 'c:\users\david\Data\CannedData_DS';

%
% On a mac running OSX the two pathnames will look something like this:
%
    fmriToolsPath = '/Users/silvia/Dropbox/MATLAB/fmriTools';
    dataPath = '/Users/silvia/Dropbox/MATLAB/fmriTools/CannedData_DS';
%
% Edit the above  lines to have the correct paths and evaluate them.
% Then evaluate the following to add the fmriTools to your matlab path so
% that matlab will know how to find these functions.
% addpath(fmriToolsPath);

%%%%%%%%%%%%%%%%%
% Load the data %
%%%%%%%%%%%%%%%%%

% This data directory contains the following files:
%
% DS030975+06+Bold_030314-header.txt	
% DS030975+06+Bold_030314.hdr		
% DS030975+06+Bold_030314.img		
% ROImask.hdr
% ROImask.img
%
% The BOLD fMRI data are in DS030975+06+Bold_030314.img
%
% The file my-mask.img contains a mask that was saved out from fsl,
% indicating the pixels that were active in the fMRI expeiment.

% Load the fMRI data:
[epi, epiHeader] = niftiread(fullfile(dataPath,'06+Bold_030314','DS030975+06+Bold_030314.img'));

% Load the mask:
[roi, roiHeader] = niftiread(fullfile(dataPath,'06+Bold_030314','ROImask.img'));

% View the mask images in sequence.
figure(1), clf
for slice = 1:size(roi,3)
    showImage(roi(:,:,slice)',[0 1]);
    title(['slice ',num2str(slice)]);
    drawnow
    pause(0.4)
end
% Which slices contain the ROI (region of interest)?

% Remember that it's always a good idea to look at your data.
% Display a movie (over time) of a single slice from the bold images:
figure(1), clf
slice = 15;
for t = 1:size(epi,4)
    showImage(epi(:,:,slice,t)');
    drawnow
end

% How many voxels are in the ROI? Note that you should type
% "help find" to learn what this matlab function does.
roiSize = length(find(roi))

% How long (number of temporal frames) are the time series?
numFrames = size(epi,4)

% Pick out the time series for those voxels. We will produce a matrix of
% time series data in which each column corresponds to the time series from
% one of the active voxels. The size of the matrix is nTimePoints by
% nVoxels. Note that you should type "help ind2sub" to learn what this
% matlab function does.
tSeries = zeros(numFrames,roiSize);
[x y z] = ind2sub(size(roi),find(roi));
for voxel = 1:roiSize
    tSeries(:,voxel) = squeeze(epi(x(voxel),y(voxel),z(voxel),:));
end

% Plot the time series from all of the voxels, superimposed on top of one
% another:
figure(2), clf
plot(tSeries)
ylabel('fMRI signal (raw image intensity');
xlabel('Frame');
title('Time series of all voxels in the ROI');

% Notice that they all have a different baseline/mean intensity. Let's fix
% that.
percentTseries = zeros(size(tSeries));
for voxel = 1:roiSize
    baseline = mean(tSeries(:,voxel));
    percentTseries(:,voxel) = 100 * (tSeries(:,voxel)/baseline - 1);
end
figure(2), clf
plot(percentTseries)
ylabel('fMRI response (% change in image intensity');
xlabel('Frame');
title('Time series of all voxels in the ROI');

% Plot the mean percent time series. Type "help mean" to learn how to take
% the mean across the rows of the matrix (i.e., mean across voxels instead
% of the mean across time).
data = mean(percentTseries,2);
plot(data);
ylabel('fMRI response (% change in image intensity');
xlabel('Frame');
title('Mean across voxels in the ROI');
 

%Simulate activity for the model
t = [1:131];
blockDuration = 6;
neuralActivity = mod(ceil(t/blockDuration),2);

tau = 2;
delta = 2;

% Plot the HIRF with these parameter values:
t = [0:1:30];
tshift = max(t-delta,0);
HIRF = (tshift/tau).^2 .* exp(-tshift/tau) / (2*tau);

modelActivity = conv(neuralActivity,HIRF);
modelActivity = modelActivity(1:length(neuralActivity));
nTime = length(neuralActivity);
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
model;

% Now estimate b (the beta weights). The equation is:
%    y = X b
% where we want to solve for b. So we compute:
%    b = pinv(X) * y
% where the 'pinv' function in matlab computes the pseudo-inverse of a
% matrix.
modelInv = pinv(model);
b = modelInv * data


% Note that you can also do this without the loop:
%   b = modelInv * data;

% % Plot the parameter estimates
% figure(2); clf;
% subplot(3,1,1)
% plot(b(1,:))
% title('Estimates of Neural Activity')
% ylabel('Amplitude of Neural Activity (arb units)')
% xlabel('Position (voxel #)')
% subplot(3,1,2)
% plot(b(2,:))
% title('Estimates of Drift')
% ylabel('Drift rate (delta image intensity/sec)')
% xlabel('Position (voxel #)')
% subplot(3,1,3)
% plot(b(3,:))
% title('Estimates of Baseline Intensity')
% ylabel('Mean image intensity (raw image intensity units)')
% xlabel('Position (voxel #)')

