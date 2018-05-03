% fmriTutorialPart1.m
% DJH 1/26/2004
% EV 08/05/2013 modified to use niftiread

%%%%%%%%%%%%%%%%%
% Preliminaries %
%%%%%%%%%%%%%%%%%

% This tutorial demonstrates how to load and view some fMRI data in matlab.
% In the process it provides some examples of matlab code that you can
% use to build your own more sophisticated analysis tools.

% To run the tutorial, you read the comments (like this one) and evaluate
% the code below it. Read the code as you go to see how it works. To
% evaluate the code, you can either: 1) copy and past it into the Matlab
% Command Window, 2) drag your mouse over a segment of code to select it
% and then choose "Evaluate Selection" from the Text menu, or 3) drag you
% mouse to select a code segment and then hit function key "f9" (for
% windows) or function key "f7" (for OSX).

% To begin, you need to download the Canned data set from the course
% website.

% Once you have downloaded the fmriTools and CannedData, you must define
% two pathnames to keep track of where they are. On a pc running windows,
% the pathnames will look something like this:
%
%     fmriToolsPath = 'c:\users\david\Matlab\fmriTools';
%     dataPath = 'c:\users\david\Data\CannedData_DS';
%
% On a mac running OSC the two pathnames will look something like this:
%
    fmriToolsPath = '/Users/silvia/Dropbox/MATLAB/fmriTools';
    dataPath = '/Users/silvia/Dropbox/Courses\:Readings/fMRI\ Lab/fmricourse/CannedData_DS';

%
% Edit the above  lines to have the correct paths and evaluate them.
% Then evaluate the following to add the fmriTools to your matlab path so
% that matlab will know how to find these functions.
addpath(fmriToolsPath);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loading and displaying MRI images %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Directory CannedData_DS has 5 subdirectories:
%
% 08+Bold_030314
% 09+Bold_030314
% 10+Bold_030314
% 11+t1_se_tra_concat_2
% 12+t1_mprage_sag
%
% Series 08,09,10 contain functional data. Series 11 contains the
% inplane anatomical images with the same slice prescription as the
% functional data. Series 12 contains high-resolution mprage anatomical
% images.

% Load the inplane images:
[inplanes, inplaneHeader] = niftiread(fullfile(dataPath,'11+t1_se_tra_concat_2',...
    'DS030975+11+t1_sra_concat_2.img'));

% The variable 'inplanes' now contains an array of image data. The size of
% the array is:
size(inplanes)
% [128 128 35]
% The first index is the x-coordinate. The second index is the
% y-coordinate. And the third index is the slice. To pull out the 10th
% slice, evaluate:
slice10 = inplanes(:,:,10);
size(slice10)
% If you are not familiar with the colon syntax in matlab, then you should
% learn about it. You'll be seeing quite a lot of it.

% Load the hires mprage anatomical images:
[mprage, mprageHeader] = niftiread(fullfile(dataPath,'12+t1_mprage_sag',...
    'DS030975+12+t1_mprage_sag.img'));

% Load the data from the first functional MRI scan (this might take a
% minute):
[bold1, bold1Header] = niftiread(fullfile(dataPath,'08+Bold_030314',...
    'DS030975+08+Bold_030314.img'));

% View the inplane images in sequence, pausing for a fraction of a second
% between presentations.
figure(1), clf
for slice = 1:size(inplanes,3)
    showImage(inplanes(:,:,slice)');
    pause(0.1)
end

% The function showImage automatically chooses how to rescale the images
% for display. You can also control the rescaling yourself, which comes in
% handy to help visualize certain aspects of the images more clearly.
% Remember that the images arrays are just filled with numbers. We display
% them by converting those numbers to grayscale intensities.

% Redisplay the inplanes clipping values below 200 and above 800.
figure(1), clf
clim = [200 800];
for slice = 1:size(inplanes,3)
    showImage(inplanes(:,:,slice)',clim);
    pause(0.1)
end

% Now display the hires mprage anatomical images:
figure(1), clf
clim = [50 800];
for slice = 1:size(mprage,3)
    showImage(mprage(:,:,slice)',clim);
    drawnow
end
% Those were axial (horizontal) slices (from bottom to top). Note that the
% drawnow function is needed to force matlab to show each slice. Otherwise
% it skips through them too quickly.

% Redisplay the hires mprage in saggital slices:
figure(1), clf
clim = [50 800];
for slice = 1:size(mprage,1)
	img = flipud(squeeze(mprage(slice,:,:))');
    showImage(img,clim);
    drawnow
end

% And coronal (from back to front):
figure(1), clf
clim = [50 800];
for slice = 1:size(mprage,2)
	img = flipud(squeeze(mprage(:,slice,:))');
    showImage(img,clim);
    drawnow
end

% The functional data are stored in a 4D array:
size(bold1)
% where the fourth dimension is time so that the indices correspond to 
% [x y z t]

% Display a movie (over time) of a single slice:
figure(1), clf
slice = 20;
for t = 1:size(bold1,4)
    showImage(bold1(:,:,slice,t)');
    drawnow
end
% Note that the image quality is pretty crummy compared to the anatomical
% images that you've been viewing up until now. That's because each a full
% volume of functional images was acquired every 2 sec whereas it took
% several minutes to acquire the anatomical images. The functional images
% are acquired in a different way to be sensitive to variations in
% bold oxygen level dependent (BOLD) signal over time.

% Plot the time series (image intensity over time) for an individual voxel.
tseries = squeeze(bold1(30,54,10,:));
figure(2)
plot(tseries);

% Try another one.
tseries = squeeze(bold1(32,12,10,:));
figure(2)
plot(tseries);

% This was a visual stimulation experiment in which a stimulus was
% displayed in the right hemifield for 10 seconds and then in the left
% hemifield for the next 10 seconds. So you should a modulation in time
% series at locations in visual cortex that are stimulated by either the
% left or right halves of the visual field.

% The function 'plotTseries' lets you do this interactively. Take a look at
% the source code for this function. Because the inplanes are registered
% with the functional images, we display an inplane slice to click on and
% then extract the corresponding time series from the functionals. Try to
% find some interesting voxels. The function 'ginput' prompts you to select
% a pixel by clicking on the figure. Click left to select a point. Click
% any other button to quit.
slice = 10;
TR = 2;
plotTseries(inplanes,bold1,slice,TR);
