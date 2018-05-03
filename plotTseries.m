function plotTseries(inplanes,bold,slice,TR)
% function plotTseries(inplanes,bold,slice,TR)
%
% Plots time series for selected voxels. Displays the selected inplane
% slice in figure 1. Left click on the image plots a graph (in figure 2) of
% the time series from the corresponding bold voxel. Right click quits.
%
% inplanes: inplane anatomical images
% bold: fMRI data assumed to be in register with the inplanes
% slice: slice that you want to display and plot
% TR: repetition time in the fMRI acquisition (needed to get the axis
% labeled correctly in seconds)
%
% DJH 1/26/2004

% compute duration in seconds for use later in making the plot
duration = TR*size(bold,4);
timepoints = [TR:TR:duration];

% display the inplane
figure(1); clf;
showImage(inplanes(:,:,slice)');

% loop while left button is clicked
cont = 1;
while cont
    
    % Choose figure 1 and get mouse click
    figure(1)
    [x,y,button]=ginput(1);
    if (button ~= 1)
        cont = 0;
    end
    
    % If left button was pressed...
    if cont
        
        % Divide x and y coordinates by 2 because inplane resolution is
        % twice that of the functional images
        x = round(x/2);
        y = round(y/2);
        
        % Extract the tseries
        tseries = squeeze(bold(x,y,slice,:));
        
        % Plot time series, adding title and axis labels
        figure(2), clf;
        plot(timepoints,tseries);
        xlim([0 duration]);
        title(['Time series for voxel [x y z]: ',...
                num2str(x),' ',num2str(y),' ',num2str(slice)]);
        ylabel('Image intensity');
        xlabel('Time (sec)');
        
    end
end






