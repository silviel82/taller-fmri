function showImage(image,clim)
% function showImage(image,[cmin cmax])
%
% image: just what you think it is
% [cmin cmax] determines how the image is rescaled for display (see
% imagesc)
%
% DJH 1/26/2004

if ~exist('clim','var')
    histThresh = length(image(:))/1000;
    [cnt, val] = hist(image(:),100);
    goodVals = find(cnt>histThresh);
    cmin = val(min(goodVals));
    cmax = val(max(goodVals));
    clim = [cmin cmax];
end

imagesc(image,clim);
colormap('gray')
axis image, axis off
