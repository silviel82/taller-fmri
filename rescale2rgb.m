function rgb = rescale2rgb(image,cmap,clip)
%function rgb = rescale2rgb(image,cmap,[clipMin,clipMax])

if ~exist('clip','var')
    % Choose clipping based on histogram
    histThresh = length(image(:))/1000;
    [cnt, val] = hist(image(:),100);
    goodVals = find(cnt>histThresh);
    clipMin = val(min(goodVals));
    clipMax = val(max(goodVals));
    clip = [clipMin,clipMax];
else
    clipMin = clip(1);
    clipMax = clip(2);
end

% Clip
result = image;
result(find(image < clipMin)) = clipMin;
result(find(image > clipMax)) = clipMax;

% Scale
indices = round(255 * (result-clipMin)/(clipMax-clipMin)) + 1;
indices = max(1,min(indices,size(cmap,1)));

% Extract r,g,b components
r = zeros(size(image));
g = zeros(size(image));
b = zeros(size(image));
r(:) = cmap(indices,1);
g(:) = cmap(indices,2);
b(:) = cmap(indices,3);

% Stuff them into rgb
dims = [size(image),3];
rgb = cat(length(dims),r,g,b);