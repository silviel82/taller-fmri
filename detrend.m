function detrended_vec = detrend(vec)
%
% function detrended_vec = detrend(vec)
% Removes linear drift
%

npts = length(vec);
X = [linspace(-1,1,npts)' ones(npts,1)];
b = pinv(X)*vec;
detrended_vec = vec - X*b;

return

