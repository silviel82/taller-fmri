function plotErrorBars(x,y,ci,color1,color2,trans,ciOnly)

%-------------------------------------------
%
% plotErrorBars(x,y,ci,color1,color2,trans,ciOnly)
%
% plot some data with shaded error bars
%
% required inputs:
% x -- vector of x values
% y -- vector of y values
% ci -- vector of error bar widths at each point
%
% optional inputs:
% color1 -- color of error bars
% color2 -- color of data points
% trans -- transparency? (1=yes, 0=no)
% ciOnly -- just plot the error bars? (1=yes,0=no)
% outputs:
%
% freeman, 1-2-2010
%-------------------------------------------

% check arguments
if nargin < 3
    error('(plotErrorBars) not enough inputs');
end
if ~exist('color1','var') | isempty(color1)
    color1 = 'k';
end
if ~exist('color2','var') | isempty(color2)
    color2 = 'k';
end 
if ~exist('trans','var') | isempty(trans)
    trans = 0;
end
if ~exist('ciOnly','var') | isempty(ciOnly)
    ciOnly = 0;
end 

% vectorize inputs
x = x(:)';
y = y(:)';
ci = ci(:)';

% if given confidence interval is a scalar,
% apply the same value at all points
if length(ci) == 1
    ci = ones(size(x))*ci;
end

ind1=1:length(y);
ind2=ind1(end:-1:1);

% plot the error bars, with different settings
% if we want transparency
% (WARNING: transparency yields bitmapped .eps files)
if trans
    h=patch([x x(end:-1:1)],[y-ci y(ind2)+ci(ind2)],.6*ones(1,3), 'facecolor',color1,'edgecolor', 'none','facealpha',0.4);
else
    h=patch([x x(end:-1:1)],[y-ci y(ind2)+ci(ind2)],.6*ones(1,3), 'facecolor',color1,'edgecolor', 'none');
end   

% plot the data itself, unless flag is set
if ~ciOnly
hold on; 
plot(x,y,'-','linewidth',2,'color',color2,'MarkerEdgeColor',color2,'MarkerFaceColor',color2,'MarkerSize',3);
hold off
end
