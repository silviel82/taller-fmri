%%%   Flickering Checkerboard LR  %%%
%
% This scripts creates and displays the visual stimuli of Flickering
% Checkerboads on the right or on the left of the display (we start with
% the left)`
%

%% Initialize variables
% scan parameters
TR = 2;                     % Repetition time (sec)
measurements = 65;          % Total number of measurements (= repetitions)
% stim characteristics
R         = 0.5;                % Radius of each checkerboard ring  (degrees)
Theta     = pi/16;              % Angle  of each checkerboard wedge (radians)
freq      = 8;                  % Frequency of flickering (in units of 1/TRs)
offAper   = 4;                  % Aperture off-center in 'x' (degrees)
dAper     = 6;                  % Aperture diameter          (degrees)
fixationW = 0.3;                % Fixation cross width       (degrees)
fixationLW= 1;                  % Fixation cross line width
contrast = 90;                 % Percent contrast
% stim paradigm
tBlock    =  5;                   % block duration (units of TRs)
nBlocks   =  measurements/tBlock; % number of blocks (per trial)
if mod(nBlocks,1)
  sprintf('The number of measurements must be a multiple of tBlock')
  return;
end
% setup
frameRate = 60;                 % (Hz)
screenDist= 57;                 % distance subject-screen (cm)
horAngle  = 16;                 % screen visual angle -horizontal (degrees)
horResol  = 1024;               % horizontal screen resolution (pixels)
verResol  = 768;                % vertical   screen resolution (pixels)
resolution= horResol/horAngle;  % screen resolution @ fixation (pixels/degree)
verAngle  = verResol/resolution;% screen visual angle -vertical (degrees)

framesPerBlock = frameRate * tBlock;


%% Create checkerboard (only half)
% In vertical, we need to create only the diameter of the aperture
% In horizontal, we need to create the aperture offset + 0.5 * diameter
% If you want to create the whole (half) screen:
% A = zeros(verResol,horResol/2);
A = zeros(dAper*resolution,(offAper+dAper/2)*resolution);
i0 = round(size(A,1)/2+1); j0 = 1;
for j=1:size(A,2)       % column index
  for i=1:size(A,1)     % row index
    if (j-j0)==0
      A(i,j) = (-1)^floor( sqrt((i-i0).^2+(j-j0).^2)./(R*resolution) );
    else
      A(i,j) = (-1)^floor( sqrt((i-i0).^2+(j-j0).^2)./(R*resolution) ) .* ...  % radius
        (-1)^floor( atan( (i-i0)/(j-j0) )./Theta );               % angle
    end
  end
end
A = contrast/100*A/2; % between -C/2:+C/2

%% Open screen & create texture
mglOpen;
mglVisualAngleCoordinates(screenDist,[horAngle verAngle]); % define coordinates
RChecker1 = mglCreateTexture(255*(0.5+A));
RChecker0 = mglCreateTexture(255*(0.5-A));  % reverse


%% Create stencils
% Stencil 1: left aperture
mglStencilCreateBegin(1);
% usage:   (        x,  y,[width height])
mglFillOval( -offAper,  0,[dAper  dAper]);
mglFixationCross(fixationW,fixationLW);
mglStencilCreateEnd;
mglClearScreen;
% Stencil 2: right aperture
mglStencilCreateBegin(2);
% usage:   (        x,  y,[width height])
mglFillOval(  offAper,  0,[dAper  dAper]);
mglFixationCross(fixationW,fixationLW);
mglStencilCreateEnd;
mglClearScreen;


%% Wait for backtick from scanner:
% Show fixation while waiting:
mglClearScreen(0.5); 
mglFixationCross(fixationW,fixationLW);
mglFlush;
% Prepare first screen (we draw on the left), so that we only need to flush
hAlignment = +1;
rotation = 180;   % rotate RCheckers so that the orientation is correct
mglClearScreen(0.5);
mglStencilSelect(1);    % left aperture
mglBltTexture(RChecker1,[0, 0], hAlignment, 0, rotation);  % positive RChecker
mglFixationCross(fixationW,fixationLW);
t  = zeros(nBlocks,tBlock*freq,2);    % to store timing
t0 = mglGetSecs;                      % time when starting stimulus
ti = t0;
% Wait for keypress, and check whether it is a backtick:
while (1);
  backtickCode = mglCharToKeycode({'`'});
  k = mglGetKeys( backtickCode );
  if (any(k))
    break;
  end;
end


%% Draw stim:
for bn = 1:nBlocks      % block number
  if mod(bn,2)
    % we draw on the left
    mglStencilSelect(1);    % left aperture
    hAlignment = +1;
    rotation = 180;   % rotate RCheckers so that the orientation is correct
  else
    mglStencilSelect(2);    % right aperture
    hAlignment = -1;
    rotation = 0;     % no rotation
  end

  for cycle = 1:freq*tBlock  % alternation cycle number
    % check time, to see if we need to flush:
    while (1);
      t1 = mglGetSecs;            % get time
      if ( (t1-t0) >= TR/(2*freq) )
        t(bn,cycle,1) = t1 - t0;  % store time since last 't0' (for debugging)
        t0 = t1;                  % update t0
        break;
      end;
    end
    % flush what we had prepared before:
    mglFlush;

    %%%%%
    % Prepare next screen:
    mglClearScreen;
    mglBltTexture(RChecker0,[0, 0], hAlignment, 0, rotation);  % negative RChecker
    mglFixationCross(fixationW,fixationLW);

    % check time, to see if we need to flush:
    while (1);
      t1 = mglGetSecs;            % get time
      if ( (t1-t0) >= TR/(2*freq) )
        t(bn,cycle,2) = t1 - t0;  % store time since last 't0' (for debugging)
        t0 = t1;                  % update t0
        break;
      end;
    end
    % flush what we had prepared before:
    mglFlush;
    %%%%%
    % Prepare next screen:
    mglClearScreen;
    mglBltTexture(RChecker1,[0, 0], hAlignment, 0, rotation);  % positive RChecker
    mglFixationCross(fixationW,fixationLW);
  end
end
mglStencilSelect(0); % Stop using stencil
mglClearScreen; mglFlush

% total stim duration, to check that it all worked correctly.
% Since the first time recording included the waiting for backtick, we
% don't count it, and replace it by an extra TR/(2*freq)
totalStimT = TR/(2*freq) + sum(t(2:end))
% this totalStimT should be equal to TR*measurements

%% Close screen:
mglPrivateClose;
% mglClose;