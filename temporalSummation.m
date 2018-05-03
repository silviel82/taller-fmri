%%%   Block Flickering Checkerboard  %%%
%
% This scripts creates and displays, for a brief moment (or two), the
% visual stimuli of Flickering Checkerboads BOTH on the right or on the
% left of the display
%

%%%%%%%%%%%%%%%%%%%%%

% Trial Layout
% tON: flicker on
% ISI: flicker off
% tON: flicker on if doubleFlag = 1
% until tTrial: flicker off
% doubleFlag controls whether one flash or two

% Some example setups follow

% Block design
%tON       = 5;
%ISI       = 0.25;
%tTrial    = 10;
%doubleFlag = 0;

% Flashes
%tON       = .125;
%ISI       = 0.25;
%tTrial    = 10;
% Single
%doubleFlag = 0;
% Double
%doubleFlag = 1;

%%%%%%%%%%%%%%%%%%%%%

%% Initialize variables
% scan parameters
TR = 0.75;                     % Repetition time (sec)
measurements = 280;          % Total number of measurements (= repetitions)

% timing characteristics
% stim paradigm
tON       = .125;               % time with stimulus ON (in units of TR)
ISI       = 0.25;               % inter-stimulus interval (when 2 presentations) (units of TR)
doubleFlag = 1;   % If you want the second flash (set to 0 if you want only 1 flash)
tTrial    = 28;                 % trial duration (units of TR)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% You should not need to modify below here %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nBlocks   =  measurements/tTrial; % number of blocks (per trial)

% a little error checking
if mod(nBlocks,1)
  sprintf('The number of measurements must be a multiple of tTrial')
  return;
end

%% stim characteristics
R         = 0.5;                % Radius of each checkerboard ring  (degrees)
Theta     = pi/16;              % Angle  of each checkerboard wedge (radians)
freq      = 16;                 % Frequency of flickering (in units of 1/TRs)
offAper   = 4;                  % Aperture off-center in 'x' (degrees)
dAper     = 6;                  % Aperture diameter          (degrees)
fixationW = 0.3;                % Fixation cross width       (degrees)
fixationLW= 1;                  % Fixation cross line width

% Screen setup
frameRate = 60;                 % (Hz)
screenDist= 57;                 % distance subject-screen (cm)
horAngle  = 16;                 % screen visual angle -horizontal (degrees)
horResol  = 1024;               % horizontal screen resolution (pixels)
verResol  = 768;                % vertical   screen resolution (pixels)
resolution= horResol/horAngle;  % screen resolution @ fixation (pixels/degree)
verAngle  = verResol/resolution;% screen visual angle -vertical (degrees)

framesPerBlock = frameRate * tTrial;

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
A = (A+1)/2;    % between 0 and 1;


%% Open screen & create texture
mglOpen;
mglVisualAngleCoordinates(screenDist,[horAngle verAngle]); % define coordinates
RChecker1 = mglCreateTexture(A*255);
RChecker0 = mglCreateTexture((1-A)*255);  % inverse


%% Create stencils
% Stencil 1: both right and left aperture
mglStencilCreateBegin(1);
% usage:   (        x,  y,[width height])
mglFillOval(  offAper,  0,[dAper  dAper]);
mglFillOval( -offAper,  0,[dAper  dAper]);
mglFixationCross(fixationW,fixationLW);
mglStencilCreateEnd;
mglClearScreen(0.5);


%% Wait for backtick from scanner:
% Show fixation while waiting:
mglClearScreen(0.5); 
mglFixationCross(fixationW,fixationLW);
mglFlush;
% Prepare first screen, so that we only need to flush
rotation = 180;   % rotate RCheckers so that the orientation is correct
mglClearScreen(0.5);
mglStencilSelect(1);
mglBltTexture(RChecker1,[0, 0], +1, 0, 180);  % positive contrast, Left
mglBltTexture(RChecker1,[0, 0], -1, 0,   0);  % positive contrast, Right
mglFixationCross(fixationW,fixationLW);
t1stON = zeros(nBlocks,freq*tON,2); % to store timing for 1st presentation
tISI   = zeros(nBlocks,1);          % to store timing for ISI
t2ndON = zeros(nBlocks,freq*tON,2); % to store timing for 2nd presentation
tBlank = zeros(nBlocks,1);          % to store timing for blank end of trial
t0 = mglGetSecs;                    % time when starting stimulus
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
  % First brief presentation:
  for cycle = 1:freq*tON  % alternation cycle number
    % check time, to see if we need to flush:
    while (1);
      t1 = mglGetSecs;            % get time
      if ( (t1-t0) >= TR/(2*freq) )
        t1stON(bn,cycle,1) = t1 - t0;  % store time since last 't0' (for debugging)
        t0 = t1;                  % update t0
        break;
      end;
    end
    % flush what we had prepared before:
    mglFlush;

    %%%%%
    % Prepare next screen:
    mglClearScreen(0.5);
    mglBltTexture(RChecker0,[0, 0], +1, 0, 180);  % positive contrast, Left
    mglBltTexture(RChecker0,[0, 0], -1, 0,   0);  % positive contrast, Right
    mglFixationCross(fixationW,fixationLW);

    % check time, to see if we need to flush:
    while (1);
      t1 = mglGetSecs;            % get time
      if ( (t1-t0) >= TR/(2*freq) )
        t1stON(bn,cycle,2) = t1 - t0;  % store time since last 't0' (for debugging)
        t0 = t1;                  % update t0
        break;
      end;
    end
    % flush what we had prepared before:
    mglFlush;
    %%%%%
    % Prepare next screen:
    mglClearScreen(0.5);
    mglBltTexture(RChecker1,[0, 0], +1, 0, 180);  % positive contrast, Left
    mglBltTexture(RChecker1,[0, 0], -1, 0,   0);  % positive contrast, Right
    mglFixationCross(fixationW,fixationLW);
  end
  
  % Wait Inter-Stimulus Interval (ISI) with just fixation
  mglClearScreen(0.5);
  mglFixationCross(fixationW,fixationLW);
  mglFlush;
  while (1);
    t1 = mglGetSecs;            % get time
    if ( (t1-t0) >= TR*ISI )
      tISI(bn) = t1 - t0;       % store time since last 't0' (for debugging)
      t0 = t1;                  % update t0
      break;
    end;
  end
  
  % Second brief presentation if required
  if ( doubleFlag==1 )
      % First brief presentation:
      for cycle = 1:freq*tON  % alternation cycle number
        % check time, to see if we need to flush:
        while (1);
          t1 = mglGetSecs;            % get time
          if ( (t1-t0) >= TR/(2*freq) )
            t1stON(bn,cycle,1) = t1 - t0;  % store time since last 't0' (for debugging)
            t0 = t1;                  % update t0
            break;
          end;
        end
        % flush what we had prepared before:
        mglFlush;

        %%%%%
        % Prepare next screen:
        mglClearScreen(0.5);
        mglBltTexture(RChecker0,[0, 0], +1, 0, 180);  % positive contrast, Left
        mglBltTexture(RChecker0,[0, 0], -1, 0,   0);  % positive contrast, Right
        mglFixationCross(fixationW,fixationLW);

        % check time, to see if we need to flush:
        while (1);
          t1 = mglGetSecs;            % get time
          if ( (t1-t0) >= TR/(2*freq) )
            t1stON(bn,cycle,2) = t1 - t0;  % store time since last 't0' (for debugging)
            t0 = t1;                  % update t0
            break;
          end;
        end
        % flush what we had prepared before:
        mglFlush;
        %%%%%
        % Prepare next screen:
        mglClearScreen(0.5);
        mglBltTexture(RChecker1,[0, 0], +1, 0, 180);  % positive contrast, Left
        mglBltTexture(RChecker1,[0, 0], -1, 0,   0);  % positive contrast, Right
        mglFixationCross(fixationW,fixationLW);
      end      
  end % end second presentation
    
  % Blank time (just fixation) until end of trial:
  mglClearScreen(0.5);
  mglFixationCross(fixationW,fixationLW);
  mglFlush;
  while (1);
    t1 = mglGetSecs;            % get time
    if ( (t1-t0) >= TR*(tTrial-tON-ISI-tON*doubleFlag) )
      tBlank(bn) = t1 - t0;     % store time since last 't0' (for debugging)
      t0 = t1;                  % update t0
      break;
    end;
  end
end
mglStencilSelect(0); % Stop using stencil
mglClearScreen(0.5); mglFlush

% total stim duration, to check that it all worked correctly.
% Since the first time recording included the waiting for backtick, we
% don't count it, and replace it by an extra TR/(2*freq)
totalStimT = TR/(2*freq) + sum(t1stON(2:end)) + sum(tISI(:)) + ...
  + sum(t2ndON(:)) + sum(tBlank(:))
% this totalStimT should be equal to TR*measurements

%% Close screen:
% mglPrivateClose;
mglClose;