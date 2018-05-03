function [conv_stim_vec] = hrfconv(stim_vec,tau,delta,TR)
%
% [conv_stim_vec] = hrfconv(stim_vec,tau,delta,TR)
%
% convolves stim_vec with a shifted gamma of width tau
% and delay delta
% stiv_vec is the expected neural response in TRs
% tau, delta and TR are in seconds
%

% Plot the HIRF with these parameter values:
t = [0:TR:30];
tshift = max(t-delta,0);
HIRF = TR*(tshift/tau).^2 .* exp(-tshift/tau) / (2*tau);

temp = conv(stim_vec,HIRF);

conv_stim_vec = temp(1:length(stim_vec));

return