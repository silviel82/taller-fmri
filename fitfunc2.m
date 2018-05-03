function [y] = fitfunc2(x,data,responseMatrix,TR)
%
% [y] = fitfunc(x,data,responseMatrix,TR)
%
% computes the expected mri signal by convolving
% response with a shifted gamma of width tau
% and delay delta and multiplying each columng by an amplitude
% where amplitudes = x(1:4), tau=x(5), delta=x(6)
%
% y is the difference between the data
% and the expected mri signal
%
% stim_vec is the expected neural response in TRs
% tau, delta and TR are in seconds
% see hrfconv for the convolution
%

tau = x(5);
delta = x(6);

designMatrix = zeros(length(data),4);
designMatrix(:,1) = hrfconv(responseMatrix(:,1),tau,delta,TR);
designMatrix(:,2) = hrfconv(responseMatrix(:,2),tau,delta,TR);
designMatrix(:,3) = hrfconv(responseMatrix(:,3),tau,delta,TR);
designMatrix(:,4) = hrfconv(responseMatrix(:,4),tau,delta,TR);

y = data - designMatrix*x(1:4)';

return