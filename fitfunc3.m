function [y] = fitfunc3(x,data,responseMatrix,TR)
%
% [y] = fitfunc(x,data,responseMatrix,TR)
%
% computes the expected mri signal by convolving
% response with a shifted gamma of width tau
% and delay delta and multiplying each column by an amplitude
% where amplitudes = x(1:4), tau=x(5), delta=x(6)
%
% y is the difference between the data
% and the expected mri signal
%
% responseMatrix is the expected neural response in TRs
% tau, delta and TR are in seconds
% see hrfconv for the convolution
%

t = [0 1 2 3]';
amp = x(1)+x(2)*exp(-x(3)*t);

tau = x(4);
delta = x(5);

designMatrix = zeros(length(data),4);
designMatrix(:,1) = hrfconv(responseMatrix(:,1),tau,delta,TR);
designMatrix(:,2) = hrfconv(responseMatrix(:,2),tau,delta,TR);
designMatrix(:,3) = hrfconv(responseMatrix(:,3),tau,delta,TR);
designMatrix(:,4) = hrfconv(responseMatrix(:,4),tau,delta,TR);

y = data - designMatrix*amp;

return