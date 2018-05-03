function [y] = fitfunc(x,data,stim_vec,TR)
%
% [y] = fitfunc(x,data,stim_vec,TR)
%
% computes the expected mri signal by convolving
% stim_vec with a shifted gamma of width tau
% and delay delta and multiplying by an amplitude
% where amplitude = x(1), tau=x(2), delta=x(3)
%
% y is the difference between the data
% and the expected mri signal
%
% stim_vec is the expected neural response in TRs
% tau, delta and TR are in seconds
% see hrfconv for the convolution
%

amp = x(1);
tau = x(2);
delta = x(3);

conv_stim_vec = hrfconv(stim_vec,tau,delta,TR);

y = data - amp*conv_stim_vec;

return