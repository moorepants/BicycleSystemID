function [data, percent] = remove_process_noise(data)
% function [data, percent] = remove_process_noise(data)
%
% Returns a data object in which the output is now the result of a
% simulation through an ARMAX derived transfer function.
%
% Parameters
% ----------
% data : iddata
%   A SISO data object.
%
% Returns
% -------
% data : iddata
%   A SISO data object with the same input as the provided data and a
%   "noise-free" output.
% percent : double
%   The percentage of the measured output variation explained by the
%   simulation output. Being that the simulation is always from zero initial
%   conditions, this will always be lower than if the initial conditions
%   were estimated.

model = armax(data, [8 8 4 0]);

s = warning('off', 'Control:analysis:LsimStartTime');
y = lsim(tf(model, 'm'), data.InputData, data.SamplingInstants);
warning(s)

percent = 1 - norm(y - data.OutputData) / norm(data.OutputData - ...
mean(data.OutputData));

data.OutputData = y;
