function [z, speed, rider] = build_id_data(runid, outputs)
% Returns a structure with an iddata object for run with PullForce as the
% input.
%
% Parameters
% ----------
% idRun : char
%   The file name of a run. (e.g. '00105.mat')
% outputs : cell array of chars
%   The desired list of outputs from the data.
%
% Returns
% -------
% z : iddata
%   The iddata for the run.
% speed : double
%   The mean speed for this run.

% load in the configuration variables
config
addpath(PATH_TO_CONTROL_MODEL)

% load all the variables into the runData structure
runData = load([PATH_TO_RUN_MAT_DIRECTORY filesep runid]);

dataInputs = {'PullForce'};
meijaardInputs = {'fB'};

meijaardOutputs = cell(size(outputs));
dataOutputs = cell(size(outputs));
for i = 1:length(outputs)
    dataOutputs{i} = convert_variable(outputs{i}, 'data');
    meijaardOutputs{i} = convert_variable(outputs{i}, 'meijaard');
end

% build the input matrix
u = zeros(length(runData.PullForce), length(dataInputs));
for j = 1:length(dataInputs)
    u(:, j) = runData.(dataInputs{j});
end

% build the output matrix
y = zeros(length(runData.RollRate), length(outputs));
for j = 1:length(outputs)
    y(:, j) = runData.(dataOutputs{j});
end

z = iddata(y, u, 1 / 200);

set(z, 'InputName', meijaardInputs)
set(z, 'InputUnit', get_units(dataInputs))
set(z, 'OutputName', meijaardOutputs)
set(z, 'OutputUnit', get_units(dataOutputs))

speed = mean(runData.ForwardSpeed);
rider = runData.Rider;

function units = get_units(signalNames)
% Returns the units for a given cell array of signal names.

unitMapping = struct('SteerAngle', 'Radian', ...
                     'RollAngle', 'Radian', ...
                     'YawAngle', 'Radian', ...
                     'LateralRearContact', 'Meter', ...
                     'SteerRate', 'Radian/Second', ...
                     'RollRate', 'Radian/Second', ...
                     'YawRate', 'Radian/Second', ...
                     'LateralRearContactRate', 'Meter/Second', ...
                     'PullForce', 'Newton', ...
                     'SteerTorque', 'Newton-Meter', ...
                     'RollTorque', 'Newton-Meter');

units = {};
for i = 1:length(signalNames)
    units{i} = unitMapping.(signalNames{i});
end
