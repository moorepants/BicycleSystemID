function [z, meta] = build_id_data(runid, outputs, inputs, directory, varargin)
% Returns an iddata object for the run.
%
% Parameters
% ----------
% idRun : char
%   The file name of a run. (e.g. '00105.mat')
% outputs : cell array of chars
%   The desired list of outputs from the data.
% inputs : cell array of chars
%   The desired list of inputs from the data.
% directory : char
%   The path to a directory containing the data files. Specificy '' if you
%   want to use the default directory in config.m.
% detrend : boolean, optional
%   If true, the data will be detrended (mean subracted). The PullForce will
%   not be detrended.
%
% Returns
% -------
% z : iddata
%   The iddata for the run.
% meta : structure
%   speed : double
%       The mean speed for this run.
%   rider : char
%       The rider of the bicycle for this run.
%   maneuver : char
%       The maneuver name.

% load in the configuration variables
config
addpath(PATH_TO_CONTROL_MODEL)

if ~isempty(directory)
    PATH_TO_RUN_MAT_DIRECTORY = directory;
end

% load all the variables into the runData structure
runData = load([PATH_TO_RUN_MAT_DIRECTORY filesep runid]);

meijaardInputs = cell(size(inputs));
dataInputs = cell(size(inputs));
for i = 1:length(inputs)
    dataInputs{i} = convert_variable(inputs{i}, 'data');
    meijaardInputs{i} = convert_variable(inputs{i}, 'meijaard');
end

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

% the sample rate has been 200 so far, but that should theorectically be
% loaded from runData also
z = iddata(y, u, 1 / 200);

if size(varargin, 2) > 0
    if varargin{1} == true
        trend = getTrend(z, 0);
        % if one of the inputs is the lateral force, don't detrend it, as it
        % may not be a very symmetric signal
        if ismember('PullForce', dataInputs)
            trend.InputOffset(find(strcmp('PullForce', dataInputs))) = 0.0;
        end
        z = detrend(z, trend);
    end
end

set(z, 'InputName', meijaardInputs)
set(z, 'InputUnit', get_units(dataInputs))
set(z, 'OutputName', meijaardOutputs)
set(z, 'OutputUnit', get_units(dataOutputs))

meta.speed = mean(runData.ForwardSpeed);
meta.rider = runData.Rider;
meta.maneuver = runData.Maneuver;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sub Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function units = get_units(signalNames)
% Returns the units for a given cell array of signal names.

unitMapping = struct('SteerAngle', 'Radian', ...
                     'RollAngle', 'Radian', ...
                     'YawAngle', 'Radian', ...
                     'LateralRearContact', 'Meter', ...
                     'LateralFrontContact', 'Meter', ...
                     'SteerRate', 'Radian/Second', ...
                     'RollRate', 'Radian/Second', ...
                     'YawRate', 'Radian/Second', ...
                     'LateralRearContactRate', 'Meter/Second', ...
                     'LateralFrontContactRate', 'Meter/Second', ...
                     'PullForce', 'Newton', ...
                     'SteerTorque', 'Newton-Meter', ...
                     'RollTorque', 'Newton-Meter');

units = {};
for i = 1:length(signalNames)
    units{i} = unitMapping.(signalNames{i});
end
