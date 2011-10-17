function z = build_id_data(idRun, valRun)
% Returns a structure with an iddata object for an identification and
% validation run.
%
% Parameters
% ----------
% idRun : char
%   The file name of a run.
% valRun : char
%   The file name of a run.
%
% Returns
% -------
% z : structure
%   id : iddata
%       The iddata for the identification run.
%   val : iddata
%       The iddata for the validation run.

% load in the configuration variables
config

run.id = idRun;
run.val = valRun;

% load the two runs
dataType = {'id', 'val'};
for i = 1:length(dataType)
    runData.(dataType{i}) = load([PATH_TO_RUN_MAT_DIRECTORY run.(dataType{i})], ...
    'RollRate', 'SteerRate', 'RollAngle', 'SteerAngle', 'PullForce', ...
    'SteerTorque', 'YawRate', 'YawAngle', 'LateralRearContact', ...
    'LateralRearContactRate', 'NISampleRate');
end

inputs = {'PullForce'};

outputs = {'LateralRearContact', ...
           'YawAngle', ...
           'RollAngle', ...
           'SteerAngle', ...
           'LateralRearContactRate', ...
           'YawRate', ...
           'RollRate', ...
           'SteerRate', ...
           'SteerTorque'};

for i = 1:length(dataType)
    % build the input matrix
    u = zeros(length(runData.(dataType{i}).PullForce), length(inputs));
    for j = 1:length(inputs)
        u(:, j) = runData.(dataType{i}).(inputs{j});
    end

    % build the output matrix
    y = zeros(length(runData.(dataType{i}).RollRate), length(outputs));
    for j = 1:length(outputs)
        y(:, j) = runData.(dataType{i}).(outputs{j});
    end

    zTmp = iddata(y, u, 1 / 200);

    set(zTmp, 'InputName', inputs)
    set(zTmp, 'InputUnit', get_units(inputs))
    set(zTmp, 'OutputName', outputs)
    set(zTmp, 'OutputUnit', get_units(outputs))

    z.(dataType{i}) = zTmp;
end

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
                     'SteerTorque', 'Newton-Meter');

units = {};
for i = 1:length(signalNames)
    units{i} = unitMapping.(signalNames{i});
end
