function [zees, ems] = bicycle_pem(idRun, checkRun)

pathToData = ['..' filesep 'BicycleDataProcessor' filesep 'exports' ...
    filesep 'mat' filesep];

dataType = {'id', 'check'};
for i = 1:length(dataType)
    runData.(dataType{i}) = load([pathToData idRun], 'RollRate', 'SteerRate', 'RollAngle', ...
        'SteerAngle', 'PullForce', 'SteerTorque', 'YawRate', 'YawAngle', ...
        'LateralRearContact', 'LateralRearContactRate');
end

idNumSamples = length(runData.id.RollRate);

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

% this the only input when looking at the whole system
inputs.system = {'PullForce'};

% these are the inputs when looking at only the bicycle
inputs.bicycle = {'SteerTorque', 'PullForce'};

outputs.measuredRates = {'SteerRate', 'RollRate', 'YawRate'};
outputs.allRates = {'SteerRate', 'RollRate', 'YawRate', ...
    'LateralRearContactRate'};
outputs.ratesCoords = {'SteerRate', 'RollRate', 'YawRate', ...
    'LateralRearContactRate', 'SteerAngle', 'RollAngle', 'YawAngle', ...
    'LateralRearContact'};

sampleSize = 1 / 200;

for ins = fieldnames(inputs)'
    for outs = fieldnames(outputs)'

        % build the input matrix
        inList = inputs.(char(ins));
        u = zeros(idNumSamples, length(inList));
        for i = 1:length(inList)
            u(:, i) = runData.id.(inList{i});
        end

        % build the output matrix
        outList = outputs.(char(outs));
        y = zeros(idNumSamples, length(outList));
        for i = 1:length(outList)
            y(:, i) = runData.id.(outList{i});
        end

        ztmp = iddata(y, u, sampleSize);

        set(ztmp, 'InputName', inList)
        set(ztmp, 'InputUnit', get_units(inList, unitMapping))
        set(ztmp, 'OutputName', outList)
        set(ztmp, 'OutputUnit', get_units(outList, unitMapping))

        zees.([char(ins) '_' char(outs)]) = ztmp;
    end
end

for model = fieldnames(zees)'
    z = zees.(char(model));
    ze = z(1850:17450);
    display(sprintf('Finding the model for %s', char(model)))
    ems.(char(model)) = pem(ze);
end

%bode(m)
%
%compare(ze, m)

function units = get_units(signalNames, unitMapping)

units = {};
for i = 1:length(signalNames)
    units{i} = unitMapping.(signalNames{i});
end
