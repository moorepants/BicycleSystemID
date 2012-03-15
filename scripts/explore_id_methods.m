% This script is just a way to get a better feel of what the system
% identification toolbox is actually doing.
PATH_TO_BICYCLE_SYSTEM_ID = '/media/Data/Documents/School/UC Davis/Bicycle Mechanics/BicycleSystemID';
addpath(PATH_TO_BICYCLE_SYSTEM_ID)
dataDir = [PATH_TO_BICYCLE_SYSTEM_ID filesep 'scripts' filesep 'exports'];

inputs = {'tDelta', 'fB'};
%inputs = {'tDelta'};
outputs = {'phi', 'delta', 'phiDot', 'deltaDot'};
states = {'phi', 'delta', 'phiDot', 'deltaDot'};

% 273 is a run to test with no pertubations, 280 has perturbations
[data, speed, rider] = build_id_data('00588.mat', outputs, inputs, dataDir, true);

bicycle = bicycle_structured(['Rigid' rider], speed, 'inputs', inputs, ...
    'outputs', outputs, 'states', states);

first = pem(data, bicycle, 'Display', 'on');

pemArgs = { ...
           'Focus', 'Prediction', ...
           'DisturbanceModel', 'Estimate', ...
           'InitialState', 'Estimate', ...
           'SearchMethod', 'lsqnonlin', ...
           'Maxiter', 1000, ...
           'LimitError', 0.0, ...
           'Tolerance', 0.005, ...
           'Criterion', 'trace', ...
           'Weighting', [0.9, 0, 0, 0;
                         0, 1, 0, 0;
                         0, 0, 1, 0;
                         0, 0, 0, 1], ...
           'Display', 'on', ...
          };

second = pem(data, first, pemArgs{:});

% defautlPem, prediction and simulation all seem to give the same result
%defaultPem = pem(data, bicycle, 'InitialState', 'Estimate', 'DisturbanceModel', 'Estimate', 'Display', 'on');
%defaultPem = pem(data, bicycle, 'Display', 'on');
%prediction = pem(data, bicycle, 'Focus', 'Prediction', 'Display', 'on')
%simulation = pem(data, bicycle, 'Focus', 'Simulation', 'Display', 'on')
%passband = pem(data, bicycle, 'Focus', [0.1, 10.])
% stability crashed
%stability = pem(data, bicycle, 'Focus', 'Stability', 'Display', 'on')
% focus can also to a frequency range

%figure(1)
%compare(data, defaultPem, prediction, simulation, passband)

%init = pem(data, bicycle, 'InitialState', 'Estimate', 'Display', 'on')
%noise = pem(data, bicycle, 'DisturbanceModel', 'Estimate', 'Display', 'on')
%initNoise = pem(data, bicycle, 'InitialState', 'Estimate', 'DisturbanceModel', 'Estimate', 'Display', 'on')
%
%figure(2)
%compare(data, init, noise, initNoise)

% the initNoiseSim model gave unstable output
%initNoiseSim = pem(data, bicycle, ...
    %'InitialState', 'Estimate', ...
    %'DisturbanceModel', 'Estimate', ...
    %'Focus', 'Simulation', ...
    %'Display', 'on')
% this one worked, the deltaDot was great but the rest were about 17% fit
%initNoisePassband = pem(data, bicycle, ...
    %'InitialState', 'Estimate', ...
    %'DisturbanceModel', 'Estimate', ...
    %'Focus', [0.1, 10.], ...
    %'Display', 'on')
%
%figure(3)
%compare(data, initNoiseSim, initNoisePassband)

% search methods
% the default search method that the default 'auto' selects seems to be
% 'lsqnonlin'. I thought that it would use one of the other solutions
% because the minimization function can be formed analytically as a least
% squares problem.

% only giving the search method option: gn and lm seem to give the same
% answer in terms of the A, B matrices, but all fit percentages are slightly
% different among the four methods.
%searchMethods = {'gn', 'gna', 'lm', 'lsqnonlin'};
%for i = 1:length(searchMethods)
    %smModels.(searchMethods{i}) = pem(data, bicycle, ...
        %'SearchMethod', searchMethods{i}, ...
        %'InitialState', 'Estimate', ...
        %'DisturbanceModel', 'Estimate', ...
        %'Display', 'on', ...
        %'MaxIter', 80)
%end
