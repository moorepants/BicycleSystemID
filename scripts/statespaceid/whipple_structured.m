function whipple_structured(directory, varargin)
% function whipple_structured(directory, varargin)
%
% Identifies the parameters of a 4th order structured black box state space
% model of a bicycle.
%
% Parameters
% ----------
% directory : char
%   Specify the directory which contains the .mat files to estimate from.
% estimateX0 : boolean, optional
%   If true, the initial states will be estimated.
% estimateK : boolean, optional
%   If true, the disturbance model (Kalman gain matrix), will be estimated.
% onesFree : boolean, optional
%   If true all non-zero parameters in the A matrix will be estimated.
% phiWeight : double, optional
%   The desired weighting relative to 1 of the roll angle signal.

addpath('..')

% parse all the optional inputs
p = inputParser;
p.addRequired('directory');
p.addParamValue('estimateX0', false);
p.addParamValue('estimateK', false);
p.addParamValue('onesFree', false);
p.addParamValue('lightenPhi', 1);
p.addParamValue('new', false);
p.parse(directory, varargin{:});
args = p.Results;

pemArgs = {'Maxiter', 100};
pemArgs = [pemArgs {'Weighting', [args.lightenPhi, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1]}];

tag = '';

if args.estimateX0 == true
    pemArgs = [pemArgs {'InitialState', 'Estimate'}];
    tag = [tag '-x0'];
end

if args.estimateK == true
    pemArgs = [pemArgs {'Disturbance', 'Estimate'}];
    tag = [tag '-K'];
end

if args.onesFree == true
    tag = [tag '-freeones'];
end

if args.lightenPhi ~= 1
    pemArgs = [pemArgs {'Criterion', 'Trace'}];
    tag = [tag '-phi'];
end

if args.new == true
    pemArgs = { ...
               'Focus', 'Prediction', ...
               'DisturbanceModel', 'Estimate', ...
               'InitialState', 'Estimate', ...
               'SearchMethod', 'lsqnonlin', ...
               'Maxiter', 1000, ...
               'LimitError', 0.0, ...
               'Tolerance', 0.005, ...
               'Criterion', 'trace', ...
               'Weighting', [0.95, 0, 0, 0;
                             0, 1, 0, 0;
                             0, 0, 1, 0;
                             0, 0, 0, 1], ...
               'Display', 'on', ...
              };
    tag = [tag '-new'];
end

w = what(args.directory);
matFiles = sort(w.mat);

% only fit the minimal states of the Whipple model
states = {'phi', 'delta', 'phiDot', 'deltaDot'};
outputs = {'phi', 'delta', 'phiDot', 'deltaDot'};

% initialize matrices for all the data we want to save
stateMatrices = zeros(length(matFiles), length(states), length(states));
inputMatrices = zeros(length(matFiles), length(states), 2);
fits = zeros(length(matFiles), length(outputs));
initialConditions = zeros(length(matFiles), length(states));
speeds = zeros(length(matFiles), 1);
durations = zeros(length(matFiles), 1);

for i = 1:length(matFiles)
    display(sprintf('Fitting: %s', matFiles{i}))
    tmp = load([args.directory filesep matFiles{i}], 'Maneuver');

    if strcmp(tmp.Maneuver(end-6:end), 'urbance')
        inputs = {'tDelta', 'fB'};
        display('Found Disturbance!')
    else
        inputs = {'tDelta'};
        display('No disturbance.')
    end

    % load the data
    [data, v, rider] = build_id_data(matFiles{i}, outputs, inputs, ...
        args.directory, true);

    speeds(i) = v;

    durations(i) = data.SamplingInstants(end);

    % build a structured continuous idss model
    whippleModel = bicycle_structured(['Rigid' rider], v, 'states', states, ...
        'inputs', inputs, 'outputs', outputs);

    if args.onesFree == true
        whippleModel.As(1, 3) = NaN;
        whippleModel.As(2, 4) = NaN;
    end

    % identify the free entries of the structured model
    if args.new == true
        first = pem(data, whippleModel, 'Display', 'on');
        identifiedModel = pem(data, first, pemArgs{:});
    else
        identifiedModel = pem(data, whippleModel, pemArgs{:});
    end

    % store the resulting state space matrices and initial conditions
    stateMatrices(i, :, :) = identifiedModel.A;
    if strcmp(tmp.Maneuver(end-6:end), 'urbance')
        inputMatrices(i, :, :) = identifiedModel.B;
    else
        inputMatrices(i, :, 1) = identifiedModel.B;
    end
    initialConditions(i, :) = identifiedModel.x0;

    % save the plots of the fits
    [~, fit, ~] = compare(data, identifiedModel);
    fits(i, :) = fit(:);
    fig = figure('Visible', 'off');
    compare(data, identifiedModel)
    saveas(fig, ['compare-plots' filesep matFiles{i}(1:end-4) tag '.png'])
    close all
end

% save all the data to a file
save(['whipple-structured-results' tag '.mat'], 'matFiles', 'stateMatrices', ...
    'inputMatrices', 'fits', 'initialConditions', 'speeds', 'durations')
