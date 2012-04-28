function compare_canonical_id_output()
%
% Simulates the 12 models from the canonical identification and the Whipple
% and arm models for up to 20 second simulations.

addpath('../../src/matlab')

canonData = '/media/Data/Documents/School/UC Davis/Bicycle Mechanics/CanonicalBicycleID/data';
idmat = load([canonData filesep 'idMatrices.mat']);
load([canonData filesep 'goodRuns.mat']);

% only fit the minimal outputs of the Whipple model
states = {'phi', 'delta', 'phiDot', 'deltaDot'};
outputs = {'phi', 'delta', 'phiDot', 'deltaDot'};

C = eye(4);
K = zeros(4, 4);
x0 = zeros(4, 1);

% runs, models, outputs
fits = zeros(length(goodRuns), 14, length(outputs));

directory = '../statespaceid/exports';

for i = 1:length(goodRuns)
    runID = pad_with_zeros(num2str(goodRuns(i)), 5);
    display(sprintf('Fitting: %s', runID))

    tmp = load([directory filesep runID '.mat'], 'Maneuver');

    if strcmp(tmp.Maneuver(end-6:end), 'urbance')
        inputs = {'tDelta', 'fB'};
        latForce = true;
        display('Found Disturbance!')
        D = zeros(4, 2);
    else
        inputs = {'tDelta'};
        latForce = false;
        display('No disturbance.')
        D = zeros(4, 1);
    end

    % load the data
    [data, v, rider] = build_id_data([runID '.mat'], outputs, inputs, ...
        '../statespaceid/exports', true);

    models = fieldnames(idmat);

    % create state space model for all 14 models
    stateSpaceModels = cell(length(models) + 2, 1);

    for j = 1:length(models)
        [A, B] = canon_to_state_space(idmat.(models{j}), v);
        if latForce == false
            B = B(:, 1);
        end
        % I found an error here. ss threw an error if you pass the wrong
        % dimensions for D, but idss(A, B, C, D, K, x0) didn't and it would
        % set the x0 to be 8x1 instead of 4x1. I initially was passing a D
        % of 4 x 2 everytime and regardless of the number of inputs.
        %sys = ss(A, B, C, D);
        stateSpaceModels{j} = idss(A, B, C, D, K, x0, 0, 'StateName', states, ...
            'OutputName', outputs, 'InputName', inputs);
    end

    % add the whipple model
    stateSpaceModels{13} = bicycle_structured(['Rigid' rider], v, 'states', states, ...
        'inputs', inputs, 'outputs', outputs);

    % add the arm model
    arm = load(['armsAB-' rider '.mat']);
    indice = round(v * 10) + 1;
    A = squeeze(arm.stateMatrices(indice, :, :));
    B = squeeze(arm.inputMatrices(indice, :, [2, 3]));
    if latForce == false
        B = B(:, 1);
    end
    stateSpaceModels{14} = idss(A, B, C, D, K, x0, 0, 'StateName', states, ...
            'OutputName', outputs, 'InputName', inputs);

    % only select data 20 second a smaller parts of data to test against
    duration = data.SamplingInstants(end);
    limit = 8;
    if duration > limit
        [numSamples, ~] = size(data.OutputData);
        start = randi(numSamples - limit * 200, 1);
        stop = start + limit * 200;
    else
        start = 1;
        % the data(start:stop) was balking if I didn't have the round() for
        % some runs, what an annoying bug, matlab sucks
        stop = round(duration * 200);
    end

    [~, fit, ~] = compare(data(start:stop), stateSpaceModels{:});
    fits(i, :, :) = squeeze(fit);
    % save the plots of the fits
    fig = figure('Visible', 'off');
    figWidth = 8;
    figHeight = 16;
    set(gcf, ...
        'Color', [1, 1, 1], ...
        'PaperOrientation', 'portrait', ...
        'PaperUnits', 'inches', ...
        'PaperPositionMode', 'manual', ...
        'OuterPosition', [424, 305 - 50, 518, 465], ...
        'PaperPosition', [0, 0, figWidth, figHeight], ...
        'PaperSize', [figWidth, figHeight])
    compare(data(start:stop), stateSpaceModels{:})
    axes = findall(fig, 'type', 'axes');
    set(axes(8), 'ylim', [-pi / 12, pi / 12]) % phi
    set(axes(6), 'ylim', [-pi / 9, pi / 9]) % delta
    set(axes(4), 'ylim', [-0.5, 0.5]) % phiDot
    set(axes(2), 'ylim', [-1.5, 1.5]) % deltaDot
    saveas(fig, ['output-compare-plots' filesep runID '.png'])
    close all
end

% save all the data to a file
models = [models; 'Whipple'; 'Arm'];
save('output-compare.mat', 'goodRuns', 'fits', 'models')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A, B] = canon_to_state_space(idMatCell, v)
% function [A, B] = canon_to_state_space(idMatCell, v)
%
% Returns the state space form of the identified model.

M = [idMatCell{1, 1} idMatCell{1, 2}]';
C1 = [idMatCell{2, 1} idMatCell{2, 2}]';
K0 = [idMatCell{3, 1} idMatCell{3, 2}]';
K2 = [idMatCell{4, 1} idMatCell{4, 2}]';
H =  [idMatCell{5, 1}; idMatCell{5, 2}];
g = 9.81;

A = zeros(4); % [phi, delta, phiDot, deltaDot]
B = zeros(4, 2); % [Tdel, Fcl]

Minv = inv(M);

A(1:2, 3:4) = eye(2);
A(3:4, 1:2) = -Minv * (g * K0 + v^2 * K2);
A(3:4, 3:4) = -Minv * v * C1;

BT = Minv;
BF = BT * H;

B(3:4, 1) = BT(:, 2);
B(3:4, 2) = BF;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function num = pad_with_zeros(num, digits)
% Adds zeros to the front of a string needed to produce the number of
% digits.
%
% Parameters
% ----------
% num : string
%   A string representation of a number (i.e. '25')
% digits : integer
%   The total number of digits desired.
%
% If digits = 4 and num = '25' then the function returns '0025'.

for i = 1:digits-length(num)
    num = ['0' num];
end
