function [A, B, C, D, K, x0] = linear_parameterization(parameters, sampleTime, auxiliary)
% Returns the continous state space for the bicycle/rider closed loop control
% model.

data = generate_data(auxiliary.bicycle, auxiliary.speed, ...
                     'gains', parameters(1:5), ...
                     'neuroFreq', parameters(6), ...
                     'stateSpace', auxiliary.stateSpace, ...
                     'simulate', 0, ... don't simulate
                     'loopTransfer', 0, ... don't find the loop tfs
                     'handlingQuality', 0, ... don't find the HQM
                     'forceTransfer', {}, ... don't find the force tfs
                     'display', 0);

% the state space for the continous model
A = data.system.A;
B = data.system.B;
C = data.system.C;
D = data.system.D;
outputs = data.system.outputs;

keepers = {'yP', ...
           'psi', ...
           'phi', ...
           'delta', ...
           'psiDot', ...
           'phiDot', ...
           'deltaDot', ...
           'tDelta'}';

% find the outputs that are not in keepers
[~, I] = setdiff(outputs, keepers);
% remove the rows in the outputs that aren't in keepers
C(I, :) = [];
D(I) = [];

K = zeros(17, 8);
x0 = zeros(17, 1);
