function [A, B, C, D, K, x0] = linear_parameterization(parameters, sampleTime, auxiliary)
% Returns the continous state space for the bicycle/rider closed loop control
% model.

sys = system_state_space(auxiliary.bicycle, parameters(1:5), ...
    parameters(6), auxiliary.inputs, auxiliary.outputs);

A = sys.A;
B = sys.B;
C = sys.C;
D = sys.D;

[numOutputs, numStates] = size(C);
K = zeros(numStates, numOutputs);
x0 = zeros(numStates, 1);
