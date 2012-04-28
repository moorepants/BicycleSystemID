function [A, B, C, D, K, x0] = linear_parameterization(parameters, sampleTime, aux)
% function [A, B, C, D, K, x0] = linear_parameterization(parameters, sampleTime, aux)
%
% Returns the continous state space for the bicycle/rider closed loop
% control model.
%
% Parameters
% ----------
% parameters : double, 1 x 6 or 1 x 7
%   The closed loop free parameters: [kDelta kPhiDot kPhi kPsi <kYq> wnm
%   zetanm].
% sampleTime : double
%   The required sample time argument. This does not effect the continous
%   definition.
% aux : structure
%   stype : char
%       Either 'heading' or 'lateral'.
%   bicycle : ss
%       The basic bicycle state space system with all the state, input and
%       output names defined in the Meijaard form.
%   inputs : cell array of char
%       The desired inputs for the closed loop system.
%   outputs : cell array of char
%       The desired outputs for the closed loop system.
%
% Returns
% -------
% A : double
%   The state matrix.
% B : double
%   The input matrix.
% C : double
%   The output matrix.
% D : double
%   The feed forward matrix.
% K : double
%   The Kalman gain matrix.
% x0 : double
%   The intial conditions.

if strcmp(aux.stype, 'lateral') && length(parameters) == 7
    gains = parameters(1:5);
    neuro = parameters(6:7);
elseif strcmp(aux.stype, 'heading') && length(parameters) == 6
    gains = parameters(1:4);
    neuro = parameters(5:6);
else
    exception = MException('VerifyParameters:LengthMismatch', ...
       'The length of the parameters does not match the system type.');
    throw(exception);
end

sys = system_state_space(aux.stype, aux.bicycle, gains, neuro, ...
    aux.inputs, aux.outputs);

A = sys.A;
B = sys.B;
C = sys.C;
D = sys.D;

[numOutputs, numStates] = size(C);

K = zeros(numStates, numOutputs);

x0 = zeros(numStates, 1);
