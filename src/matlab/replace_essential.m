function bicycle = replace_essential(bicycle, M, C1, K0, K2, H, v, g)
% function bicycle = replace_essential(bicycle, M, C1, K0, K2, H, v, g)
%
% Returns a bicycle state space in which the essential equations of motion
% have been replaced with the ones derived from the provided canonical form.
%
% Parameters
% ----------
% bicycle : ss
%   A bicycle state space with at least the states [phi delta phiDot
%   delatDot] and inputs [tPhi tDelta fB].
% M : double, 2 x 2
%   The mass matrix.
% C1 : double, 2 x 2
%   The damping matrix proportional to speed.
% K0 : double, 2 x 2
%   The stiffness matrix proportional to gravity.
% K2 : double, 2 x 2
%   The stiffness matrix proportional to the square of speed.
% H : double, 2 x 1
%   The vector which gives the propotion of the lateral force to the roll
%   and steer torques.
% v : double
%   The speed at which the bicycle state space was evaluated at.
% g : double
%   The acceleration due to gravity.
%
% Returns
% -------
% bicycle : ss
%   The augmented bicycle state space.

row = [find(strcmp(bicycle.StateName, 'phiDot')), ...
       find(strcmp(bicycle.StateName, 'deltaDot'))];
aCol = [find(strcmp(bicycle.StateName, 'phi')), ...
        find(strcmp(bicycle.StateName, 'delta')), ...
        find(strcmp(bicycle.StateName, 'phiDot')), ...
        find(strcmp(bicycle.StateName, 'deltaDot'))];
bCol = [find(strcmp(bicycle.InputName, 'tPhi')), ...
        find(strcmp(bicycle.InputName, 'tDelta')), ...
        find(strcmp(bicycle.InputName, 'fB'))];

invM = inv(M);
bicycle.A(row, aCol) = [-invM * [K0 * g + K2 * v^2], -invM * C1 * v];
bicycle.B(row, bCol) = [invM, invM * H];
