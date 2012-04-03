function grey = bicycle_structured(bicycle, speed, varargin)
% function grey = bicycle_structured(bicycle, speed, varargin)
%
% Returns a structured idss model of the Whipple model linearized about the
% nominal configuration.
%
% Parameters
% ----------
% bicycle : char
%   The name of the bicycle.
% speed : double
%   The speed of the bicycle.
% varargin : char/cell array pairs, optional
%   Specify a subset of states, inputs or outputs by setting one of the
%   following: `states`, `inputs`, `outputs` as a cell array of
%   chars which include the subset variable names. Beaware that not all
%   state, input and output combinations are necessarily possible.
%
%   Valid state names: 'xP', 'yP', 'psi', 'phi', 'thetaP', 'thetaR', 'delta',
%       'thetaF', 'phiDot', 'thetaRDot', 'deltaDot'
%   Valid input names: 'tPhi', 'tDelta', 'fB'
%   Valid output names: 'xP', 'yP', 'psi', 'phi', 'thetaP', 'thetaR', 'delta',
%       'thetaF', 'xPDot', 'yPDot', 'psiDot', 'phiDot', 'thetaPDot',
%       'thetaRDot', 'deltaDot', 'thetaFDot', 'xQ', 'yQ'
%
% Returns
% -------
% grey : idss
%   The structured idss model.

% bring in the control model functions
config
addpath(PATH_TO_CONTROL_MODEL)

% The following definition of the structured model could potentially set a
% an entry in the A, B, C, D matrices that should be free to structured
% because in only checks for zeros and ones. There is the possibily that
% with a certain combination of parameters one of the entries would equal
% zero or one. It would be smarter to explicity define the nan for every
% entry in the model based on the known results. This is the quicker way,
% that should work 99% of the time (or more).

% generate the bicycle model based on the Whipple model for the base model
bicycleSS = bicycle_state_space(bicycle, speed, varargin{:});

grey = idss(bicycleSS, 'SSParameterization', 'Structured');

% for every entry in A, B, C, D that is not zero or 1, set the structured to
% nan
for m = {'A', 'B', 'C', 'D'}
    mat = bicycleSS.(m{:});
    matS = bicycleSS.(m{:});
    % finds all non zero entries
    nonZero = abs(mat) > 1e-10;
    % finds all non one entries
    nonOne = abs(abs(mat) - ones(size(mat))) > 1e-10;
    matS(nonZero & nonOne) = nan;
    grey.([m{:} 'S']) = matS;
end
