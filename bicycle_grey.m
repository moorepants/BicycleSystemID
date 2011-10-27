function model = bicycle_grey(bicycle, speed, outputs, gains, neuroFreq)
% function [data, model] = bicycle_grey(runFile, bicycle, speed, outputs, gains, neuroFreq)
% Returns the data for a run and a compatible grey box model for the Hess
% control model.
%
% Parameters
% ----------
% bicycle : char
%   Name of the bicycle.
% speed : double
%   The bicycle speed.
% outputs : cell array of chars
%   The desired outputs of the model in Meijaard form.
% gains : matrix, 5 x 1
%   Either give and empty matrix [] for the gains to be calculated or
%   provide the gains as [kDelta, kPhiDot, kPhi, kPsi, kY].
% neuroFreq : double
%   The default value should be 30, otherwise specify the desired value.
%
% Returns
% -------
% model : idgrey
%   A grey box model of the bicycle/rider system with a single input,
%   lateral force, and the outputs specified in the arguments. The initial
%   guess for the gains and neuromuscular frequencies are the one's
%   predicted from the Hess method otherwise they are as specified in the
%   arguments.

% bring in the control model functions
config
addpath(PATH_TO_CONTROL_MODEL)

% calculate the bicycle state space and the gains (unless specified)
hess = generate_data(bicycle, speed, ...
                     'gains', gains, ...
                     'neuroFreq', neuroFreq, ...
                     'simulate', 0, ...
                     'loopTransfer', 0, ...
                     'handlingQuality', 0, ...
                     'forceTransfer', {}, ...
                     'display', 0);

% get the state space model so we don't have to calculate it during the gain
% matching algorithm
aux.bicycle = hess.bicycle; % loads the state, input and output names
aux.bicycle.A = hess.modelPar.A;
aux.bicycle.B = hess.modelPar.B;
aux.bicycle.C = hess.modelPar.C;
aux.bicycle.D = hess.modelPar.D;

% set the desired outputs and inputs for the resulting model and data set
aux.outputs = outputs;
aux.inputs = {'fB'};

% build the grey box model
guess = [hess.modelPar.kDelta,
         hess.modelPar.kPhiDot,
         hess.modelPar.kPhi,
         hess.modelPar.kPsi,
         hess.modelPar.kY,
         neuroFreq];

model = idgrey('linear_parameterization', guess, 'c', aux);
model.InputName = aux.inputs;
model.OutputName = aux.outputs;
model.StateName = [aux.bicycle.states, 'tDelta', 'tDeltaDot'];
