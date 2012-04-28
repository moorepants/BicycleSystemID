function model = bicycle_grey(stype, bicycle, speed, outputs, gains, neuro)
% function model = bicycle_grey(stype, bicycle, speed, outputs, gains, neuro)
%
% Returns grey box model for the Hess control model.
%
% Parameters
% ----------
% stype : char
%   Either 'heading' or 'lateral' to specify the desired system to return.
% bicycle : char or ss
%   Name of the bicycle or a bicycle state space system.
% speed : double
%   The bicycle speed.
% outputs : cell array of chars
%   The desired outputs of the model in Meijaard form.
% gains : double, 1 x 4 or 1 x 5
%   Either give an empty matrix [] for the gains to be calculated or provide
%   the gains as [kDelta, kPhiDot, kPhi, kPsi, kY] for 'lateral' and
%   [kDelta, kPhiDot, kPhi, kPsi] for 'heading'. These will be used as the
%   initial parameter guesses.
% neuro : double, 1 x 2
%   The neuromuscular natural frequency and damping ratio, [wnm, zetanm]. If
%   an empty matrix is supplied the default values will be set to 30 and
%   0.707. This will be used as an initial parameter guess.
%
% Returns
% -------
% model : idgrey
%   A grey box model of the bicycle/rider system with a single input,
%   lateral force, and the outputs specified in the arguments. The initial
%   guess for the gains and neuromuscular frequencies are the one's
%   predicted from the Hess method unless specified in the arguments.

% bring in the control model functions
config
addpath(PATH_TO_CONTROL_MODEL)

% If a bicycle state space is supplied, set the bicycle to RigidJason for
% the possible gain calculations.
if ischar(bicycle)
    bicycleName = bicycle;
else
    bicycleName = 'RigidJason';
end

if isempty(neuro)
    neuro = [30.0, 0.707];
end

if isempty(gains)
    % calculate the bicycle state space and the gains
    hess = generate_data(bicycleName, speed, ...
                         'gains', gains, ...
                         'neuroFreq', neuro(1), ...
                         'simulate', 0, ...
                         'loopTransfer', 0, ...
                         'handlingQuality', 0, ...
                         'forceTransfer', {}, ...
                         'display', 0);
    % build the grey box model
    guess = [hess.modelPar.kDelta,
             hess.modelPar.kPhiDot,
             hess.modelPar.kPhi,
             hess.modelPar.kPsi,
             hess.modelPar.kY,
             neuro(1),
             neuro(2)];
    if strcmp(stype, 'heading')
        guess(5) = [];
    end
    % get the state space model so we don't have to calculate it during the gain
    % matching algorithm
    if ischar(bicycle)
        aux.bicycle = ss(hess.modelPar.A, hess.modelPar.B, hess.modelPar.C, ...
            hess.modelPar.D);
        aux.bicycle.StateName = hess.bicycle.states;
        aux.bicycle.OutputName = hess.bicycle.outputs;
        aux.bicycle.InputName = hess.bicycle.inputs;
    else
        aux.bicycle = bicycle;
    end
else
    guess = [gains, neuro];
    if ischar(bicycle)
        par = par_text_to_struct([PATH_TO_CONTROL_MODEL filesep 'parameters' ...
            filesep bicycle 'Par.txt']);
        aux.bicycle = whipple_pull_force_abcd(par, speed);
    else
        aux.bicycle = bicycle;
    end
end

% set the desired outputs and inputs for the resulting model and data set
aux.outputs = outputs;
aux.inputs = {'fB'};
aux.stype = stype;

model = idgrey('linear_parameterization', guess, 'c', aux);
model.InputName = aux.inputs;
model.OutputName = aux.outputs;
model.StateName = [aux.bicycle.StateName; 'tDelta'; 'tDeltaDot'];
PName = {'kDelta', 'kPhiDot', 'kPhi', 'kPsi', 'kYQ', 'wnm', 'zetanm'};
if strcmp(stype, 'heading')
    PName(5) = [];
end
model.PName = PName;
