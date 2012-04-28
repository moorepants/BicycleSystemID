function [parameters, fit, speed] = id_rider_siso(runFile, output, ...
    gainGuess, neuroGuess, removeNoise, ignoreHeading, bicycle, varargin)
% function [parameters, fit, speed] = id_rider_siso(runFile, output, ...
%     gainGuess, neuroGuess, removeNoise, ignoreHeading, varargin)
%
% Parameters
% ----------
% runFile : char
%   The complete path to the run file.
% output : char
%   The single output to fit.
% gainGuess : double, 1 x 5 or 1 x 4
%   Guesses for the five or four gains [kDelta, kPhiDot, kPhi, kPsi, kYq].
% neuroGuess : double, 1 x 2
%   Guesses for the two neuromuscular parameters [wnm, zetanm].
% removeNoise : boolean
%   If true the process and output noise will be removed with an ARMAX
%   model.
% ignoreHeading : boolean
%   If true the "Balance With Disturbance" runs will include the lateral
%   deviation loop.
% bicycle : cell array
%   If true the essential bicycle dynamics will be replaced with the
%   provided canonical form. An empty cell array leaves the Whipple model
%   definition otherwise specify {M, C1, K0, K2, H}.
% varargin : char value pairs
%   Any remaining arguments will be passed to the pem function.
%
% Returns
% -------
% parameters : double, 1 x 7 or 1 x 6
%   The optimal gains and neuromuscular parameters.
% fit : double
%   The percent of the output variation explained by the model. This is with
%   respect to the noise-free output if removeNoise is true.
% speed : double
%   The mean speed during the run.

[directory, fileName, ext] = fileparts(runFile);
[data, meta] = build_id_data([fileName ext], {output}, {'fB'}, directory, true);

if removeNoise
    [data, percent] = remove_process_noise(data);
end

if strcmp(meta.maneuver, 'Balance With Disturbance')
    if ignoreHeading
        stype = 'lateral';
    else
        stype = 'heading';
        gainGuess = gainGuess(1:4);
    end
elseif strcmp(meta.maneuver, 'Track Straight Line With Disturbance')
    stype = 'lateral';
else
    err = MException('InvalidRun:NotDisturbance', ...
        'The run must be one with lateral disturbances');
    throw(err)
end

if isempty(bicycle)
    grey = bicycle_grey(stype, ['Rigid' meta.rider], meta.speed, {output}, ...
        gainGuess, neuroGuess);
else
    whipple = bicycle_state_space(['Rigid' meta.rider], meta.speed);
    augmented = replace_essential(whipple, bicycle{:}, meta.speed, 9.81);
    grey = bicycle_grey(stype, augmented, meta.speed, {output}, ...
        gainGuess, neuroGuess);
end

optimal = pem(data, grey, varargin{:});
parameters = optimal.ParameterVector;

[~, fit, ~] = compare(data, optimal);

speed = meta.speed;
