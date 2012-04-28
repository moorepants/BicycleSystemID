function guess = guess_from_loop_closure(bicycle, rider, speed, calculation, varargin)
% function guess = guess_from_loop_closure(bicycle, rider, speed, calculation, varargin)
%
% Parameters
% ----------
% bicycle : char
%   The name of the bicycle, first letter capitalized and the rest
%   lowercase.
% rider : char
%   The name of the rider, first letter capitalized and the rest
%   lowercase.
% speed : double
%   The speed for the desired gain guess.
% calculation : char
%   This should be either 'full' or 'estimate', where 'full' finds the
%   optimal solution with the loop closure technique and 'estimate'
%   interpolates from previous found gain values. 'full' takes longer but
%   may be more accurate.
% gainGuess : double, 1 x 5
%   A starting guess for the gain search.
%
% Returns
% -------
% guess : vector, size(1, 6)
%   The gains and the neuromuscular frequency: [kDelta, kPhiDot, kPhi, kPsi,
%   kY, wnm]

config
addpath(PATH_TO_CONTROL_MODEL)

if size(varargin, 2) > 0
    extra = {'gainGuess', varargin{1};
else
    extrz = {};
end

bicycleRider = [bicycle rider];

wnm = 30;

if strcmp(calculation, 'full')
    data = generate_data(bicycleRider, speed, ...
                         'simulate', 0, ...
                         'loopTransfer', 0, ...
                         'handlingQuality', 0, ...
                         'forceTransfer', {}, ...
                         'fullSystem', 0, ...
                         'neuroFreq', wnm, ..
                         extra{:});
    guess = [data.modelPar.kDelta,
             data.modelPar.kPhiDot,
             data.modelPar.kPhi,
             data.modelPar.kPsi,
             data.modelPar.kY,
             sqrt(data.modelPar.neuroNum)]';
elseif strcmp(calculation, 'estimate')
    pathToGains = [PATH_TO_CONTROL_MODEL filesep 'gains' filesep ...
        bicycleRider 'SteerGains.txt'];
    [kDelta, kPhiDot, kPhi, kPsi, kY] = lookup_gains(pathToGains, speed);
    guess = [kDelta, kPhiDot, kPhi, kPsi, kY, wnm];
else
    err = MException('CalculationChk:IncorrectValue', ...
        '%s is not a valid option', calculation);
    throw(err)
end
