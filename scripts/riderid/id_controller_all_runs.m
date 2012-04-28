addpath('../../src/matlab')

% Here I get a list of all the "good" disturbance runs. These have been
% loaded from disturbance_runs.py.
disturbanceRunDir = '../../data/riderid/disturbance-runs';
contents = what(disturbanceRunDir);
matFiles = contents.mat;

% Load in the default guesses for the eight speed bins.
speedBin = importdata('../../data/riderid/speed-bin-guesses.txt');

% This is the model identified from Luke's pavilion runs.
M = [129.3615, 2.5592;
     2.5592, 0.2505];
C1 = [ 0,   33.5263;
      -0.5486,    2.0997];
K0 = [-115.7074, -4.5261;
      -4.5261,   -0.4889];
K2 = [0, 103.9425;
      0, 2.6034];
H = [0.9017;
     0.0111];

zetanm = 0.707;
output = 'delta';

if exist('../../data/riderid/bestControllerIdResults.mat') == 2
    goodResults = load('../../data/riderid/bestControllerIdResults.mat');
    parameters = goodResults.parameters;
    fits = goodResults.fits;
    speeds = goodResults.speeds;
else
    parameters = zeros(length(matFiles), 7);
    fits = zeros(length(matFiles), 1);
    speeds = zeros(length(matFiles), 1);
    save('../../data/riderid/bestControllerIdResults.mat', ...
        'parameters', 'fits', 'speeds', 'matFiles')
end

%subset = matFiles(143:end);
subset = matFiles;
for i = 1:length(subset)

    % Load the speed and metadata.
    [~, meta] = build_id_data(subset{i}, {'delta'}, {'fB'}, ...
        disturbanceRunDir);
    speed = meta.speed;

    display(repmat('=', 1, 79))
    display(sprintf('Run file: %s %s at %1.2f m/s.', subset{i}, meta.maneuver, speed))
    display(repmat('=', 1, 79))

    % Load in the previously found results.
    goodResults = load('../../data/riderid/bestControllerIdResults.mat');

    if any(isnan(goodResults.parameters(strcmp(matFiles, subset{i}), :))) || abs(sum(goodResults.parameters(strcmp(matFiles, subset{i}), :)) - 0) < 1e-10
        % First compute the best fit using the speed bin seed values.
        speedBins = speedBin.data(:, 1);
        [~, speedBinIndex] = min(abs(speedBins - speed));
        speedBinGains = speedBin.data(speedBinIndex, 2:6);
        speedBinNeuro = [speedBin.data(speedBinIndex, 7), zetanm];
        display(sprintf('Trying gains: %1.1f, %1.2f, %1.2f, %1.2f, %1.2f and neuro: %1.1f, %1.3f from speed bin %1.2f m/s.', speedBinGains, speedBinNeuro, speedBin.data(speedBinIndex, 1)))

        try
            [par, fit, ~] = id_rider_siso([disturbanceRunDir filesep ...
                subset{i}], output, speedBinGains, speedBinNeuro, false, true, ...
                {M, C1, K0, K2, H}, 'FixedParameter', 'zetanm', 'Focus', 'stability');
            display(sprintf('The fit from the speed bin guess is %1.0f %% with parameters %1.4f, %1.4f, %1.4f, %1.4f, %1.4f, %1.4f, %1.4f', fit, par))
        catch err
            display(err.identifier)
            fit = nan;
            par = ones(1, 7) * nan;
        end
    else
        display('Loading the best result from previous computations.')
        fit = goodResults.fits(strcmp(matFiles, subset{i}));
        par = goodResults.parameters(strcmp(matFiles, subset{i}), :);
        display(sprintf('The loaded fit is %1.0f %% with parameters %1.4f, %1.4f, %1.4f, %1.4f, %1.4f, %1.4f, %1.4f', fit, par))
    end

    if fit < 70 || isnan(fit)
        display('The speed bin fit is less that 70 %, so more guesses will be tried')
        % Select some guesses from similar speeds in previous good results.
        bool = (goodResults.speeds < (speed + 0.5)) & ...
               (goodResults.speeds > (speed - 0.5));
        trialIndices = find(bool);
        notTrialIndices = find(~bool);

        if isempty(trialIndices)
            trialGuesses = [];
            notTrialGuesses = goodResults.parameters;
        else
            trialGuesses = goodResults.parameters(trialIndices, :);
            notTrialGuesses = goodResults.parameters(notTrialIndices, :);
        end

        if size(trialGuesses, 1) > 5
            r = randperm(size(trialGuesses, 1));
            trialGuesses = trialGuesses(r(1:5), :);
        end

        r = randperm(size(notTrialGuesses, 1));
        trialGuesses = [trialGuesses; notTrialGuesses(r(1:3), :)];

        % Compute the SISO fit for all of the guesses.
        trialFits = zeros(size(trialGuesses, 1), 1);
        trialParameters = zeros(size(trialGuesses, 1), 7);

        for j = 1:size(trialGuesses, 1);
            gains = trialGuesses(j, 1:5);
            neuro = [trialGuesses(j, 6), 0.707];
            display(sprintf('Trying gains: %1.1f, %1.2f, %1.2f, %1.2f, %1.2f and neuro: %1.1f, %1.3f', gains, neuro))
            try
                [tmpPar, tmpFit, ~] = id_rider_siso([disturbanceRunDir filesep ...
                    subset{i}], output, gains, neuro, false, true, ...
                    {M, C1, K0, K2, H}, 'FixedParameter', 'zetanm', 'Focus', 'stability');

                trialParameters(j, :) = tmpPar';
                trialFits(j) = tmpFit;
                display(sprintf('Found gains: %1.1f, %1.2f, %1.2f, %1.2f, %1.2f and neuro: %1.1f, %1.3f for a %1.1f%% fit.', tmpPar, tmpFit))
            catch err
                display(err.identifier)
                trialParameters(j, :) = nan * ones(1, 7);
                trialFits(j) = nan;
            end
        end

        % choose the best fit
        [bestFit, bestIndex] = max(trialFits);
        if isnan(fit) || bestFit > fit
            fit = bestFit;
            par = trialParameters(bestIndex, :);
        end
    end

    if fit < 20 || isnan(fit)
        display('Fit was less than 20%, trying a guess from the Hess method.')
        try
            ronsGuess = guess_from_loop_closure('Rigid', rider, ...
                speed, 'full', par(1:5));
        catch err
            display('Rons guess did not work.')
            display(err.identifier)
            ronsGuess = [];
        end
        if ~isempty(ronsGuess)
            try
                [tmpPar, tmpFit, ~] = id_rider_siso([disturbanceRunDir filesep ...
                    subset{i}], output, ronsGuess, [30, 0.707], false, true, ...
                    {M, C1, K0, K2, H}, 'FixedParameter', 'zetanm', 'Focus', 'stability');
                display(sprintf('Found gains: %1.1f, %1.2f, %1.2f, %1.2f, %1.2f and neuro: %1.1f, %1.3f for a %1.1f%% fit.', tmpPar, tmpFit))
            catch err
                display(err.identifier)
            end
        else
            tmpFit = nan;
        end
        if tmpFit > fit
            fit = tmpFit;
            par = tmpPar;
        end
    end

    display(sprintf('The best fit is %1.0f %% with parameters %1.4f, %1.4f, %1.4f, %1.4f, %1.4f, %1.4f, %1.4f', fit, par))

    speeds(strcmp(matFiles, subset{i})) = speed;
    fits(strcmp(matFiles, subset{i})) = fit;
    parameters(strcmp(matFiles, subset{i}), :) = par;

    % append these results to file
    save('../../data/riderid/bestControllerIdResults.mat', ...
        'parameters', 'fits', 'speeds', 'matFiles', '-append')

    clear fit speed par trialGuesses notTrialGuesses
end
