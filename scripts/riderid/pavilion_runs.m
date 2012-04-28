contents = what('../../data/riderid/pavilion');
matFiles = contents.mat;
outputs = {'phiDot'};

load('../../guesses.mat')

% speed, kDelta, kPhiDot, kPhi, kPsi, kY, wnm
%3.85, 11.1857, -1.6166, 3.7822, 0.4393, 0.3056, 45.0976
%4.959, 11.1048, -1.0180, 3.8664, 0.4166, 0.1014, 42.9452
%5.667, 4.0937, -3.2782, 1.9031, 0.7296, 0.1174, 42.2266
%5.585, 3.0943, -4.5908, 1.5155, 0.9361, 0.0462, 40.9320

speedGuesses = [speedGuesses; 3.85; 4.959; 5.667; 5.585; 1.0];

moreGuesses = [11.1857 -1.6166 3.7822 0.4393 0.3056 45.0976
               11.1048 -1.0180 3.8664 0.4166 0.1014 42.9452
               4.0937 -3.2782 1.9031 0.7296 0.1174 42.2266
               3.0943 -4.5908 1.5155 0.9361 0.0462 40.9320
               37.5 -0.08 47.2 0.0951 0.7475 35];

parameterGuesses = [parameterGuesses; moreGuesses];

%21.7853   -0.3561    6.6855    0.3814    0.0006   40.9405
%7.0581   -1.5155    2.4920    0.5435    0.0937   41.7825
%2.3386   -6.1365    1.3539    1.0447    0.0407   40.6199
%1.7623   -7.9026    1.2503    1.1118    0.0431   40.4697
%1.8206   -8.1939    1.7427    0.8872    0.0268   49.6578
%1.3674  -10.2213    1.6922    0.9339    0.0610   50.4560
%9.2945   -1.9182    3.4331    0.5301    0.2304   41.6667
%2.8433   -4.5740    2.8468    0.5117    0.3707   76.3847
%2.5817   -7.1793    1.5493    1.3446    0.0000   55.0865
%2.0122   -8.0565    2.7209    0.5536    0.2326   67.1306];

parameters = zeros(length(matFiles), 6);
fits = zeros(length(matFiles), 1);
speeds = zeros(length(matFiles), 1);

for i = 1:length(matFiles)
    % build the identification data
    [z, speed, rider] = build_id_data(matFiles{i}, outputs, 'pavilion/mat');

    % select some guesses from similar speeds
    trialIndices = find((speedGuesses < (speed + 0.5)) & (speedGuesses > (speed - 0.5)));
    notTrialIndices = find(~((speedGuesses < (speed + 0.5)) & (speedGuesses > (speed - 0.5))));

    if isempty(trialIndices)
        trialGuesses = [37.5 -0.08 47.2 0.0951 0.7475 35];
        notTrialGuesses = parameterGuesses;
        r = randperm(size(parameterGuesses, 1));
    else
        trialGuesses = parameterGuesses(trialIndices, :);
        % add a couple of random guesses from other speeds
        notTrialGuesses = parameterGuesses(notTrialIndices, :);
        r = randperm(size(notTrialGuesses, 1));
    end
    trialGuesses = [trialGuesses; notTrialGuesses(r(1:2), :)];
    % add the guess generated from Ron's loop closure method
    try
        ronsGuess = guess_from_loop_closure('Rigid', rider, speed, 'full');
        trialGuesses = [trialGuesses; ronsGuess];
    catch
        display('Rons guess did not work.')
    end

    % initialize
    trialFits = zeros(size(trialGuesses, 1), 1);
    trialParameters = zeros(size(trialGuesses, 1), 6);

    for j = 1:size(trialGuesses, 1);
        guess = trialGuesses(j, :);
        try
            display('Building grey model')
            m = bicycle_grey(['Rigid' rider], speed, outputs, guess(1:5), guess(6));
            display('Solving model')
            s = pem(z, m);
            [~, fit, ~] = compare(z, m);
            display(sprintf('%s, %s, %1.3f m/s, %1.0f', matFiles{i}, rider, speed, fit))
            display(s.ParameterVector')
            trialParameters(j, :) = s.ParameterVector';
            trialFits(j) = fit;
        catch
            display('Error')
            trialParameters(j, :) = [nan, nan, nan, nan, nan, nan];
            trialFits(j) = nan;
        end
    end
    % choose the best fit
    [bestFit, bestIndex] = max(trialFits);
    parameters(i, :) = trialParameters(bestIndex, :);
    fits(i) = bestFit;
    speeds(i) = speed;
    display(sprintf('The best fit is %1.0f %% with parameters %1.4f, %1.4f, %1.4f, %1.4f, %1.4f, %1.4f', fits(i), parameters(i, :)))
end
