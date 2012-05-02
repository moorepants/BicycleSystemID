% Jason K. Moore, January 22, 2012
%
% This script identifies the parameters of the control system for several
% line tracking runs at about 5 m/s for the three riders. It plots the
% variation of the parameters in a box plot and creates bode plots for each
% of the lateral force transfer functions with the mean frequency response
% for each rider. We looked at these to see if the frequency response varied
% much even though the parameters varied.

addpath('..')

% these are runs from Jason, Charlie and Luke at the same speed
jasonRuns = 527:535;
charlieRuns = 602:611;
lukeRuns = [693:700 704];

runs = [jasonRuns, charlieRuns, lukeRuns];

%% load the data from the runs
speeds = zeros(length(runs), 1);
riders = cell(length(runs), 1);
data = cell(length(runs), 1);

outputs = {'phiDot', 'delta', 'psiDot', 'phi', 'deltaDot'};
inputs = {'fB'};

for i = 1:length(runs)
    [z, speed, rider] = build_id_data(['00' num2str(runs(i)) '.mat'], ...
        outputs, inputs, '../pavilion/mat');
    speeds(i) = speed;
    riders{i} = rider;
    data{i} = z;
end

uniqueRiders = unique(riders);

display(speeds)

%% find the parameters that give the best fit

models = cell(length(runs), 1);
parameters = zeros(length(runs), 6);

for i = 1:length(runs)
    display(sprintf('Identifying run %s', num2str(runs(i))))
    guess = [9.7676, -1.5008, 3.0039, 0.7115, 0.1670, 43.3128];
    m = bicycle_grey(['Rigid' rider], speeds(i), outputs, guess(1:5), ...
        guess(6));
    s = pem(data{i}, m);
    parameters(i, :) = s.par';
    models{i} = s;
end

display(parameters)

%% plot the output comparisons for each model
for i = 1:length(runs)
    figure()
    compare(data{i}, models{i})
end

%% create box plots comparing the parameter variation for each rider
parameterNames = {'kDelta', 'kPhiDot', 'kPhi', 'kPsi', 'kYq', 'wnm'};

parVar = figure();
for i = 1:size(parameters, 2)
    subplot(3, 2, i)
    perRider = nan * ones(size(parameters, 1), length(uniqueRiders));
    labels = cell(length(uniqueRiders), 1);
    for j = 1:length(uniqueRiders)
        riderInd = find(strcmp(uniqueRiders{j}, riders));
        perRider(1:length(riderInd), j) = parameters(riderInd, i);
        labels{j} = [uniqueRiders{j} ' n=' num2str(length(riderInd))];
    end
    boxplot(perRider, labels)
    ylabel(parameterNames{i})
end
saveas(parVar, '../plots/rider-par-var-single-speed.png')

%% compute the input to outputs frequency responses for each model

n = 200;
w = logspace(-1, 1.3, n);

magnitudes = zeros(length(runs), length(outputs), n);
phases = zeros(length(runs), length(outputs), n);

for i = 1:length(runs)
    [mag, phase] = bode(models{i}, w);
    magnitudes(i, :, :) = squeeze(mag);
    phases(i, :, :) = squeeze(phase);
end

%% create Bode plots for each transfer function

colors = {'b', 'r', 'g'}; % only had three riders

for i = 1:length(uniqueRiders)

    riderMags = magnitudes(find(strcmp(riders, uniqueRiders{i})), :, :);
    magMean = squeeze(mean(20 * log10(riderMags)));
    magStd = squeeze(std(20 * log10(riderMags)));

    riderPhases = phases(find(strcmp(riders, uniqueRiders{i})), :, :);
    phaseMean = squeeze(mean(riderPhases));
    phaseStd = squeeze(std(riderPhases));

    for j = 1:length(outputs)

        if i > 1
            figure(plots.(outputs{j}))
        else
            plots.(outputs{j}) = figure();
        end

        subplot(211)
        if i > 1
            hold on
        end
        magLines = semilogx(w, magMean(j, :), colors{i}, ...
            w, magMean(j, :) + magStd(j, :), [colors{i} '--'], ...
            w, magMean(j, :) - magStd(j, :), [colors{i} '--']);
        xlim([w(1), w(end)])
        ylabel('Magnitude [dB]')
        title(sprintf('%s / LateralForce', outputs{j}))
        if i > 1
            hold off
        end

        subplot(212)
        if i > 1
            hold on
        end
        phaseLines = semilogx(w, phaseMean(j, :), colors{i}, ...
            w, phaseMean(j, :) + phaseStd(j, :), [colors{i} '--'], ...
            w, phaseMean(j, :) - phaseStd(j, :), [colors{i} '--']);
        xlim([w(1), w(end)])
        ylabel('Phase [deg]')
        xlabel('Frequency [rad/s]')
        if i > 1
            hold off
        end

        set(magLines(1), 'Linewidth', 2)
        set(phaseLines(1), 'Linewidth', 2)

    end
end

for i = 1:length(outputs)
    saveas(plots.(outputs{i}), ['../plots/freq-compare-single-speed' outputs{i} '.png'])
end

    %transferFunctions = tf(s, 'm');
%
    %for j = 1:length(outputs)
%
        %bodePlots.(outputs{j}) = figure(length(runs) + j);
        %hold all
        %if strcmp(rider, 'Charlie')
            %bode(transferFunctions(j), 'r', {1, 20})
        %elseif strcmp(rider, 'Luke')
            %bode(transferFunctions(j), 'b', {1, 20})
        %end
        %hold off
    %if i == length(runs)
        %legend(strread(num2str(runs),'%s'))
    %end
%
    %end

%end
%saveas(fig, ['plots/frequency-compare.png'])
