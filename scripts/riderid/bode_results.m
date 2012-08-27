% This loads the many tranfer functions for each run and plots a ton of Bode
% plots.

results = load('../../data/riderid/bestControllerIdResults.mat');
transferFuncs = load('../../data/riderid/transferFunctions.mat');

tfNames = fieldnames(transferFuncs);
w = logspace(-1, 2, 500);
%speedBins = [2.0, 3.0, 4.0, 4.92, 5.8, 7.0, 9.0];
speedBins = [2.25, 3.375, 4.0, 4.5, 5.0, 5.625, 7.375];

bopt = bodeoptions;
bopt.PhaseMatching = 'on';
bopt.PhaseMatchingFreq = 0.1;
bopt.PhaseMatchingValue = 0;
bopt.Title.Interpreter = 'None';

speedBinRange = 0.5;

for k = 1:length(speedBins)
    v = speedBins(k);
    indices = find(v - speedBinRange / 2 < results.speeds & ...
        results.speeds < v + speedBinRange / 2);

    speeds = results.speeds(indices);
    matFiles = results.speeds(indices);
    numRuns = length(matFiles);

    for j = 1:length(tfNames)
        trfc = transferFuncs.(tfNames{j});
        trfc = trfc(indices);

        magnitudes = zeros(numRuns, length(w));
        phases = zeros(numRuns, length(w));

        for i = 1:numRuns
            try
                [mag, phase] = bode(trfc{i}, w, bopt);
                magnitudes(i, :) = squeeze(mag);
                phases(i, :) = squeeze(phase);
            catch
                magnitudes(i, :) = ones(length(w), 1) * nan;
                phases(i, :) =  ones(length(w), 1) * nan;
            end
        end

        meanMag = nanmean(20 * log10(magnitudes));
        stdMag = nanstd(20 * log10(magnitudes));

        f1 = figure();
        set(f1, 'Visible', 'off')
        set(f1,'PaperUnits','inches','PaperPosition',[0 0 4 3])
        subplot(2, 1, 1)
        semilogx(w, meanMag, w, meanMag + stdMag, 'b:', w, meanMag - stdMag, 'b--')
        tit = ['Loop: ' tfNames{j} ', Speed: ' num2str(v) ' +/- ' ...
            num2str(speedBinRange / 2) ' m/s, n = ' num2str(length(indices))];
        title(tit, 'Interpreter', 'None')

        meanPhase = nanmean(phases, 1);
        stdPhase = nanstd(phases, 1);

        subplot(2, 1, 2)
        semilogx(w, meanPhase, w, meanPhase + stdPhase, 'b:', w, meanPhase - stdPhase, 'b--')

        print(f1, '-dpng', ['../../plots/riderid/bode/mean-' tfNames{j} '-' strrep(num2str(v), '.', 'p') '.png'], '-r150')

        % when I have the visibility off the plot has no lines, I can't
        % figure that out
        %f2 = figure('Visible', 'off');
        f2 = figure();
        set(f2,'PaperUnits','inches','PaperPosition',[0 0 4 3])
        bopt.Title.String = tit;
        hold on
        bode(f2, trfc{:}, w, bopt)
        lines = findobj(f2, 'type', 'line');
        set(lines, 'color', [0.7, 0.7, 0.7])
        ax = findobj(f2, 'type', 'axes');
        plot(ax(2), w, meanMag, 'k-', 'LineWidth', 2)
        plot(ax(2), w, meanMag + stdMag, 'k--', 'LineWidth', 1.5)
        plot(ax(2), w, meanMag - stdMag, 'k--', 'LineWidth', 1.5)
        phaseLines = findobj(ax(1), 'type', 'line');
        phaseYData = get(phaseLines, 'YData');
        phases2 = zeros(size(phases, 1), length(w));
        m = 1;
        for n = 1:length(phaseYData)
            if length(phaseYData{n}) == 500
                phases2(m, :) = phaseYData{n};
                m = m + 1;
            end
        end
        meanPhases2 = mean(phases2, 1);
        stdPhases2 = std(phases2, 1);
        plot(ax(1), w, meanPhases2, 'k-', 'LineWidth', 2)
        plot(ax(1), w, meanPhases2 + stdPhases2, 'k--', 'LineWidth', 1.5)
        plot(ax(1), w, meanPhases2 - stdPhases2, 'k--', 'LineWidth', 1.5)
        hold off
        print(f2, '-dpng', ['../../plots/riderid/bode/all-' tfNames{j} '-' strrep(num2str(v), '.', 'p') '.png'], '-r150')
        close all
    end
end
