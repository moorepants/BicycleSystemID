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

        print(f1, '-dpng', ['../../plots/riderid/bode/mean-' tfNames{j} '-' num2str(v) '.png'], '-r150')

        f2 = figure();
        set(f2,'PaperUnits','inches','PaperPosition',[0 0 4 3])
        %set(f2, 'Visible', 'off') % if this is on the plot doesn't seem to
        %show up
        bopt.Title.String = tit;
        bode(f2, trfc{:}, w, bopt)
        print(f2, '-dpng', ['../../plots/riderid/bode/all-' tfNames{j} '-' num2str(v) '.png'], '-r150')
        close all
    end
end
