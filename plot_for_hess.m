% load the iddata for run 105
z105 = build_id_data('00105.mat', {'phiDot'});

% take data after the bump
id = z105(8.455 * 200:end);

time = id.SamplingInstants;

% get the transfer function with Ron's ARX function
[g, gd, regmax, err] = gettf1(id.InputData, id.OutputData, [6 7 0], time, 0);

% use Ron's manually found gains with the control model
config
addpath(PATH_TO_CONTROL_MODEL)
kagress = 4.0;
manualGains = [37.5, -0.08, 9.44 * 1.2 * kagress, 0.0951, 0.1495 * kagress]
omega = 35;
hess = generate_data('Rigid', 3.0, ...
                     'gains', manualGains, ...
                     'neuroFreq', omega, ...
                     'simulate', 0, ...
                     'loopTransfer', 0, ...
                     'handlingQuality', 0);

% find the optimal set of gains using the best fit I originally found with
% the transfer function fitting
gainGuess = [39.3, -0.018, 209.895, 0.081, 0.799];
omegaGuess = 37.0;
m = bicycle_grey('Rigid', 3.0, {'phiDot'}, gainGuess, omegaGuess);
s = pem(id, m);

% plot the Bode diagrams
hold all
bodeDiagrams = figure(1);
blackbox = bodeplot(g);
hessTF = tf(hess.forceTF.PhiDot.num, hess.forceTF.PhiDot.den)
%manual = bodeplot(hessTF);
greybox = bodeplot(ss(s.A, s.B, s.C, s.D));
hold off
opts = getoptions(greybox)
opts.PhaseMatching = 'on';
opts.PhaseMatchingValue = 0;
opts.XLim = {[1, 20]};
opts.YLim(1) = {[-80, -30]};
opts.YLim(2) = {[-720, 0]};
setoptions(greybox, opts)
%legend('Black Box', 'Model - Manual', 'Grey Box - Optimal')
legend('Black Box', 'Grey Box')
saveas(bodeDiagrams, 'plots/bode-for-delft-talk.png')

% plot the time histories for each model versus the experimental
timeHistories = figure(2);
yBlack = lsim(g, id.InputData, time);
yHess = lsim(hessTF, id.InputData, time);
yGrey = lsim(ss(s.A, s.B, s.C, s.D), id.InputData, time);
%plot(time, id.OutputData, 'k', time, yBlack, time, yHess, time, yGrey);
plot(time, id.OutputData, 'k', time, yBlack, time, yGrey);
xlim([20, 40])
legend('Experiment', 'Black Box', 'Grey Box');
saveas(timeHistories, 'plots/time-history-for-delft-talk.png')
