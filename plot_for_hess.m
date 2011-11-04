% load the iddata for run 105
z105 = build_id_data('00105.mat', {'phiDot'});

% take data after the bump
id = z105(8.455 * 200:end); % skips the bump at the beginning
%id = z105(22.44 * 200:25.05 * 200); % just the disturbance

time = id.SamplingInstants;

% get the transfer function with Ron's ARX function
[g, gd, regmax, err] = gettf1(id.InputData, id.OutputData, [6 7 0], time, 0);

config
addpath(PATH_TO_CONTROL_MODEL)
% use Ron's manually found gains with the control model
%kagress = 4.0;
%manualGains = [37.5, -0.08, 9.44 * 1.2 * kagress, 0.0951, 0.1495 * kagress]
% new set of parameters provided by Ron on Nov 1 2011
manualGains = [37.5, -0.08, 47.2, 0.0951, 0.7475];
omega = 35;
hess = generate_data('Rigid', 3.17, ...
                     'gains', manualGains, ...
                     'neuroFreq', omega, ...
                     'simulate', 0, ...
                     'loopTransfer', 0, ...
                     'handlingQuality', 0);

% find the optimal set of gains using the best fit I originally found with
% the transfer function fitting
gainGuess = [39.3, -0.018, 209.895, 0.081, 0.799];
omegaGuess = 37.0;

% use Ron's guess as a starting point
gainGuess = manualGains;
omegaGuess = omega;

m = bicycle_grey('Rigid', 3.17, {'phiDot'}, gainGuess, omegaGuess);
s = pem(id, m);

% now lets id on the output associated with the transfer function found from
% gettf1
id2 = id;
id2.OutputData = lsim(g, id.InputData, time);
s2 = pem(id2, m);

% plot the Bode diagrams
hold all
bodeDiagrams = figure(1);
blackbox = bodeplot(g);
hessTF = tf(hess.forceTF.PhiDot.num, hess.forceTF.PhiDot.den)
manual = bodeplot(hessTF);
grey2 = bodeplot(ss(s2.A, s2.B, s2.C, s2.D));
greybox = bodeplot(ss(s.A, s.B, s.C, s.D));
hold off
opts = getoptions(greybox)
opts.PhaseMatching = 'on';
opts.PhaseMatchingValue = 0;
opts.XLim = {[1, 20]};
opts.YLim(1) = {[-80, -30]};
opts.YLim(2) = {[-720, 0]};
setoptions(greybox, opts)
legend('Black Box', 'Model - Manual', 'Grey from clean', 'Grey Box - Optimal')
saveas(bodeDiagrams, 'plots/bode-for-delft-talk.png')

% plot the time histories for each model versus the experimental
timeHistories = figure(2);
yBlack = lsim(g, id.InputData, time);
yHess = lsim(hessTF, id.InputData, time);
yGrey = lsim(ss(s.A, s.B, s.C, s.D), id.InputData, time);
yGrey2 = lsim(ss(s2.A, s2.B, s2.C, s2.D), id.InputData, time);
plot(time, id.OutputData, 'k', time, yBlack, time, yHess, time, yGrey2, time, yGrey);
xlim([20, 40])
legend('Experiment', 'Black Box', 'Hess', 'Grey2', 'Grey Box');
saveas(timeHistories, 'plots/time-history-for-delft-talk.png')
