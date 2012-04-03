% this file attempts a double ID process. First, find a black box model for
% the SISO data and secondly use the output from the black box model and the
% experimental input to do a grey box identification.

% load the iddata for run 105
z105 = build_id_data('00105.mat', {'phiDot'});

% take data after the bump
id = z105(8.455 * 200:end); % skips the bump at the beginning
%id = z105(22.44 * 200:25.05 * 200); % just the disturbance

% find the optimal set of gains using the best fit I originally found with
% the transfer function fitting
%gainGuess = [39.3, -0.018, 209.895, 0.081, 0.799];
%omegaGuess = 37.0;

% find the optimal set of gains using Ron's guesses
gainGuess = [37.5, -0.08, 47.2, 0.0951, 0.7475];
omegaGuess = 35;
m = bicycle_grey('Rigid', 3.17, {'phiDot'}, gainGuess, omegaGuess);

time = id.SamplingInstants;

% get the transfer function with Ron's ARX function
[gettf670TF, gd, regmax, err] = gettf1(id.InputData, id.OutputData, [6 7 0], time, 0);
gettf670id = id;
gettf670id.OutputData = lsim(gettf670TF, id.InputData, time);
gettf670s = pem(gettf670id, m);
gettf670output = sim(gettf670s, id.InputData);
gettf670time = figure();
plot(time, id.OutputData, 'b', ...
     time, gettf670id.OutputData, 'r--', ...
     time, gettf670output, 'c', ...
     'LineWidth', 2);
legend('Experimental', 'ARX 670 (gettf1)', 'Grey Box (on ARX 670)')
xlim([20, 40])
saveas(gettf670time, 'plots/double_id_gettf670time.png')

gettf670bode = figure();
h = bodeplot(gettf670TF, 'k--', tf(gettf670s, 'm'), 'k', {1, 20});
set(findall(gcf, 'type', 'line'), 'linewidth', 2)
legend('ARX 670 (gettf1)', 'Grey Box (on ARX 670)')
opts = getoptions(h);
opts.PhaseMatching = 'on';
opts.PhaseMatchingFreq = 1;
opts.PhaseMatchingValue = 180;
opts.Title.String = 'gettf1 6 7 0';
setoptions(h, opts)
saveas(gettf670bode, 'plots/double_id_gettf670bode.png')

[gettf880TF, gd, regmax, err] = gettf1(id.InputData, id.OutputData, [8 8 0], time, 0);
gettf880id = id;
gettf880id.OutputData = lsim(gettf880TF, id.InputData, time);
gettf880s = pem(gettf880id, m);
gettf880output = sim(gettf880s, id.InputData);
gettf880time = figure();
plot(time, id.OutputData, 'b', ...
     time, gettf880id.OutputData, 'r--', ...
     time, gettf880output, 'c', ...
     'LineWidth', 2);
legend('Experimental', 'ARX 880 (gettf1)', 'Grey Box (on ARX 880)')
xlim([20, 40])
saveas(gettf880time, 'plots/double_id_gettf880time.png')

gettf880bode = figure();
h = bodeplot(gettf880TF, 'k--', tf(gettf880s, 'm'), 'k', {1, 20});
set(findall(gcf, 'type', 'line'), 'linewidth', 2)
legend('ARX 880 (gettf1)', 'Grey Box (on ARX 880)')
opts = getoptions(h);
opts.PhaseMatching = 'on';
opts.PhaseMatchingFreq = 1;
opts.PhaseMatchingValue = 180;
opts.Title.String = 'gettf1 8 8 0';
setoptions(h, opts)
saveas(gettf880bode, 'plots/double_id_gettf880bode.png')

%%%% use matlab's toolbox for some models at the "right" order
%%%% output error model
%%%oe881 = oe(id, [8 8 1]);
%%%output = sim(oe881, id);
%%%oe881id = id;
%%%oe881id.OutputData = output.OutputData;
%%%oe881s = pem(oe881id, m);
%%%compare(id, oe881, oe881s)
%%%h = bodeplot(oe881, oe881s);
%%%opts = getoptions(h);
%%%opts.PhaseMatching = 'on';
%%%opts.PhaseMatchingFreq = 1;
%%%opts.PhaseMatchingValue = 180;
%%%setoptions(h, opts)
%%%
%%%% arx model
%%%arx881 = arx(id, [8 8 1]);
%%%output = sim(arx881, id);
%%%arx881id = id;
%%%arx881id.OutputData = output.OutputData;
%%%arx881s = pem(arx881id, m);
%%%compare(id, arx881, arx881s)
%%%h = bodeplot(arx881, arx881s);
%%%opts = getoptions(h);
%%%opts.PhaseMatching = 'on';
%%%opts.PhaseMatchingFreq = 1;
%%%opts.PhaseMatchingValue = 180;
%%%setoptions(h, opts)
