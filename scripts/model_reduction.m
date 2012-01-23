outputs = {'delta', 'phiDot'};

[z, speed] = build_id_data('00531.mat', outputs);
id = detrend(z(1.75 * 200:6.11 * 200), 1);

%speed = 3.85
%z = build_id_data('00516.mat', outputs);
%id = detrend(z(3.21 * 200:10.55 * 200), 1);

% black box model
pem8 = pem(id, 8);

% reduced model
for i = 1:length(outputs)
    zSISO = build_id_data('00516.mat', outputs(i));
    idSISO = detrend(zSISO(3.58 * 200:10.55 * 200), 1);
    eval(['m' num2str(i) ' = idss(oe(idSISO, [8 8 1]));'])
end
reducedModel = balred([m1; m2], 8);
reduced = pem(id, reducedModel);

% grey box model
% Ron's guess
%gainGuess = [37.5, -0.08, 47.2, 0.0951, 0.7475];
%omegaGuess = 35;
% previous solution 516
gainGuess = [11.1857, -1.6166, 3.7822, 0.4393, 0.3056];
omegaGuess = 45.0976;

m = bicycle_grey('Rigid', speed, outputs, gainGuess, omegaGuess);
grey = pem(id, m);

figure()
compare(id, pem8, reduced, grey)
figure()
bode(pem8, reduced, grey, {1, 20})
