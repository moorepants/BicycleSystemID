outputs = {'delta', 'phiDot'};
z = build_id_data('00516.mat', outputs);
id = detrend(z(3.58 * 200:10.55 * 200), 1);

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
gainGuess = [37.5, -0.08, 47.2, 0.0951, 0.7475];
omegaGuess = 35;
m = bicycle_grey('Rigid', 3.85, outputs, gainGuess, omegaGuess);
grey = pem(id, m);

compare(id, pem8, reduced, grey)
bode(pem8, reduced, grey, {1, 20})
