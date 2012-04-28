outputs = {'phiDot'};

d = '../../scripts/statespaceid/exports';
[data, speed, rider] = build_id_data('00700.mat', outputs, {'fB'}, d, 1);

bicycleName = ['Rigid' rider];

bStates = {'yP', 'psi', 'phi', 'delta', 'phiDot', 'deltaDot'};
bOutputs = {'psi', 'phi', 'delta', 'phiDot', 'deltaDot', 'yQ'};
bInputs = {'tDelta', 'fB'};
bicycle = bicycle_state_space(bicycleName, speed, 'states', bStates, ...
    'inputs', bInputs, 'outputs', bOutputs);

load('../../../CanonicalBicycleID/data/cid-L-P.mat')
invM = inv(M);
bicycle.A(5:6, 3:6) = [-invM * [K0 * 9.81 + K2 * speed^2], -invM * C1 * speed];
bicycle.B(5:6, :) = [invM(:, 2), invM * H];

bicycle = bicycle_state_space(bicycleName, speed);
A = [-invM * [K0 * 9.81 + K2 * speed^2], -invM * C1 * speed];
B = [invM(:, 2), invM * H];
bicycle.A([9, 11], [4, 7, 9, 11]) = A;
bicycle.B([9, 11], [2, 3]) = B;
gainGuess = [121.0000   -0.0150   28.4447    0.0869    0.0805];
hess = generate_data(bicycleName, speed, ...
                     'stateSpace', {bicycle.A, bicycle.B, bicycle.C, bicycle.D}, ...
                     'gainGuess', gainGuess, ...
                     'crossover', [1.5, 0.75, 0.375], ...
                     'loopTransfer', 0, ...
                     'simulate', 0, ...
                     'handlingQuality', 0, ...
                     'forceTransfer', {}, ...
                     'fullSystem', 0);
%gains = [9.7676, -1.5008, 3.0039, 0.7115, 0.1670];
%neuro = [43.3128, 0.707];
gains = [hess.modelPar.kDelta, hess.modelPar.kPhiDot, hess.modelPar.kPhi, hess.modelPar.kPsi, hess.modelPar.kY];
neuro = [30.0, 0.707];
grey = bicycle_grey('lateral', bicycle, speed, outputs, gains, neuro);

s = pem(data, grey, 'FixedParameter', {'wnm', 'zetanm'});
s2 = pem(data, s, 'FixedParameter', {'zetanm'});

compare(data, grey, s, s2)
