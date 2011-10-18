config
addpath(PATH_TO_CONTROL_MODEL)

% build the grey box model
aux.bicycle = 'Rigid';
aux.speed = 7.0;
% this is the answer we were getting with the roll rate arx
%guess = [39.3, -0.018, 209.895, 0.081, 0.799, 37];
% thes are the gains from Ron's method
guess = [76.3808, -0.0516, 7.2456, 0.2632, 0.0708, 30];
data = generate_data(aux.bicycle, aux.speed, ...
                     'simulate', 0, ...
                     'loopTransfer', 0, ...
                     'handlingQuality', 0, ...
                     'gains', guess(1:5), ...
                     'neuroFreq', guess(6), ...
                     'forceTransfer', {}, ...
                     'display', 0);
% get the state space model so we don't have to calculate it during the gain
% matching algorithm
aux.stateSpace = {data.modelPar.A,
                  data.modelPar.B,
                  data.modelPar.C,
                  data.modelPar.D};

m = idgrey('linear_parameterization', guess, 'c', aux, ...
    'DisturbanceModel', 'Estimate');

% load the data
z = build_id_data('00264.mat', '00265.mat');
