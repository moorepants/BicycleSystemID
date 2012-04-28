PATH_TO_BICYCLE_SYSTEM_ID = '/media/Data/Documents/School/UC Davis/Bicycle Mechanics/BicycleSystemID';
addpath([PATH_TO_BICYCLE_SYSTEM_ID '/src/matlab'])
dataDir = [PATH_TO_BICYCLE_SYSTEM_ID filesep 'scripts' filesep ...
    'statespaceid' filesep 'exports'];

inputs = {'tDelta', 'fB'};
states = {'phi', 'delta', 'phiDot', 'deltaDot'};
%outputs = {'phi', 'delta', 'phiDot', 'deltaDot'};
%outputs = {'phiDot', 'deltaDot'};
outputs = {'phiDot'};

runID = '00699';
[data, v, rider] = build_id_data([runID '.mat'], outputs, inputs, ...
    dataDir, true);

whippleModel = bicycle_structured(['Rigid' rider], v, 'states', states, ...
    'inputs', inputs, 'outputs', outputs);

whippleModel.As(3, 2) = whippleModel.A(3, 2);

% the arm model
arm = load(['../canonicalid/armsAB-' rider '.mat']);

C = zeros(4, 19);
C(1, 4) = 1;
C(2, 7) = 1;
C(3, 17) = 1;
C(4, 19) = 1;
indice = round(v * 10) + 1
armSS = ss(squeeze(arm.stateMatrices(indice, :, :)), ...
    squeeze(arm.inputMatrices(indice, :, [3, 4])), ...
    C, 0, 'InputName', {'tDelta', 'fB'}, ...
    'OutputName', {'phi', 'delta', 'phiDot', 'deltaDot'});
armModel = idss(armSS);

% this model was computed with the benchmark canonical identification for
% all riders and both environments @ 2.0092 m/s
identifiedModel = whippleModel;

M = [129.86277082, 2.28375768;
       2.28375768, 0.20003257];
C1 = [  0.        , 23.94008131;
       -0.88845094,  1.73368327];
K0 = [-114.59019659, -3.91797216;
        -3.91797216, -0.66780731];
K2 = [0., 102.9447477;
      0.,   2.33973324];
H = [ 0.91545833;
      0.0086155 ];

MI = inv(M);
Ar = -MI * (9.81 * K0 + v^2 * K2);
Al = -MI * v * C1;

BT = MI;
BF = BT * H;

identifiedModel.A = [zeros(2, 2), eye(2);
                    Ar, Al];
identifiedModel.B = [zeros(2, 2);
                    BT(:, 2), BF];
pemArgs = { ...
           'Focus', 'Prediction', ...
           'DisturbanceModel', 'Estimate', ...
           'InitialState', 'Estimate', ...
           'SearchMethod', 'lm', ... 'lsqnonlin', ...
           'Maxiter', 1000, ...
           'LimitError', 0.0, ...
           'Tolerance', 0.005, ...
           'Display', 'on', ...
          };

           %'Criterion', 'trace', ...
           %'Weighting', [1, 0, 0, 0;
                         %0, 1, 0, 0;
                         %0, 0, 1, 0;
                         %0, 0, 0, 1], ...
new = pem(data, identifiedModel, pemArgs{:})

compare(data, identifiedModel, whippleModel, armModel, new);
