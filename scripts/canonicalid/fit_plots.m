PATH_TO_BICYCLE_SYSTEM_ID = '/media/Data/Documents/School/UC Davis/Bicycle Mechanics/BicycleSystemID';
addpath(PATH_TO_BICYCLE_SYSTEM_ID)
dataDir = [PATH_TO_BICYCLE_SYSTEM_ID filesep 'scripts' filesep 'exports'];

inputs = {'tDelta', 'fB'};
states = {'phi', 'delta', 'phiDot', 'deltaDot'};
outputs = {'phi', 'delta', 'phiDot', 'deltaDot'};

runID = '00311';
[data, v, rider] = build_id_data([runID '.mat'], outputs, inputs, ...
    dataDir, true);

v

whippleModel = bicycle_structured(['Rigid' rider], v, 'states', states, ...
    'inputs', inputs, 'outputs', outputs);

% the arm model
arm = load(['armsAB-' rider '.mat']);

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

% run 588
%identifiedModel.A = [0, 0, 1, 0;
                     %0, 0, 0, 1;
                     %6.60296243, -3.31536843, -1.96360195e-01, -8.02745331e-02;
                     %1.16759687e+02, 2.33834646e+01, 1.11657551e+01, -1.64972573e+01];
% run 670
%identifiedModel.A = [  0.        ,    0.        ,    1.        ,    0.;
          %0.        ,    0.        ,    0.        ,    1.;
          %6.60296243,  -16.34541958,   -0.45606932,   -0.18644691;
        %116.75968666,  -35.35707903,   25.93376111,  -38.31679315];

% run 253
%identifiedModel.A = [   0.        ,    0.        ,    1.        ,    0.;
                       %0.        ,    0.        ,    0.        ,    1.;
          %6.60296243,  -13.61579371,   -0.41533524,   -0.1697943;
        %116.75968666,  -23.0517008 ,   23.61747267,  -34.89450725];
%data = data(25*200:200*45);

identifiedModel.A = [   0.        ,    0.        ,    1.        ,    0.;
          0.        ,    0.        ,    0.        ,    1.;
          6.60296243,  -39.77382541,   -0.71600323,   -0.29271118;
        116.75968666, -140.97426193,   40.71454892,  -60.155214  ];
data = data(32*200:200*62);

identifiedModel.B = [ 0, 0;
                      0, 0;
                     -0.11000105, 0.00787264;
                      6.25505994, -0.04681092];
[YH, FIT, X0] = compare(data, identifiedModel, whippleModel, armModel);

time = data.SamplingInstants;

fig = figure('Visible', 'off');
figWidth = 6.;
figHeight = figWidth * 1.25;
set(gcf, ...
    'Color', [1, 1, 1], ...
    'PaperOrientation', 'portrait', ...
    'PaperUnits', 'inches', ...
    'PaperPositionMode', 'manual', ...
    'OuterPosition', [424, 305 - 50, 518, 465], ...
    'PaperPosition', [0, 0, figWidth, figHeight], ...
    'PaperSize', [figWidth, figHeight])

axesHandles = tight_subplot(6, 1, [0.05, 0.0], [0.1, 0.1], [0.17, 0.05]);

ax = axesHandles(1);
lh = plot(ax, time, data.InputData(:, 1), 'g');
title(ax, [runID ' ' rider ' ' num2str(v) ' m/s'])
ylabel(ax, '\(T_\delta\) [N-m]', 'Interpreter', 'Latex')

ax = axesHandles(2);
lh = plot(ax, time, data.InputData(:, 2), 'g');
ylabel(ax, '\(F_{c_l}\) [N]', 'Interpreter', 'Latex')

ylabels = {'\(\phi\) [rad]', '\(\delta\) [rad]',
           '$\dot{\phi}$ [rad/s]', '$\dot{\delta}$ [rad/s]'};

f = squeeze(FIT);
for i = 1:4
    ax = axesHandles(i + 2);
    lh = plot(ax, ...
              time, data.OutputData(:, i), 'g', ...
              time, YH{1}.OutputData(:, i), 'k', ...
              time, YH{2}.OutputData(:, i), 'b', ...
              time, YH{3}.OutputData(:, i), 'r');
    ylabel(ax, ylabels{i}, 'Interpreter', 'Latex')
    idLeg = sprintf('I %1.0f%%', f(1, i));
    whipLeg = sprintf('W %1.0f%%', f(2, i));
    armLeg = sprintf('A %1.0f%%', f(3, i));
    leg = legend(ax, 'M', idLeg, whipLeg, armLeg);
    set(leg, 'FontSize', 5)
end

% add the the time axis to the last bottom plot
xlabel('Time [s]')

saveas(fig, ['canonical-id-plots/example-fit-' runID '.pdf'])
%print(fig, '-dpng', 'canonical-id-plots/example-fit.png', '-r200')
