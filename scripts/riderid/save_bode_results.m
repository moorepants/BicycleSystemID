addpath('../../src/matlab')
config
addpath(PATH_TO_CONTROL_MODEL)

results = load('../../data/riderid/bestControllerIdResults.mat');

% This is the model identified from Luke's pavilion runs.
M = [129.3615, 2.5592;
     2.5592, 0.2505];
C1 = [ 0,   33.5263;
      -0.5486,    2.0997];
K0 = [-115.7074, -4.5261;
      -4.5261,   -0.4889];
K2 = [0, 103.9425;
      0, 2.6034];
H = [0.9017;
     0.0111];

% I want to store the transfer functions (num and den) for the closed loop
% delta and phidot loops, the open loop phi, psi and yq loops [gain
% applied], crossover freqencies and slopes for the three outer
% loops.

numRuns = length(results.matFiles);

% Matlab SUCKS! I had 'closed' as my variable name, but the cell
% initialization was being overwritten in the loop. WTF? Changed it to 'cl'
% and things seem fine. I think I know what it was. It has to do with my
% generate_data function with writes a 'closed' variable to the main
% workspace because linmod can't run in a function. WTF!! This is why Matlab
% needs to figure out namespaces, how fucked up.
cl.delta_deltac = cell(numRuns, 1);
cl.phiDot_phiDotc = cell(numRuns, 1);
cl.yQ_yQc = cell(numRuns, 1);
cl.delta_fB = cell(numRuns, 1);
cl.phiDot_fB = cell(numRuns, 1);
cl.phi_fB = cell(numRuns, 1);
cl.psi_fB = cell(numRuns, 1);
cl.yQ_fB = cell(numRuns, 1);
cl.tDelta_fB = cell(numRuns, 1);

open.phi_ePhi = cell(numRuns, 1);
open.psi_ePsi = cell(numRuns, 1);
open.yQ_eyQ = cell(numRuns, 1);

crossover.phi = zeros(numRuns, 1);
crossover.psi = zeros(numRuns, 1);
crossover.yQ = zeros(numRuns, 1);

w = logspace(-1, 10, 1000);

for i = 1:numRuns
    display(results.matFiles{i})
    if any(isnan(results.parameters(i, :)))
        cl.delta_deltac{i} = nan;
        cl.phiDot_phiDotc{i} = nan;
        cl.yQ_yQc{i} = nan;
        cl.delta_fB{i} = nan;
        cl.phiDot_fB{i} = nan;
        cl.phi_fB{i} = nan;
        cl.psi_fB{i} = nan;
        cl.yQ_fB{i} = nan;
        cl.tDelta_fB{i} = nan;

        open.phi_ePhi{i} = nan;
        open.psi_ePsi{i} = nan;
        open.yQ_eyQ{i} = nan;

        crossover.phi(i) = nan;
        crossover.psi(i) = nan;
        crossover.yQ(i) = nan;
    else
        speed = results.speeds(i);
        parameters = results.parameters(i, :);
        gains = parameters(1:5);
        wnm = parameters(6);

        bicycleName = 'RigidLuke';

        bicycle = bicycle_state_space(bicycleName, speed);
        bicycle = replace_essential(bicycle, M, C1, K0, K2, H, speed, 9.81);

        modelData = generate_data('RigidLuke', speed, ...
            'gains', gains, ...
            'neuroFreq', wnm, ...
            'simulate', false, ...
            'handlingQuality', false, ...
            'forceTransfer', {'Delta', 'PhiDot', 'Phi', 'Psi', 'Y', 'Tdelta'}, ...
            'fullSystem', false, ...
            'display', false, ...
            'plot', false, ...
            'stateSpace', {bicycle.A, bicycle.B, bicycle.C, bicycle.D});

        cl.delta_deltac{i, 1} = ...
            minreal(tf(modelData.closedLoops.Delta.num, ...
                       modelData.closedLoops.Delta.den));
        cl.phiDot_phiDotc{i, 1} = ...
            minreal(tf(modelData.closedLoops.PhiDot.num, ...
                       modelData.closedLoops.PhiDot.den));
        cl.yQ_yQc{i, 1} = minreal(tf(modelData.closedLoops.Y.num, ...
            modelData.closedLoops.Y.den));
        cl.delta_fB{i, 1} = minreal(tf(modelData.forceTF.Delta.num, ...
            modelData.forceTF.Delta.den));
        cl.phiDot_fB{i, 1} = minreal(tf(modelData.forceTF.PhiDot.num, ...
            modelData.forceTF.PhiDot.den));
        cl.phi_fB{i, 1} = minreal(tf(modelData.forceTF.Phi.num, ...
            modelData.forceTF.Phi.den));
        cl.psi_fB{i, 1} = minreal(tf(modelData.forceTF.Psi.num, ...
            modelData.forceTF.Psi.den));
        cl.yQ_fB{i, 1} = minreal(tf(modelData.forceTF.Y.num, ...
            modelData.forceTF.Y.den));
        cl.tDelta_fB{i, 1} = minreal(tf(modelData.forceTF.Tdelta.num, ...
            modelData.forceTF.Tdelta.den));

        open.phi_ePhi{i, 1} = minreal(tf(modelData.openLoops.Phi.num, ...
            modelData.openLoops.Phi.den));
        [mag, ~, freq] = bode(open.phi_ePhi{i}, w);
        crossover.phi(i) = freq(find(20 * log10(mag) < 0, 1));

        open.psi_ePsi{i, 1} = minreal(tf(modelData.openLoops.Psi.num, ...
            modelData.openLoops.Psi.den));
        [mag, ~, freq] = bode(open.psi_ePsi{i}, w);
        crossover.psi(i) = freq(find(20 * log10(mag) < 0, 1));

        open.yQ_eyQ{i, 1} = minreal(tf(modelData.openLoops.Y.num, ...
            modelData.openLoops.Y.den));
        [mag, ~, freq] = bode(open.yQ_eyQ{i}, w);
        crossover.yQ(i) = freq(find(20 * log10(mag) < 0, 1));
    end
end

save('../../data/riderid/transferFunctions.mat', '-struct', 'open')
save('../../data/riderid/transferFunctions.mat', '-append', '-struct', 'cl')
save('../../data/riderid/crossoverFrequencies.mat', '-struct', 'crossover')
