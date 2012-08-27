% This script tests fitting different numbers of outputs for a single run.
addpath('../../src/matlab')

allOutputs = {'phiDot', 'delta', 'psiDot', 'phi', 'deltaDot', 'tDelta', 'yQ'};
inputs = {'fB'};
parameters = zeros(6, 7);
runNum = '605';
for i = 6 %= 1:length(allOutputs)
    outputs = allOutputs(1:i);
    [z, meta] = build_id_data(['00' runNum '.mat'], outputs, inputs, '../../data/riderid/pavilion/mat', true);
    gains = [9.7676, -1.5008, 3.0039, 0.7115, 0.1670];
    neuro = [43.3128, 0.707];
    m = bicycle_grey('lateral', ['Rigid' meta.rider], meta.speed, outputs, gains, neuro);
    s = pem(z, m, 'FixedParameter', {'zetanm'});
    parameters(i, :) = s.par';
    compare(z, s)
    saveas(gcf, ['../../plots/riderid/run' runNum '-' num2str(i) '.png'])
    %close all
end
display(parameters)
