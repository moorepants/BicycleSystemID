% This script tests fitting different numbers of outputs for a single run.

allOutputs = {'phiDot'}; %, 'delta', 'psiDot', 'phi', 'deltaDot', 'tDelta', 'yP'};
inputs = {'fB'};
parameters = zeros(6, 6);
for i = 1:length(allOutputs)
    outputs = allOutputs(1:i);
    [z, speed, rider] = build_id_data('00700.mat', outputs, inputs, 'pavilion/mat');
    guess = [9.7676, -1.5008, 3.0039, 0.7115, 0.1670, 43.3128];
    m = bicycle_grey(['Rigid' rider], speed, outputs, guess(1:5), guess(6));
    s = pem(z, m);
    parameters(i, :) = s.par';
    compare(z, s)
    saveas(gcf, ['plots/run700-' num2str(i) '.png'])
    close all
end
display(guess)
display(parameters)
