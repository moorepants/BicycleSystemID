
allOutputs = {'phiDot', 'delta', 'psiDot', 'phi', 'deltaDot', 'tDelta', 'yP'};
parameters = zeros(length(allOutputs), 6);
for i = 1:length(allOutputs)
    outputs = allOutputs(1:i);
    [z, speed, rider] = build_id_data('00700.mat', outputs, 'pavilion/mat');
    %guess = [9.7676, -1.5008, 3.0039, 0.7115, 0.1670, 43.3128];
    guess = [15.0387, -0.8800, 4.2130, 0.4739, 0.2715, 43.3185];
    m = bicycle_grey(['Rigid' rider], speed, outputs, guess(1:5), guess(6));
    grey = pem(z, m, 'DisturbanceModel', 'Estimate', 'InitialState', 'Estimate');
    greyKMatrices.(['K' num2str(i)]) = grey.K;
    parameters(i, :) = grey.par';
    black = pem(z, 8, 'DisturbanceModel', 'Estimate', 'InitialState', 'Estimate');
    blackKMatrices.(['K' num2str(i)]) = black.K;
    compare(z, grey, black)
    saveas(gcf, ['plots/run700-estimate-k-x0-' num2str(i) '.png'])
    close all
end
display(guess)
display(parameters)
