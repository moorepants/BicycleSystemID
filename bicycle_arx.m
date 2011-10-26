% go through each input/output pair and get the best parameters for an
% arx model

outputs = {'yP', 'psi', 'phi', 'delta', 'yPDot', 'psiDot', 'phiDot', 'deltaDot', 'tDelta'};

nanbnk = zeros(length(outputs), 3);

for i = 1:length(outputs)
    z = build_id_data('00105.mat', outputs(i));
    naGuess = 6:10;
    nbGuess = 6:10;
    nkGuess = delayest(z(200 * 8.26:200 * 44.48));
    trials = struc(naGuess, nbGuess, nkGuess);
    arxPar = selstruc(arxstruc(z(200 * 8.26:200 * 44.48), z(200 * 8.26:200 * 44.48), trials), 0)
    nanbnk(i, :) = arxPar;
    %arxMod = arx(z, arxPar);
    %bode(arxMod)
    %pause
end

%naVals = 10 * ones(size(y, 2))
%nbVals
%nkVals
%
%arxModel = arx(z, 'na', naVals, 'nb', nbVals, 'nk', nbVals)
