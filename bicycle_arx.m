% go through each input/output pair and get the best parameters for an
% arx model
%naVals = zeros(size(y, 2), 1);
%nbVals = zeros(size(y, 2), 1);
%nkVals = zeros(size(y, 2), 1);
%
%for i = 1:size(y, 2);
    %zs = iddata(y(:, i), u, 1 / 200);
    %set(zs, 'OutputName', outputNames{i}, 'OutputUnit', outputUnits{i})
    %set(zs, 'InputName', 'Lateral Force', 'InputUnit', 'Newtons')
    %naGuess = 6:10;
    %nbGuess = 6:10;
    %nkGuess = delayest(zs(2000:4000));
    %trials = struc(naGuess, nbGuess, nkGuess);
    %arxPar = selstruc(arxstruc(zs(2000:4000), zs, trials), 0);
    %arxMod = arx(zs, arxPar);
    %bode(arxMod)
    %pause
    %naVals(i) = arxPar(1);
    %nbVals(i) = arxPar(2);
    %nkVals(i) = arxPar(3);
%end

%naVals = 10 * ones(size(y, 2))
%nbVals
%nkVals
%
%arxModel = arx(z, 'na', naVals, 'nb', nbVals, 'nk', nbVals)
