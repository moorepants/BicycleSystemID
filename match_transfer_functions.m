function bestParameters = match_transfer_functions(bicycle, speed)

% get the initial guess for the gains using Ron's technique and the state
% space model
%guess = [39.3, -0.018, 209.895, 0.081, 0.799, 37];
%data = generate_data(bicycle, speed, ...
                     %'simulate', 0, ...
                     %'loopTransfer', 0, ...
                     %'handlingQuality', 0, ...
                     %'gains', guess(1:5), ...
                     %'neuroFreq', guess(6), ...
                     %'forceTransfer', {'Delta', 'PhiDot', 'Tdelta'});

data = generate_data(bicycle, speed, ...
                     'simulate', 0, ...
                     'loopTransfer', 0, ...
                     'handlingQuality', 0, ...
                     'forceTransfer', {'Delta', 'PhiDot', 'Tdelta'});

% the initial guess for gains
guess = [data.modelPar.kDelta,
         data.modelPar.kPhiDot,
         data.modelPar.kPhi,
         data.modelPar.kPsi,
         data.modelPar.kY,
         sqrt(data.modelPar.neuroNum)];

% get the state space model so we don't have to calculate it during the gain
% matching algorithm
stateSpace = {data.modelPar.A,
              data.modelPar.B,
              data.modelPar.C,
              data.modelPar.D};

load('00105.mat', 'PullForce', 'SteerTorque', 'SteerAngle', 'RollRate')

outputs.PhiDot = RollRate;
outputs.Delta = SteerAngle;
outputs.Tdelta = SteerTorque;

regressors.PhiDot = [6, 7, 0];
regressors.Delta = [7, 6, -1];
regressors.Tdelta = [6, 2, -3];

time = linspace(0, (length(PullForce) - 1) / 200, length(PullForce));
w = logspace(0, 1.3, 200);

tfOutputs = {'Tdelta'};

for i = 1:length(tfOutputs)
    % compute the transfer function from the experimental data
    [expTransFuncs.(tfOutputs{i}), ~, ~, ~] = ...
        gettf1(PullForce, outputs.(tfOutputs{i}), ...
               regressors.(tfOutputs{i}), time, 0);

    % find the parameters that force the model to best fit the experimental data
    bestParameters.(tfOutputs{i}) = match_gains(bicycle, speed, ...
        stateSpace, tfOutputs{i}, expTransFuncs.(tfOutputs{i}), guess)

    display(sprintf(['The lateral force to %s transfer function has a ' ...
                     'best match with these parameters:'], tfOutputs{i}))
    % now find the transfer function with the best guess
    bestPar = bestParameters.(tfOutputs{i})
    % generate the model with the best parameters
    data = generate_data(bicycle, speed, ...
                         'gains', bestPar(1:5), ...
                         'neuroFreq', bestPar(6), ...
                         'simulate', 0, ...
                         'loopTransfer', 0, ...
                         'handlingQuality', 0, ...
                         'stateSpace', stateSpace);
    % plot the fit
    figure(i)
    b1 = bodeplot(expTransFuncs.(tfOutputs{i}), w);
    p = getoptions(b1);
    p.PhaseMatching = 'on';
    p.PhaseMatchingFreq = 1;
    p.PhaseMatchingValue = 0;
    setoptions(b1, p);
    hold all
    modelTF = tf(data.forceTF.(tfOutputs{i}).num, ...
                 data.forceTF.(tfOutputs{i}).den);
    b2 = bodeplot(modelTF, w);
    p = getoptions(b2);
    p.PhaseMatching = 'on';
    p.PhaseMatchingFreq = 1;
    p.PhaseMatchingValue = 0;
    setoptions(b2, p);
    hold off
end
