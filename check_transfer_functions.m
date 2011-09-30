function check_transfer_functions(bicycle, speed, parameters)

tfOutputs = {'Delta', 'PhiDot', 'Tdelta'}

data = generate_data(bicycle, speed, ...
                     'simulate', 0, ...
                     'loopTransfer', 0, ...
                     'handlingQuality', 0, ...
                     'gains', parameters(1:5), ...
                     'neuroFreq', parameters(6), ...
                     'forceTransfer', tfOutputs);

% load the signals from the experiment
load('00105.mat', 'PullForce', 'SteerTorque', 'SteerAngle', 'RollRate', 'RollAngle', 'YawAngle')

% link the data to the variable name
outputs.PhiDot = RollRate;
outputs.Delta = SteerAngle;
outputs.Tdelta = SteerTorque;

regressors.PhiDot = [6, 7, 0];
regressors.Delta = [7, 6, -1];
regressors.Tdelta = [6, 2, -3];

time = linspace(0, (length(PullForce) - 1) / 200, length(PullForce));
w = logspace(0, 1.3, 200);

for i = 1:length(tfOutputs)
    % compute the transfer function from the experimental data
    [expTransFuncs.(tfOutputs{i}), ~, ~, ~] = ...
        gettf1(PullForce, outputs.(tfOutputs{i}), ...
               regressors.(tfOutputs{i}), time, 0);
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
    p.Title.String = tfOutputs{i};
    setoptions(b2, p);
    hold off
end
