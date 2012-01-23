% only fit the minimal states of the Whipple model
states = {'phi', 'delta', 'phiDot', 'deltaDot'};
inputs = {'tDelta'};
outputs = {'phi', 'delta', 'phiDot', 'deltaDot'};

% load the data
[data, v] = build_id_data('00638.mat', outputs, inputs, '');

% build a structured idss model
whippleModel = bicycle_structured('RigidLuke', v, 'states', states, ...
    'inputs', inputs, 'outputs', outputs);

% identify the free entries of the structured model
identifiedModel = pem(data, whippleModel, 'InitialState', 'Estimate');

% this is from my first attempt to calculate the arms model, i think it is
% correct
armsA = [0         0    1.0000         0
         0         0         0    1.0000
         8.7171  -18.6499   -0.0368   -1.4557
         4.3115   -1.3594    2.4701   -7.0037];

armsB = [0
         0
         -0.1019
         5.5687];

armsModel = whippleModel;
armsModel.A = armsA;
armsModel.B = armsB;

% plot the comparison of the models
compare(data, whippleModel, identifiedModel, armsModel);
