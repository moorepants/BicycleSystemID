addpath('..')

runFile = '../../../data/riderid/disturbance-runs/00699.mat';

gains = [9.7676, -1.5008, 3.0039, 0.7115, 0.1670];
neuro = [43.3128, 0.707];

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

gains = [50, -2.5008, 3.0039, 0.7115, 0.1670];
neuro = [43.3128, 0.707];

gains = [45.4460, -0.5284, 4.7458, 0.3923, 0.1670];
neuro = [27.3638, 0.707];

% 2.0 m/s gains from Ron's method
%gains = [28.0000, -1.1300, 3.7908, 0.0932, 0.2237];
%neuro = [30, 0.707];

% 9.0 m/s gains from Ron's method
%gains = [39.4000   -0.4250    2.4514    1.1023    0.0568];
%neuro = [30, 0.707];

[parameters, fit, speed] = id_rider_siso(runFile, 'delta', gains, ...
    neuro, false, {M, C1, K0, K2, H}, 'FixedParameter', 'zetanm', ...
    'Focus', 'stability')

whipple = bicycle_state_space(['RigidLuke'], speed);
augmented = replace_essential(whipple, M, C1, K0, K2, H, speed, 9.81);

%hess = generate_data('RigidLuke', speed, ...
                     %'gains', parameters(1:5), ...
                     %'neuroFreq', parameters(6), ...
                     %'plot', 1, ...
                     %'stateSpace', {augmented.A, augmented.B, augmented.C, augmented.D});

outputs = {'phi', 'delta', 'phiDot' 'deltaDot', 'psi', 'yQ', 'tDelta'};
grey = bicycle_grey('lateral', augmented, speed, ...
    outputs, parameters(1:5)', parameters(6:7)');

[directory, fileName, ext] = fileparts(runFile);
[data, meta] = build_id_data([fileName ext], outputs, {'fB'}, directory, true);

% see if results improve when fitting to all the outputs
all = pem(data, grey, 'FixedParameter', 'zetanm');

figure()
compare(data, grey, all)
