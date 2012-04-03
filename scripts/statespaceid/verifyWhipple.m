config
addpath(PATH_TO_CONTROL_MODEL)

runData = load([PATH_TO_RUN_MAT_DIRECTORY filesep '00644.mat']);

u = runData.SteerTorque;
y = [runData.RollAngle',
     runData.SteerAngle',
     runData.RollRate',
     runData.SteerRate'];

z = iddata(y', u, 1 / 200);

set(z, 'InputName', {'Steer Torque'})
set(z, 'InputUnit', {'N-m'})
set(z, 'OutputName', {'Roll Angle', 'Steer Angle', 'Roll Rate', 'Steer Rate'})
set(z, 'OutputUnit', {'rad', 'rad', 'rad/s', 'rad/s'})

pathToFile = [PATH_TO_CONTROL_MODEL filesep 'parameters' filesep 'RigidLukePar.txt']
par = par_text_to_struct(pathToFile);
[A, B, C, D] = whipple_pull_force_abcd(par, 4.0);

A([1, 2, 3, 5, 6, 8, 10], :) = [];
A(:, [1, 2, 3, 5, 6, 8, 10]) = [];

B([1, 2, 3, 5, 6, 8, 10], :) = [];
B(:, [1, 3]) = [];

C([1, 2, 3, 5, 6, 8, 9, 10, 11, 13, 14, 16, 17, 18], :) = [];
C(:, [1, 2, 3, 5, 6, 8, 10]) = [];

D = 0;

con = ss(A, B, C, D);
set(con, 'InputName', {'Steer Torque'})
set(con, 'OutputName', {'Roll Angle', 'Steer Angle', 'Roll Rate', 'Steer Rate'})
time = linspace(0, (length(u) - 1) / 200, length(u))';
figure()
lsim(con, u, time)%, y(:, 1)');

bicycle = idss(con);
set(bicycle, 'InputName', {'Steer Torque'})
set(bicycle, 'InputUnit', {'N-m'})
set(bicycle, 'OutputName', {'Roll Angle', 'Steer Angle', 'Roll Rate', 'Steer Rate'})
set(bicycle, 'OutputUnit', {'rad', 'rad', 'rad/s', 'rad/s'})

figure()
[ysim, fit, x0sim] = compare(z, bicycle); % compare estimates x0 by default

set(bicycle, 'InitialState', 'Fixed', 'x0', y(:, 1))
% this makes the model unstable if we use the given initial condition from
% the data, or if you use zero
%compare(z, bicycle, 'InitialState', 'm')
