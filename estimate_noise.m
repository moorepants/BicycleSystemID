config

load([PATH_TO_RUN_MAT_DIRECTORY '00246.mat'], ...
     'YawAngle', ...
     'RollAngle', ...
     'SteerAngle', ...
     'YawRate', ...
     'RollRate', ...
     'SteerRate', ...
     'LateralRearContact', ...
     'SteerTorque')

start = 2000;
stop = 16000;
var(SteerAngle(start:stop))
sigma = std(SteerAngle(start:stop) - mean(SteerAngle(start:stop)))
hold all
plot(SteerAngle - mean(SteerAngle))
plot(start:stop, sigma * ones(stop-start+1, 1), 'k')
plot(start:stop, -sigma * ones(stop-start+1, 1), 'k')