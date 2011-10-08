function expTransFuncs = sysid_toolbox();

globals

pathToData = PATH_TO_RUN_MAT_DIRECTORY;

load([pathToData '00264.mat'], 'RollRate', 'SteerRate', 'RollAngle', ...
    'SteerAngle', 'PullForce', 'SteerTorque', 'YawRate', 'YawAngle', ...
    'LateralRearContact')

% build the iddata structure
u = PullForce;

y = [YawAngle, RollAngle, SteerAngle, YawRate, RollRate, SteerRate, ...
    LateralRearContact, SteerTorque];

outputNames = {'Yaw Angle', ...
               'Roll Angle', ...
               'Steer Angle', ...
               'Yaw Rate', ...
               'Roll Rate', ...
               'Steer Rate', ...
               'Lateral Rear Contact', ...
               'Steer Torque'};

outputVariables = {'Psi', 'Phi', 'Delta', 'PsiDot', 'PhiDot', 'DeltaDot', 'Y', 'Tdelta'};

outputUnits = {'Radian', ...
               'Radian', ...
               'Radian', ...
               'Radian/Second', ...
               'Radian/Second', ...
               'Radian/Second', ...
               'Meter', ...
               'Newton-Meter'};

sampleSize = 1 / 200;

z = iddata(y, u, sampleSize);

set(z, 'InputName', 'Lateral Force')
set(z, 'InputUnit', 'Newtons')

set(z, 'OutputName', outputNames)
set(z, 'OutputUnit', outputUnits)

% id on a subset of the run
ze = z(2000:4000);

% find a general model
m = pem(ze);

% get the discrete transfer functions
withoutNoise = ss(m.A, m.B, m.C, m.D, sampleSize);
% get the continous tranfer functions of the lateral force
con = d2c(withoutNoise, 'tustin');

[b, a] = ss2tf(con.a, con.b, con.c, con.d, 1);

for i = 1:length(outputNames)
    expTransFuncs.(outputVariables{i}) = tf(b(i, :), a);
end
