function [dx, y] = parameterized_bicycle(t, x, u, kDelta, kPhiDot, kPhi, kPsi, kY, fnm, varargin)

gains = [kDelta, kPhiDot, kPhi, kPsi, kY];

data = generate_data(bicycle, speed,
                     'gains', gains,
                     'neuroFreq', fnm,
                     'simulate', 0,
                     'loopTransfer', 0,
                     'handlingQuality', 0,
                     'forceTransfer', 0,
                     'stateSpace', {A, B, C, D});

A = data.system.A;
B = data.system.B;
C = data.system.C;
D = data.system.D;
outputs = data.system.outputs;

dx = A * x + B * u;
y = C * x + D * u;
