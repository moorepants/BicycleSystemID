config
addpath(PATH_TO_CONTROL_MODEL)
v = 5.0;
bicycle = bicycle_state_space('Benchmark', v);
% This is the model identified from Luke's pavilion runs with the canonical
% id method.
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
g = 9.81;
newBicycle = replace_essential(bicycle, M, C1, K0, K2, H, v, g);
