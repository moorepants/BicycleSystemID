#!/usr/bin/env python

# This is a demonstration of a very simple discrete black box system
# identification for a run.

import numpy as np
from matplotlib.pyplot import figure
import bicycledataprocessor as bdp

dataset = bdp.DataSet()
trial = bdp.Run(569, dataset, filterFreq=15.0)

states = np.vstack((trial.taskSignals['RollAngle'],
                    trial.taskSignals['SteerAngle'],
                    trial.taskSignals['RollRate'],
                    trial.taskSignals['SteerRate']))
inputs = np.vstack((trial.taskSignals['SteerTorque'],
                    trial.taskSignals['PullForce']))

numRows = states.shape[0] * states.shape[1]
# build B in Ax=B
# B is simply a vector of the "current" states, np x 1
b = states.T.reshape(numRows, 1)

# build A
numCols = states.shape[0]**2 + states.shape[0] * inputs.shape[0]
a = np.zeros((numRows - 1, numCols))

N = states.shape[1]
for i in range(N - 1):
    subCols = np.hstack((states[:, i - 1], inputs[:, i - 1]))
    for j in range(states.shape[0]):
        a[i * 4 + j, j * 6:j * 6 + len(subCols)] = subCols

# find the entries of A and B
x, resid, rank, s = np.linalg.lstsq(a, b[:-1])

A = np.hstack((x[:4], x[6:10], x[12:16], x[18:22])).T
B = np.hstack((x[4:6], x[10:12], x[16:18], x[22:])).T

# plot the results
xhat = np.dot(A, states) + np.dot(B, inputs)

fig = figure()
ax = fig.add_subplot(1, 1, 1)
time = trial.taskSignals['RollAngle'].time()
ax.plot(time,states.T, 'k.', time, xhat.T, 'b')
fig.show()
