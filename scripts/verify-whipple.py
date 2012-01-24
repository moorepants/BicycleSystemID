#!/usr/bin/env python

# This script generates the time history of the steer and roll torques given
# the steer/roll coordinates, rates and accelerations for both the Whipple
# model and the Whipple model with the inertial effects of the arms and
# compares it to the measured steer torque.

import numpy as np
import matplotlib.pyplot as plt
from bicycledataprocessor import main as bdp

# load the experimental data
dataset = bdp.DataSet(fileName='../../BicycleDataProcessor/InstrumentedBicycleData.h5')
dataset.open(mode='r')

pathToPar = '../../BicycleParameters/data'

trial = bdp.Run('00638', dataset.database, pathToPar, filterFreq=30.0)

speed = trial.taskSignals['ForwardSpeed'][50:-50].mean()

steerTorque = trial.taskSignals['SteerTorque'].subtract_mean()

rollAngle = trial.taskSignals['RollAngle'].subtract_mean()
steerAngle = trial.taskSignals['SteerAngle'].subtract_mean()
q = np.vstack((rollAngle, steerAngle))

rollRate = trial.taskSignals['RollRate'].subtract_mean()
steerRate = trial.taskSignals['SteerRate'].subtract_mean()
qd = np.vstack((rollRate, steerRate))

rollAccel = rollRate.time_derivative().filter(30.).subtract_mean()
steerAccel = steerRate.time_derivative().filter(30.).subtract_mean()
qdd = np.vstack((rollAccel, steerAccel))

# compute the inputs from the Whipple model
lukeRigid = trial.bicycle

M, C1, K0, K2 = lukeRigid.canonical(nominal=True)

C = speed * C1
K = K0 * trial.bicycleRiderParameters['g'] + K2 * speed ** 2

u = np.dot(M, qdd) + np.dot(C, qd) + np.dot(K, q)

# compute the inputs from the Whipple model with rider arms
# this is only valid for run 638 (i.e. Luke on Rigidcl at the measured speed)

armA = np.array([[8.7171, -18.8791, -0.0370, -1.4646], # phi
                 [4.3115, -1.5062, 2.4852, -7.0466]]) # delta
# rows: phi, delta, columns: Tphi, Tdelta
armB = np.array([[0.0081, -0.0241],
                 [-0.0241, 1.6196]])
xd = qdd
x = np.vstack((q, qd))

armU = np.dot(np.linalg.inv(armB), (xd - np.dot(armA, x)))

# plot the results
fig = plt.figure()

time = steerTorque.time()

ax1 = fig.add_subplot(2, 1, 1)
ax1.plot(time, u[0, :],
         time, armU[0, :])
ax1.set_ylabel('Roll Torque [Nm]')
ax1.set_xlabel('Time [s]')
ax1.legend(('Whipple Model Prediction', 'Arm Model Prediction'))
leg1 = ax1.get_legend()
plt.setp(leg1.get_texts(), fontsize='small')

ax2 = fig.add_subplot(2, 1, 2)
ax2.plot(time, steerTorque, 'k',
         time, u[1, :],
         time, armU[1, :])
ax2.set_ylabel('Steer Torque [Nm]')
ax2.set_xlabel('Time [s]')
ax2.legend(('Experimental Measurment', 'Whipple Model Prediction',
    'Arm Model Prediction'))
leg2 = ax2.get_legend()
plt.setp(leg2.get_texts(), fontsize='x-small')

fig.savefig('../plots/torque-comparison.png')
fig.show()
