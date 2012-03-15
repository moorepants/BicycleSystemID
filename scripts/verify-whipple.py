#!/usr/bin/env python

# This script generates the time history of the steer and roll torques given
# the steer/roll coordinates, rates and accelerations for both the Whipple
# model and the Whipple model with the inertial effects of the arms and
# compares it to the measured steer torque.

import numpy as np
import matplotlib.pyplot as plt
import bicycledataprocessor as bdp

# load the experimental data
dataset = bdp.DataSet(fileName='../../BicycleDataProcessor/InstrumentedBicycleData.h5')

pathToPar = '../../BicycleParameters/data'

trial = bdp.Run('00638', dataset, pathToPar, filterFreq=10.0)

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

# this is the model that I identified for this run (except for the roll torque
# entries in the B matrix, they are copied from the Whipple model)
# I'm not sure it is a good idea to use the roll entries from the Whipple
# model, it might introduce some unwanted coupling. The resulting trace doesn't
# look so good.
idA = np.array([[8.5007, -18.5512, 0.2315, -2.2580],
                [40.3916, -50.5260, -2.5057, -11.2178]])
idB = np.array([[0.00942706974495, 0.1652],
                 [-0.0971149781346, 1.4777]])
idU = np.dot(np.linalg.inv(idB), (xd - np.dot(idA, x)))

# this is a model derived from the identified data using a mixed effects
# regression
meA = np.array([[7.530644, 1.3192618 - 0.8444925 * speed**2,
                 0.08954696 - 0.20096398 * speed, -3.0624926 + 0.5778556 * speed],
                [71.20166, 9.538313 - 3.196099 * speed**2,
                 4.766327 - 1.508476 * speed, -6.7510187  -0.4292608 * speed]])
meB = np.array([[ 0.00942706974495, 0.07229608],
                [ -0.0971149781346, 1.727447]])
meU = np.dot(np.linalg.inv(meB), (xd - np.dot(meA, x)))

# plot the results
fig = plt.figure()

time = steerTorque.time()

ax1 = fig.add_subplot(2, 1, 1)
ax1.plot(time, u[0, :],
         time, armU[0, :],
         time, idU[0, :],
         time, meU[0, :])
ax1.set_ylabel('Roll Torque [Nm]')
ax1.set_xlabel('Time [s]')
ax1.legend(('Whipple Model Prediction', 'Arm Model Prediction', 'Identified', 'Mixed Effects'))
leg1 = ax1.get_legend()
plt.setp(leg1.get_texts(), fontsize='xx-small')

ax2 = fig.add_subplot(2, 1, 2)
ax2.plot(time, steerTorque, 'k',
         time, u[1, :],
         time, armU[1, :],
         time, idU[1, :],
         time, meU[1, :])
ax2.set_ylabel('Steer Torque [Nm]')
ax2.set_xlabel('Time [s]')
ax2.legend(('Experimental Measurment', 'Whipple Model Prediction',
    'Arm Model Prediction', 'Identified Model', 'Mixed Effects'))
leg2 = ax2.get_legend()
plt.setp(leg2.get_texts(), fontsize='xx-small')

fig.savefig('../plots/torque-comparison.png')
fig.show()
