# This is an attempt to use the measured steer torque input from a run with no
# external disturbances and see if the Whipple model predicts the steer and
# roll.
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
from uncertainties import unumpy
import bicycleparameters as bp
from bicycledataprocessor import main as bdp

dataset = bdp.DataSet(fileName='../BicycleDataProcessor/InstrumentedBicycleData.h5')
dataset.open(mode='r')

pathToPar = '../BicycleParameters/data'

run = bdp.Run('00644', dataset.database, pathToPar, filterFreq=30.0)

lukeRigid = bp.Bicycle('Rigidcl', pathToData=pathToPar, forceRawCalc=True)
lukeRigid.add_rider('Luke', reCalc=True)

# roll rate, steer rate, roll angle, steer angle
# roll torque, steer torque
A, B = lukeRigid.state_space(4.0)

A = unumpy.nominal_values(A)
B = unumpy.nominal_values(B)
#C = np.eye(4)
#D = np.zeros((4, 2))
C = np.array([1., 0., 0., 0.])
D = np.array([0., 0.])

steerTorque = run.taskSignals['SteerTorque']
rollTorque = np.zeros_like(steerTorque)
u = np.vstack((rollTorque, steerTorque))

#whipple = signal.lti(A, B, C, D)

rollRate = run.taskSignals['RollRate']
steerRate = run.taskSignals['SteerRate']
rollAngle = run.taskSignals['RollAngle']
steerAngle = run.taskSignals['SteerAngle']

x0 = np.array([rollRate[0],
               steerRate[0],
               rollAngle[0],
               steerAngle[0]])

time = steerTorque.time()

#ysim = signal.lsim2(whipple, u, time, x0)
tsim, ysim, xsim = signal.lsim((A, B, C, D), u.T, time, X0=x0)

fig = plt.figure()
ax = fig.add_subplot(4, 1, 1)
ax.plot(time, ysim, time, rollRate)
ax = fig.add_subplot(4, 1, 2)
#ax.plot(time, ysim[1], time, steerRate)
#ax = fig.add_subplot(4, 1, 3)
#ax.plot(time, ysim[2], time, rollAngle)
#ax = fig.add_subplot(4, 1, 4)
#ax.plot(time, ysim[3], time, steerAngle)
fig.show()

