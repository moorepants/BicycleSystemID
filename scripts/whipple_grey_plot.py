import os

import numpy as np
from scipy.io import loadmat
import matplotlib.pyplot as plt

# debugging
from IPython.core.debugger import Tracer
set_trace = Tracer()

import bicycledataprocessor.main as bdp
from bicycledataprocessor.database import get_row_num
import bicycleparameters as bp

a = loadmat('whipple_grey_fit.mat', squeeze_me=True)

runs = [os.path.splitext(str(r))[0] for r in a['matFiles']]

from config import PATH_TO_DATABASE, PATH_TO_H5, PATH_TO_CORRUPT, PATH_TO_PARAMETERS

runNums = np.array([int(x) for x in runs])

thresholds = [0, 20, 40, 60, 70]

def plot_coefficients(fitThreshold):

    meanFits = np.mean(a['fits'], 1)
    print(runNums[np.nanargmax(meanFits)])
    indices = np.nonzero(meanFits > fitThreshold)
    stateMats = a['stateMatrices'][indices[0], :, :]
    inputMats = a['inputMatrices'][indices[0], :]
    speeds = a['speeds'][indices[0]]

    dataset = bdp.DataSet(fileName=PATH_TO_DATABASE, pathToH5=PATH_TO_H5,
            pathToCorruption=PATH_TO_CORRUPT)
    dataset.open()

    # find all the runs for Jason that were without disturbances
    table = dataset.database.root.runTable

    riders = []
    maneuvers = []
    for r in runNums[indices[0]]:
        i = get_row_num(r, table)
        riders.append(table[i]['Rider'])
        maneuvers.append(table[i]['Maneuver'])
    riders = np.array(riders)

    dataset.close()

    modelSpeeds = np.linspace(0., 10., num=50)
    modelStateMats = {}
    modelInputMats = {}

    for rider in set(riders):
        modelStateMats[rider] = np.zeros((len(modelSpeeds), 4, 4))
        modelInputMats[rider] = np.zeros((len(modelSpeeds), 4))
        # create a Bicycle object for the rider/bicycle
        if rider == 'Jason':
            bicycle = 'Rigid'
        else:
            bicycle = 'Rigidcl'
        bicycle = bp.Bicycle(bicycle, pathToData=PATH_TO_PARAMETERS, forceRawCalc=True)
        bicycle.add_rider(rider)
        # compute the A and B matrices as a function of speed for the Whipple
        # model
        for i, v in enumerate(modelSpeeds):
            A, B = bicycle.state_space(v, nominal=True)
            modelStateMats[rider][i, :, :] = A
            # the second column is for steer torque
            modelInputMats[rider][i, :] = B[:, 1]

    fig = plt.figure()

    equations = [r'\dot{\phi}', r'\dot{\delta}', r'\ddot{\phi}', r'\ddot{\delta}']
    states = [r'\phi', r'\delta', r'\dot{\phi}', r'\dot{\delta}']

    aXlims = np.array([[0., 20.],
                       [-80., 0.],
                       [-1.5, 4.],
                       [-8., 0.],
                       [-0.2, 1.2],
                       [0., 200.],
                       [-175., 40.],
                       [-40., 60.],
                       [-100., 0.],
                       [0., 20.]])

    onlyRiders = set(riders)

    for i in range(2, 4):
        for j in range(4):
            # make a row for for the roll acceleration and steer acceleration
            # equations and a colum for the four states and one for the input
            ax = fig.add_subplot(2, 5, 5 * (i - 2) + j + 1)
            ax.set_title('$a_{' + equations[i] + states[j] + '}$')

            for rider in onlyRiders:

                modelLine = ax.plot(modelSpeeds, modelStateMats[rider][:, i, j])

                riderIndices = np.nonzero(riders == rider)
                expLine = ax.plot(speeds[riderIndices[0]], stateMats[riderIndices[0], i, j], '.')
                expLine[0].set_color(modelLine[0].get_color())
                print(rider, modelLine[0].get_color())

            ax.set_ylim(aXlims[5 * (i - 2) + j])

    for i, p in zip(range(2, 4), [5, 10]):
        ax = fig.add_subplot(2, 5, p)
        ax.set_title('$b_{' + equations[i] + r'T_\delta' + '}$')
        for rider in onlyRiders:
            riderIndices = np.nonzero(riders == rider)
            modelLine = ax.plot(modelSpeeds, modelInputMats[rider][:, i])
            expLine = ax.plot(speeds[riderIndices[0]], inputMats[riderIndices[0], i], '.')
            expLine[0].set_color(modelLine[0].get_color())
        ax.set_ylim(aXlims[p - 1])

    fig.suptitle('n = {}, Mean fit threshold: {}%'.format(len(speeds), fitThreshold))

    return fig

for thresh in thresholds:
    fig = plot_coefficients(thresh)
    fig.savefig('whipple-thresh-' + str(thresh) + '.png')

