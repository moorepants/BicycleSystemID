#!/usr/bin/env python
import cPickle
import numpy as np
import matplotlib.pyplot as plt
import canonical_system_id as csi
from dtk import bicycle, control

allRiders = ['Charlie', 'Jason', 'Luke']
# create the comparison models
canon = csi.load_benchmark_canon(allRiders)
H = csi.lateral_force_contribution(allRiders)

with open('idMatrices.p') as f:
    idMatrices = cPickle.load(f)

def plot(canon, H, riders, environments, idMats):
    filename = ''
    for rider in riders:
        filename += '-' + rider
    for env in environments:
        filename += '-' + env.replace(' ', '')

    filename = 'canonical-id-plots/' + filename[1:]

    print filename

    v0 = 0.
    vf = 10.
    num = 100

    mM, mC1, mK0, mK2, mH = csi.mean_canon(riders, canon, H)
    speeds, mAs, mBs = bicycle.benchmark_state_space_vs_speed(mM, mC1, mK0, mK2,
            v0=v0, vf=vf, num=num)
    w, v = control.eigen_vs_parameter(mAs)
    mEigenvalues, mEigenvectors = control.sort_modes(w, v)

    iM, iC1, iK0, iK2, iH = idMats
    speeds, iAs, iBs = bicycle.benchmark_state_space_vs_speed(iM, iC1, iK0, iK2,
            v0=v0, vf=vf, num=num)
    w, v = control.eigen_vs_parameter(iAs)
    iEigenvalues, iEigenvectors = control.sort_modes(w, v)

    aAs, aBs, aSpeed = csi.mean_arm(riders)
    w, v = control.eigen_vs_parameter(aAs)
    aEigenvalues, aEigenvectors = control.sort_modes(w, v)

    rlfig = plt.figure()
    ax = rlfig.add_subplot(1, 1, 1)
    ax.plot(speeds, iEigenvalues.real, 'k-')
    ax.plot(speeds, abs(iEigenvalues.imag), 'k--')
    ax.plot(speeds, mEigenvalues.real, 'b-')
    ax.plot(speeds, abs(mEigenvalues.imag), 'b--')
    ax.plot(aSpeed, aEigenvalues.real, 'r-')
    ax.plot(aSpeed, abs(aEigenvalues.imag), 'r--')
    ax.set_ylim((-15, 10))
    ax.set_xlabel('Speed [m/s]')
    ax.set_ylabel('Real and Imaginary Parts of the Eigenvalues [1/s]')
    rlfig.savefig(filename + '-eig.png')
    rlfig.clf()

    # bode plots
    speeds = [2.0, 3.0, 5.8, 9.0]
    null, mAs, mBs = bicycle.benchmark_state_space_vs_speed(mM, mC1, mK0, mK2,
            speeds)
    null, iAs, iBs = bicycle.benchmark_state_space_vs_speed(iM, iC1, iK0, iK2,
            speeds)
    C = np.array([[1.0, 0.0, 0.0, 0.0]])
    D = np.array([[0.0, 0.0]])
    systems = []
    inputNames = ['$T_\phi$', '$T_\delta$']
    stateNames = ['$\phi$', '$\delta$', '$\dot{\phi}$', '$\dot{\delta}$']
    outputNames = ['$\phi$']
    for A, B in zip(mAs, mBs):
        systems.append(control.StateSpace(A, B, C, D, name='Whipple',
                stateNames=stateNames, inputNames=inputNames,
                outputNames=outputNames))

    for A, B in zip(iAs, iBs):
        systems.append(control.StateSpace(A, B, C, D, name='Canon ID',
                stateNames=stateNames, inputNames=inputNames,
                outputNames=outputNames))

    C = np.zeros((1, 19))
    C[0, 3] = 1
    D = 0.
    indices = [20, 30, 58, 90]
    for A, B in zip(aAs[indices], aBs[indices]):
        systems.append(control.StateSpace(A, B[:, [0, 2]], C, D, name='Arms',
                inputNames=inputNames, outputNames=outputNames))

    w = np.logspace(-1, 2)
    linestyles = ['-', '--', '-.', ':'] * 3
    colors = ['k'] * len(speeds) + ['b'] * len(speeds) + ['r'] * len(speeds)
    bode = control.Bode(w, *tuple(systems), linestyles=linestyles, colors=colors)
    bode.plot()
    bode.figs[0].savefig(filename + '-Tphi.png')
    bode.figs[1].savefig(filename + '-Tdel.png')
    bode.figs[0].clf()
    bode.figs[1].clf()

for k, v in idMatrices.items():
    if k == 'All':
        riders = ['Charlie', 'Jason', 'Luke']
        environments = ['Horse Treadmill', 'Pavillion Floor']
        plot(canon, H, riders, environments, v)
    elif k == 'Horse Treadmill' or k == 'Pavillion Floor':
        riders = ['Charlie', 'Jason', 'Luke']
        environments = [k]
        plot(canon, H, riders, environments, v)
    elif k in allRiders:
        riders = [k]
        for envName, envVal in v.items():
            if envName == 'All':
                environments = ['Horse Treadmill', 'Pavillion Floor']
                plot(canon, H, riders, environments, envVal)
            else:
                environments = [envName]
                plot(canon, H, riders, environments, envVal)
