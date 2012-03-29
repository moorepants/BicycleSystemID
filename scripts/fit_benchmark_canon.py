import numpy as np
import matplotlib.pyplot as plt
import bicycledataprocessor as bdp
import canonical_system_id as csi
from dtk import bicycle


riders = ['Charlie', 'Jason', 'Luke']
#riders = ['Luke']
environments = ['Horse Treadmill', 'Pavillion Floor']
maneuvers = ['Balance',
             'Track Straight Line',
             'Balance With Disturbance',
             'Track Straight Line With Disturbance']

runs = csi.select_runs(riders, maneuvers, environments)

# This gives the proportion of the lateral force which should be added to the
# steer torque and roll torque equations in the canonical equations.
H = {}
for rider in riders:
    A, B, F = csi.whipple_state_space(rider, 1.0)
    H[rider] = np.dot(np.linalg.inv(B[2:]), F[2:])

#rollParams = ['Mpp', 'Mpd', 'C1pp', 'C1pd',
               #'K0pp', 'K0pd', 'K2pp', 'K2pd'] # all
#rollParams = ['Mpp', 'Mpd', 'C1pd', 'K0pp', 'K0pd', 'K2pd'] # all but zeros
#rollParams = ['C1pd', 'K0pp', 'K0pd', 'K2pd'] # all but zeros and mass matrix
rollParams = ['Mpd', 'C1pd', 'K0pd']

trials = {}
errors = []
dataset = bdp.DataSet()

for i, r in enumerate(runs):
    try:
        trial = bdp.Run(r, dataset, filterFreq=15.)
    except bdp.bdpexceptions.TimeShiftError:
        errors.append(r)
    except IndexError:
        errors.append(r)
    else:
        trials[r] = trial

for runNum, trial in trials.items():
    timeSeries = csi.benchmark_time_series(trial, subtractMean=True)
    M, C1, K0, K2 = trial.bicycle.canonical(nominal=True)
    fixedValues = csi.benchmark_canon_to_dict(M, C1, K0, K2,
            H[trial.metadata['Rider']])
    rollA, rollB = csi.benchmark_lstsq_matrices(rollParams, timeSeries, fixedValues)

    try:
        totRollA = np.vstack((totRollA, rollA))
        totRollB = np.hstack((totRollB, rollB))
    except NameError:
        totRollA = rollA
        totRollB = rollB

rollSol = np.linalg.lstsq(totRollA, totRollB)[0]

#steerParams = ['Mdp', 'Mdd', 'C1dp', 'C1dd',
                #'K0dp', 'K0dd', 'K2dp', 'K2dd'] # all
#steerParams = ['Mdp', 'Mdd', 'C1dp', 'C1dd', 'K0dp', 'K0dd', 'K2dd'] # all but zeros
#steerParams = ['C1dp', 'C1dd', 'K0dp', 'K0dd', 'K2dd'] # all but zeros and mass matrix
#steerParams = ['Mdp', 'Mdd', 'C1dp', 'C1dd',
                #'K0dp', 'K0dd', 'K2dd']
steerParams = ['Mdd', 'C1dp', 'C1dd',
                'K0dd', 'K2dd', 'HdF']

for runNum, trial in trials.items():
    timeSeries = csi.benchmark_time_series(trial, subtractMean=True)
    M, C1, K0, K2 = trial.bicycle.canonical(nominal=True)
    M_mod = M.copy()
    K0_mod = K0.copy()
    M_mod[1, 0] = rollSol[rollParams.index('Mpd')]
    K0_mod[1, 0] = rollSol[rollParams.index('K0pd')]
    fixedValues = csi.benchmark_canon_to_dict(M_mod, C1, K0_mod, K2,
            H[trial.metadata['Rider']])
    #fixedValues['Mdp'] = rollSol[rollParams.index('Mpd')]
    #fixedValues['K0dp'] = rollSol[rollParams.index('K0pd')]
    steerA, steerB = csi.benchmark_lstsq_matrices(steerParams, timeSeries, fixedValues)
    try:
        totSteerA = np.vstack((totSteerA, steerA))
        totSteerB = np.hstack((totSteerB, steerB))
    except NameError:
        totSteerA = steerA
        totSteerB = steerB

steerSol = np.linalg.lstsq(totSteerA, totSteerB)[0]

canon = (M, C1, K0, K2)
M_id, C1_id, K0_id, K2_id, H_id = csi.benchmark_identified_matrices((rollParams,
    steerParams), (rollSol, steerSol), M_mod, C1, K0_mod, K2, H[trial.metadata['Rider']])

speeds = np.linspace(0, 10, num=50)
eVals = np.zeros((len(speeds), 4), dtype=np.complex)
for i, v in enumerate(speeds):
        A, B = bicycle.abMatrix(M_id, C1_id, K0_id, K2_id, v, 9.8)
        eVals[i] = np.linalg.eig(A)[0]

eVals2 = np.zeros((len(speeds), 4), dtype=np.complex)
for i, v in enumerate(speeds):
    A, B = bicycle.abMatrix(M, C1, K0, K2, v, 9.8)
    eVals2[i] = np.linalg.eig(A)[0]

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
ax.plot(speeds, eVals2, '.k', speeds, eVals, '.b')
ax.set_ylim((-10, 10))

print 'M'
print M
print M_id
print 'C1'
print C1
print C1_id
print 'K0'
print K0
print K0_id
print 'K2'
print K2
print K2_id
print 'H'
print H[trial.metadata['Rider']]
print H_id
