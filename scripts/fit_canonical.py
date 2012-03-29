import numpy as np
import matplotlib.pyplot as plt
import bicycledataprocessor as bdp
import canonical_system_id as csi

# This gives the proportion of the lateral force which should be added to the
# steer torque and roll torque equations in the canonical equations.
F = {}
for rider in ['Charlie', 'Jason', 'Luke']:
    F[rider] = csi.whipple_state_space(rider, 1.0)[2][2:]

# find the runs that we want to id
dataset = bdp.DataSet()
dataset.open()

table = dataset.database.root.runTable

runs = []
for row in table.iterrows():
    con = []
    con.append(row['Rider'] in ['Jason', 'Charlie', 'Luke'])
    con.append(row['Maneuver'] in ['Balance',
                                   'Track Straight Line',
                                   'Balance With Disturbance',
                                   'Track Straight Line With Disturbance'])
    con.append(row['Environment'] == 'Horse Treadmill')
    con.append(row['corrupt'] is not True)
    con.append(int(row['RunID']) > 100)
    if False not in con:
        runs.append(row['RunID'])

dataset.close()

idMassMats = np.zeros((len(runs), 2, 2))
idDampMats = np.zeros((len(runs), 2, 2))
idStifMats = np.zeros((len(runs), 2, 2))
speeds = np.nan * np.ones(len(runs))

thetaDelta = ['m21', 'm22', 'c21', 'c22', 'k21', 'k22']

errors = []
for i, r in enumerate(runs):
    try:
        trial = bdp.Run(r, dataset, filterFreq=15.)
    except bdp.bdpexceptions.TimeShiftError:
        errors.append(r)
    except IndexError:
        errors.append(r)
    else:
        if trial.metadata['Maneuver'].endswith('Disturbance'):
            thetaPhi = ['m11', 'm12', 'c11', 'c12', 'k11', 'k12']
        else:
            thetaPhi = ['c11', 'c12', 'k11', 'k12']

        v = trial.taskSignals['ForwardSpeed'].mean()
        speeds[i] = v
        g = trial.bicycleRiderParameters['g']

        M, C1, K0, K2 = trial.bicycle.canonical(nominal=True)
        C = C1 * v
        K = K0 * g + K2 * v**2
        canon = (M, C, K)

        timeSeries = csi.time_series(trial, F)
        M_id, C_id, K_id = csi.compute_unknowns(thetaPhi, thetaDelta,
                timeSeries, canon)
        idMassMats[i] = M_id
        idDampMats[i] = C_id
        idStifMats[i] = K_id

        #forces_id = np.dot(M_id, accels) + np.dot(C_id, rates) + np.dot(K_id,
                #coordinates)
#
        #time = trial.taskSignals['ForwardSpeed'].time()

        #fig = plt.figure()
        #for i in range(2):
            #ax = fig.add_subplot(2, 1, i + 1)
            #ax.plot(time, forces[i], time, forces_id[i])
            #ax.legend(('Experimental', 'Identified'))
        #fig.show()

fig = plt.figure()
for i in range(2):
    ax = fig.add_subplot(2, 6, 1 + i * 6)
    ax.plot(speeds, idMassMats[:, i, 0], '.')

    ax = fig.add_subplot(2, 6, 2 + i * 6)
    ax.plot(speeds, idMassMats[:, i, 1], '.')

    ax = fig.add_subplot(2, 6, 3 + i * 6)
    ax.plot(speeds, idDampMats[:, i, 0], '.')

    ax = fig.add_subplot(2, 6, 4 + i * 6)
    ax.plot(speeds, idDampMats[:, i, 1], '.')

    ax = fig.add_subplot(2, 6, 5 + i * 6)
    ax.plot(speeds, idStifMats[:, i, 0], '.')

    ax = fig.add_subplot(2, 6, 6 + i * 6)
    ax.plot(speeds, idStifMats[:, i, 1], '.')
fig.show()
