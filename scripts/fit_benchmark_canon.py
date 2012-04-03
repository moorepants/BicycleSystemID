#!/usr/bin/env python
import cPickle
import canonical_system_id as csi

riders = ['Charlie', 'Jason', 'Luke']
environments = ['Horse Treadmill', 'Pavillion Floor']
maneuvers = ['Balance',
             'Track Straight Line',
             'Balance With Disturbance',
             'Track Straight Line With Disturbance']

runs = csi.select_runs(riders, maneuvers, environments)

# This gives the proportion of the lateral force which should be added to the
# steer torque and roll torque equations in the canonical equations.
H = csi.lateral_force_contribution(riders)

# try to load in all of the runs for the given run numbers.
trials, errors = csi.load_trials(runs, H)

rollParams = ['Mpd', 'C1pd', 'K0pd']
steerParams = ['Mdd', 'C1dp', 'C1dd',
               'K0dd', 'K2dd', 'HdF']
canon = csi.load_benchmark_canon(riders)

idMatrices = {}

runs = list(set(runs).difference(errors))
means = csi.mean_canon(riders, canon, H)
idMat = csi.enforce_symmetry(runs, trials, rollParams,
        steerParams, *means)
idMatrices['All'] = idMat

for rider in riders:
    idMatrices[rider] = {}

    # this find the model for each rider in both environments
    runs = csi.select_runs([rider], maneuvers, environments)
    runs = list(set(runs).difference(errors))

    means = csi.mean_canon(riders, canon, H)
    idMat = csi.enforce_symmetry(runs, trials, rollParams,
            steerParams, *means)
    idMatrices[rider]['All'] = idMat

    # this finds the model for each rider in each environment
    for env in environments:
        runs = csi.select_runs([rider], maneuvers, [env])
        runs = list(set(runs).difference(errors))

        means = csi.mean_canon(riders, canon, H)
        idMat = csi.enforce_symmetry(runs, trials, rollParams,
                steerParams, *means)
        idMatrices[rider][env] = idMat

# for each environment calculate the model for all riders
for env in environments:
    runs = csi.select_runs(riders, maneuvers, [env])
    runs = list(set(runs).difference(errors))

    means = csi.mean_canon(riders, canon, H)
    idMat = csi.enforce_symmetry(runs, trials, rollParams,
            steerParams, *means)
    idMatrices[env] = idMat

# save all of the identified matrices to file
with open('idMatrices.p', 'w') as f:
    cPickle.dump(idMatrices, f)

#print 'M'
#print M
#print M_id
#print abs(M_id - M) / M * 100
#print 'C1'
#print C1
#print C1_id
#print abs(C1_id - C1) / C1 * 100
#print 'K0'
#print K0
#print K0_id
#print abs(K0_id - K0) / K0 * 100
#print 'K2'
#print K2
#print K2_id
#print abs(K2_id - K2) / K2 * 100
#print 'H'
#print HMean
#print H_id
#print abs(H_id - HMean) / HMean * 100
