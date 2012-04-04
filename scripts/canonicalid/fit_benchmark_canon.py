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
covarMatrices = {}

print('Computing the estimate for all runs.')
runs = list(set(runs).difference(errors))
means = csi.mean_canon(riders, canon, H)
idMat, rollCovar, steerCovar = csi.enforce_symmetry(runs, trials, rollParams,
        steerParams, *means)
idMatrices['All'] = idMat
covarMatrices['All'] = (rollCovar, steerCovar)
print('Done.')

for rider in riders:
    idMatrices[rider] = {}
    covarMatrices[rider] = {}

    # this find the model for each rider in both environments
    print('Computing the estimate for {} in both environments.'.format(rider))
    runs = csi.select_runs([rider], maneuvers, environments)
    runs = list(set(runs).difference(errors))

    means = csi.mean_canon(riders, canon, H)
    idMat, rollCovar, steerCovar = csi.enforce_symmetry(runs, trials,
            rollParams, steerParams, *means)
    idMatrices[rider]['All'] = idMat
    covarMatrices[rider]['All'] = (rollCovar, steerCovar)
    print('Done.')

    # this finds the model for each rider in each environment
    for env in environments:
        print('Computing the estimate for {} on {}.'.format(rider, env))
        runs = csi.select_runs([rider], maneuvers, [env])
        runs = list(set(runs).difference(errors))

        means = csi.mean_canon(riders, canon, H)
        idMat, rollCovar, steerCovar = csi.enforce_symmetry(runs, trials,
                rollParams, steerParams, *means)
        idMatrices[rider][env] = idMat
        covarMatrices[rider][env] = (rollCovar, steerCovar)
        print('Done.')

# for each environment calculate the model for all riders
for env in environments:
    print('Computing the estimate for all riders on {}.'.format(env))
    runs = csi.select_runs(riders, maneuvers, [env])
    runs = list(set(runs).difference(errors))

    means = csi.mean_canon(riders, canon, H)
    idMat, rollCovar, steerCovar = csi.enforce_symmetry(runs, trials,
            rollParams, steerParams, *means)
    idMatrices[env] = idMat
    covarMatrices[env] = (rollCovar, steerCovar)
    print('Done.')

# save all of the identified matrices to file
with open('idMatrices.p', 'w') as f:
    cPickle.dump(idMatrices, f)

with open('covarMatrices.p', 'w') as f:
    cPickle.dump(covarMatrices, f)
