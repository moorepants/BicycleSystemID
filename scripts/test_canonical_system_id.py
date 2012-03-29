import numpy as np
from numpy import testing
import bicycledataprocessor as bdp
import canonical_system_id as csi

def test_benchmark_time_series():
    dataset = bdp.DataSet()
    trial = bdp.Run('700', dataset)

    timeSeries = csi.benchmark_time_series(trial)

    testing.assert_allclose(trial.taskSignals['RollAngle'], timeSeries['p'])
    testing.assert_allclose(trial.taskSignals['RollRate'], timeSeries['pD'])
    testing.assert_allclose(trial.taskSignals['RollRate'].time_derivative(),
        timeSeries['pDD'])
    testing.assert_allclose(trial.taskSignals['SteerAngle'], timeSeries['d'])
    testing.assert_allclose(trial.taskSignals['SteerRate'], timeSeries['dD'])
    testing.assert_allclose(trial.taskSignals['SteerRate'].time_derivative(),
        timeSeries['dDD'])
    testing.assert_allclose(trial.taskSignals['ForwardSpeed'], timeSeries['v'])
    testing.assert_allclose(trial.bicycleRiderParameters['g'], timeSeries['g'])
    testing.assert_allclose(trial.taskSignals['PullForce'], timeSeries['F'])

def test_benchmark_lstsq_matrices():
    dataset = bdp.DataSet()
    trial = bdp.Run('700', dataset)

    A, B, F = csi.whipple_state_space(trial.metadata['Rider'], 1.0)
    H = np.dot(np.linalg.inv(B[2:]), F[2:])

    timeSeries = csi.benchmark_time_series(trial, subtractMean=False)
    M, C1, K0, K2 = trial.bicycle.canonical(nominal=True)
    fixedValues = csi.benchmark_canon_to_dict(M, C1, K0, K2, H)

    rollParams = ['Mpp', 'Mpd',
                  'C1pp', 'C1pd',
                  'K0pp', 'K0pd',
                  'K2pp', 'K2pd',
                  'HpF']

    A, B = csi.benchmark_lstsq_matrices(rollParams, timeSeries, fixedValues)

    testing.assert_allclose(trial.taskSignals['RollRate'].time_derivative(),
            A[:, 0])
    testing.assert_allclose(trial.taskSignals['SteerRate'].time_derivative(),
            A[:, 1])
    testing.assert_allclose(trial.taskSignals['ForwardSpeed'] *
            trial.taskSignals['RollRate'], A[:, 2])
    testing.assert_allclose(trial.taskSignals['ForwardSpeed'] *
            trial.taskSignals['SteerRate'], A[:, 3])
    testing.assert_allclose(9.81 *
            trial.taskSignals['RollAngle'], A[:, 4])
    testing.assert_allclose(9.81 *
            trial.taskSignals['SteerAngle'], A[:, 5])
    testing.assert_allclose(trial.taskSignals['ForwardSpeed']**2 *
            trial.taskSignals['RollAngle'], A[:, 6])
    testing.assert_allclose(trial.taskSignals['ForwardSpeed']**2 *
            trial.taskSignals['SteerAngle'], A[:, 7])
    testing.assert_allclose(-trial.taskSignals['PullForce'], A[:, 8])
    testing.assert_allclose(np.zeros_like(trial.taskSignals['PullForce']), B)

    rollParams = ['Mpp', 'Mpd',
                  'C1pp', 'C1pd',
                  'K0pp', 'K0pd',
                  'K2pp', 'K2pd']

    A, B = csi.benchmark_lstsq_matrices(rollParams, timeSeries, fixedValues)

    testing.assert_allclose(trial.taskSignals['RollRate'].time_derivative(),
            A[:, 0])
    testing.assert_allclose(trial.taskSignals['SteerRate'].time_derivative(),
            A[:, 1])
    testing.assert_allclose(trial.taskSignals['ForwardSpeed'] *
            trial.taskSignals['RollRate'], A[:, 2])
    testing.assert_allclose(trial.taskSignals['ForwardSpeed'] *
            trial.taskSignals['SteerRate'], A[:, 3])
    testing.assert_allclose(9.81 *
            trial.taskSignals['RollAngle'], A[:, 4])
    testing.assert_allclose(9.81 *
            trial.taskSignals['SteerAngle'], A[:, 5])
    testing.assert_allclose(trial.taskSignals['ForwardSpeed']**2 *
            trial.taskSignals['RollAngle'], A[:, 6])
    testing.assert_allclose(trial.taskSignals['ForwardSpeed']**2 *
            trial.taskSignals['SteerAngle'], A[:, 7])
    testing.assert_allclose(H[0] * trial.taskSignals['PullForce'], B)

    steerParams = ['Mdp', 'Mdd', 'C1dp', 'C1dd', 'K0dp', 'K0dd', 'K2dp',
            'K2dd', 'HdF']

    A, B = csi.benchmark_lstsq_matrices(steerParams, timeSeries, fixedValues)

    testing.assert_allclose(trial.taskSignals['RollRate'].time_derivative(),
            A[:, 0])
    testing.assert_allclose(trial.taskSignals['SteerRate'].time_derivative(),
            A[:, 1])

    testing.assert_allclose(trial.taskSignals['ForwardSpeed'] *
            trial.taskSignals['RollRate'], A[:, 2])
    testing.assert_allclose(trial.taskSignals['ForwardSpeed'] *
            trial.taskSignals['SteerRate'], A[:, 3])

    testing.assert_allclose(9.81 *
            trial.taskSignals['RollAngle'], A[:, 4])
    testing.assert_allclose(9.81 *
            trial.taskSignals['SteerAngle'], A[:, 5])

    testing.assert_allclose(trial.taskSignals['ForwardSpeed']**2 *
            trial.taskSignals['RollAngle'], A[:, 6])
    testing.assert_allclose(trial.taskSignals['ForwardSpeed']**2 *
            trial.taskSignals['SteerAngle'], A[:, 7])

    testing.assert_allclose(-trial.taskSignals['PullForce'], A[:, 8])
    testing.assert_allclose(trial.taskSignals['SteerTorque'], B)

    steerParams = ['Mdp', 'Mdd', 'C1dp', 'C1dd', 'K0dp', 'K0dd', 'K2dp',
            'K2dd']

    A, B = csi.benchmark_lstsq_matrices(steerParams, timeSeries, fixedValues)

    testing.assert_allclose(trial.taskSignals['RollRate'].time_derivative(),
            A[:, 0])
    testing.assert_allclose(trial.taskSignals['SteerRate'].time_derivative(),
            A[:, 1])

    testing.assert_allclose(trial.taskSignals['ForwardSpeed'] *
            trial.taskSignals['RollRate'], A[:, 2])
    testing.assert_allclose(trial.taskSignals['ForwardSpeed'] *
            trial.taskSignals['SteerRate'], A[:, 3])

    testing.assert_allclose(9.81 *
            trial.taskSignals['RollAngle'], A[:, 4])
    testing.assert_allclose(9.81 *
            trial.taskSignals['SteerAngle'], A[:, 5])

    testing.assert_allclose(trial.taskSignals['ForwardSpeed']**2 *
            trial.taskSignals['RollAngle'], A[:, 6])
    testing.assert_allclose(trial.taskSignals['ForwardSpeed']**2 *
            trial.taskSignals['SteerAngle'], A[:, 7])

    testing.assert_allclose(trial.taskSignals['SteerTorque'] +
            H[1] * trial.taskSignals['PullForce'], B)

    steerParams = ['Mdp', 'C1dd', 'K0dp', 'K2dp', 'K2dd']

    A, B = csi.benchmark_lstsq_matrices(steerParams, timeSeries, fixedValues)

    testing.assert_allclose(trial.taskSignals['RollRate'].time_derivative(),
            A[:, 0])
    testing.assert_allclose(trial.taskSignals['ForwardSpeed'] *
            trial.taskSignals['SteerRate'], A[:, 1])
    testing.assert_allclose(9.81 *
            trial.taskSignals['RollAngle'], A[:, 2])
    testing.assert_allclose(trial.taskSignals['ForwardSpeed']**2 *
            trial.taskSignals['RollAngle'], A[:, 3])
    testing.assert_allclose(trial.taskSignals['ForwardSpeed']**2 *
            trial.taskSignals['SteerAngle'], A[:, 4])

    testing.assert_allclose(
            H[1] * trial.taskSignals['PullForce'] +
            trial.taskSignals['SteerTorque'] -
            fixedValues['Mdd'] * trial.taskSignals['SteerRate'].time_derivative() -
            fixedValues['C1dp'] * trial.taskSignals['ForwardSpeed'] * trial.taskSignals['RollRate'] -
            fixedValues['K0dd'] * 9.81 * trial.taskSignals['SteerAngle'], B)

def test_benchmark_canon_to_dict():
    M = np.random.rand(2, 2)
    C1 = np.random.rand(2, 2)
    K0 = np.random.rand(2, 2)
    K2 = np.random.rand(2, 2)
    H = np.random.rand(2, 1)

    cor = {'Mpp': M[0, 0],
           'Mpd': M[0, 1],
           'C1pp': C1[0, 0],
           'C1pd': C1[0, 1],
           'K0pp': K0[0, 0],
           'K0pd': K0[0, 1],
           'K2pp': K2[0, 0],
           'K2pd': K2[0, 1],
           'Mdp': M[1, 0],
           'Mdd': M[1, 1],
           'C1dp': C1[1, 0],
           'C1dd': C1[1, 1],
           'K0dp': K0[1, 0],
           'K0dd': K0[1, 1],
           'K2dp': K2[1, 0],
           'K2dd': K2[1, 1],
           'HpF': H[0, 0],
           'HdF': H[1, 0]}

    entries = csi.benchmark_canon_to_dict(M, C1, K0, K2, H)

    for k, v in cor.items():
        assert v == entries[k]
