import numpy as np
import bicycleparameters as bp
import bicycledataprocessor as bdp
from dtk.bicycle import benchmark_to_moore, pitch_from_roll_and_steer
from LateralForce import LinearLateralForce

def select_runs(riders, maneuvers, environments):

    dataset = bdp.DataSet()
    dataset.open()

    table = dataset.database.root.runTable

    runs = []
    for row in table.iterrows():
        con = []
        con.append(row['Rider'] in riders)
        con.append(row['Maneuver'] in maneuvers)
        con.append(row['Environment'] in environments)
        con.append(row['corrupt'] is not True)
        con.append(int(row['RunID']) > 100)
        if False not in con:
            runs.append(row['RunID'])

    dataset.close()

    return runs

def benchmark_identified_matrices(freeParams, solutions, M, C1, K0, K2, H):
    """Returns the benchmark canonical matrices given the results from the
    linear regression.

    Parameters
    ----------
    freeParams : tuple
        The names of the parameters that were solved for, (roll, steer).
    solutions : tuple
        The values of the parameters that were solved for, (roll, steer).
    canon : tuple
        The M, C1, K0, K2, H matrices computed from first principles.

    Returns
    -------
    M : ndarray, shape(2,2)
    C1 : ndarray, shape(2,2)
    K0 : ndarray, shape(2,2)
    K2 : ndarray, shape(2,2)
    H : ndarray, shape(2,1)

    """

    rollFree, steerFree = freeParams
    free = {'p': rollFree, 'd': steerFree}

    rollSol, steerSol = solutions
    sol = {'p': rollSol, 'd': steerSol}

    mat = {'M': M.copy(), 'C1': C1.copy(), 'K0': K0.copy(), 'K2': K2.copy(),
            'H': H.copy().reshape(2, 1)}

    cols = ['p', 'd']
    for eq, params in free.items():
        for par in params:
            if par.startswith('M') or par.startswith('H'):
                matName = par[0]
            else:
                matName = par[:2]

            i = cols.index(eq)
            if par[-1] == 'F':
                j = 0
            else:
                j = cols.index(par[-1])

            mat[matName][i, j] = sol[eq][params.index(par)]

    return mat['M'], mat['C1'], mat['K0'], mat['K2'], mat['H']

def benchmark_canon_to_dict(M, C1, K0, K2, H):
    """Returns the entries of the benchmark model in a keyword dictionary form.

    Parameters
    ----------
    M : array_like, shape(2,2)
        The mass matrix.
    C1 : array_like, shape(2,2)
        The damping matrix which is proportional to speed.
    K0 : array_like, shape(2,2)
        The stiffness matrix which is proportional to the acceleration due to
        gravity.
    K2 : array_like, shape(2,2)
        The stiffness matrix which is proportional to the square of the speed.
    H : array_like, shape(2)
        The lateral force proportions.

    Returns
    -------
    entries : dictionary
        All of the values in the input arrays labeled

    Notes
    -----

    M * q'' + v * C1 * q' + [g * K0 + v^2 * K2] * q = T + H * F

    States:
    q = [phi,
         delta]

    Input torques:
    T = [Tphi,
         Tdelta]

    Later force:
    F = Fcl

    """
    rollParams = ['Mpp', 'Mpd',
                  'C1pp', 'C1pd',
                  'K0pp', 'K0pd',
                  'K2pp', 'K2pd',
                  'HpF']
    steerParams = ['Mdp', 'Mdd',
                   'C1dp', 'C1dd',
                   'K0dp', 'K0dd',
                   'K2dp', 'K2dd',
                   'HdF']

    params = [rollParams, steerParams]

    total = np.hstack((M, C1, K0, K2, H.reshape(2, 1)))

    entries = {}
    for i, row in enumerate(total):
        for j, val in enumerate(row):
            entries[params[i][j]] = val

    return entries

def benchmark_lstsq_matrices(freeParams, timeSeries, fixedValues):
    """Returns the matrices formed for least squares computation.

    freeParams : list
        A list of the free parameters for either the roll equation or the steer
        equation. These can be a subset but there needs to be at least one item
        in the list.

        Roll : ['Mpp', 'Mpd', 'C1pp', 'C1pd', 'K0pp', 'K0pd', 'K2pp', 'K2pd',
        'HpF']
        Steer : ['Mdp', 'Mdd', 'C1dp', 'C1dd', 'K0dp', 'K0dd', 'K2dp', 'K2dd',
        'HdF']

    timeSeries : dictionary
        Contains the coordinates, rates, accels of phi and delta, forces T
        and F, and forward speed v.

    fixedValues : dictionary
        The default parameters for the entries.

    """
    if len(freeParams) == 0:
        raise ValueError('You must provide at least one free parameter.')

    rollParams = ['Mpp', 'Mpd',
                  'C1pp', 'C1pd',
                  'K0pp', 'K0pd',
                  'K2pp', 'K2pd',
                  'HpF']

    steerParams = ['Mdp', 'Mdd',
                   'C1dp', 'C1dd',
                   'K0dp', 'K0dd',
                   'K2dp', 'K2dd',
                   'HdF']

    if set(freeParams).issubset(rollParams):
        params = rollParams
        B = np.zeros_like(timeSeries['F'])
    elif set(freeParams).issubset(steerParams):
        params = steerParams
        B = timeSeries['T'].copy()
    else:
        raise ValueError('The free parameters must be a subset of the roll or\
                steer parameters.')

    signals = ['pDD', 'dDD',
               'v*pD', 'v*dD',
               'g*p', 'g*d',
               'v2*p', 'v2*d',
               'F']

    A = np.zeros((len(timeSeries['p']), len(freeParams)))

    def split_sig(i, signals, timeSeries):
        if signals[i].startswith('v2*'):
            v, s = signals[i].split('2*')
            sig = timeSeries[v]**2 * timeSeries[s]

        elif signals[i].startswith('v*'):
            v, s = signals[i].split('*')
            sig = timeSeries[v] * timeSeries[s]

        elif signals[i].startswith('g*'):
            g, s = signals[i].split('*')
            sig = timeSeries[g] * timeSeries[s]

        elif signals[i] == 'F':
            sig = -timeSeries[signals[i]]

        else:
            sig = timeSeries[signals[i]]

        return sig

    for i, par in enumerate(params):
        if par in freeParams:
            A[:, freeParams.index(par)] = split_sig(i, signals, timeSeries)
        else:
            B -= fixedValues[par] * split_sig(i, signals, timeSeries)

    return A, B

def benchmark_time_series(trial, subtractMean=False):
    """Returns a dictionary of signals for the benchmark canonical equations.

    Parameters
    ----------
    trial : Run
        A run object with computed task signals.
    subtracMean : boolean, optional.
        If true the mean will be subtracted from all signals except the lateral
        force and forward speed.

    Returns
    -------
    timeSeries : dictionary
        Contains all of the signals needed to compute:

        M * [pDD, + v * C1 * [pD, + [g * K0 + v^2 * K2] * [p, = [0, + H * F
             dDD]             dD]                          d]    T]

    """
    timeSeries = {'p' : trial.taskSignals['RollAngle'],
                  'd' : trial.taskSignals['SteerAngle'],
                  'pD' : trial.taskSignals['RollRate'],
                  'dD' : trial.taskSignals['SteerRate'],
                  'pDD' : trial.taskSignals['RollRate'].time_derivative(),
                  'dDD' : trial.taskSignals['SteerRate'].time_derivative(),
                  'T' : trial.taskSignals['SteerTorque'],
                  'v' : trial.taskSignals['ForwardSpeed'],
                  'g' : trial.bicycleRiderParameters['g'],
                  }

    if trial.metadata['Maneuver'].endswith('Disturbance'):
        timeSeries['F'] = trial.taskSignals['PullForce']
    else:
        timeSeries['F'] = np.zeros_like(trial.taskSignals['PullForce'])

    # subtract the mean of the angles, rates, accelerations, and steer torque
    if subtractMean is True:
        for k, v in timeSeries.items():
            if k.startswith('p') or k.startswith('d') or k.startswith('T'):
                timeSeries[k] = v.subtract_mean()

    # disconnect the arrays from the trials
    for k, v in timeSeries.items():
        try:
            timeSeries[k] = v.copy()
        except AttributeError:
            pass

    return timeSeries

def whipple_state_space(rider, speed):
    """ x' = Ax + Bu + Fv

    x = [phi,
         delta,
         phiDot,
         deltaDot]
    u = [Tphi,
         Tdel]
    v = [Fcl]

    """

    bicycleModel = LinearLateralForce()

    pathToData='/media/Data/Documents/School/UC Davis/Bicycle Mechanics/BicycleParameters/data/'
    if rider == 'Jason':
        bicycleName = 'Rigid'
    elif rider == 'Charlie' or rider == 'Luke':
        bicycleName = 'Rigidcl'
    bicycle = bp.Bicycle(bicycleName, pathToData)
    bicycle.add_rider(rider)

    # set the model parameters
    benchmarkPar = bp.io.remove_uncertainties(bicycle.parameters['Benchmark'])
    benchmarkPar['xcl'] = bicycle.parameters['Measured']['xcl']
    benchmarkPar['zcl'] = bicycle.parameters['Measured']['zcl']
    moorePar = benchmark_to_moore(benchmarkPar)
    bicycleModel.set_parameters(moorePar)

    # set the default equilibrium point
    pitchAngle = pitch_from_roll_and_steer(0., 0., moorePar['rf'],
            moorePar['rr'], moorePar['d1'], moorePar['d2'], moorePar['d3'])
    wheelAngSpeed = -speed / moorePar['rr']
    equilibrium = np.zeros(len(bicycleModel.stateNames))
    equilibrium[bicycleModel.stateNames.index('q5')] = pitchAngle
    equilibrium[bicycleModel.stateNames.index('u6')] = wheelAngSpeed
    bicycleModel.linear(equilibrium)

    states = ['q4', 'q7', 'u4', 'u7']
    inputs = ['Fcl', 'T4', 'T7']
    outputs = ['q4', 'q7', 'u4', 'u7']

    A, B, C, D = bicycleModel.reduce_system(states, inputs, outputs)

    F = B[:, 0]
    B = B[:, 1:]

    return A, B, F

def time_series(trial, F):
    coordinates = np.vstack((trial.taskSignals['RollAngle'],
                        trial.taskSignals['SteerAngle']))
    rates = np.vstack((trial.taskSignals['RollRate'],
                        trial.taskSignals['SteerRate']))
    accels = np.vstack((trial.taskSignals['RollRate'].time_derivative(),
                        trial.taskSignals['SteerRate'].time_derivative()))
    if trial.metadata['Maneuver'].endswith('Disturbance'):
        latForce = trial.taskSignals['PullForce']
    else:
        latForce = np.zeros_like(trial.taskSignals['PullForce'])

    forces = np.vstack((F[trial.metadata['Rider']][0] * latForce,
        F[trial.metadata['Rider']][1] * latForce +
        trial.taskSignals['SteerTorque']))

    return coordinates, rates, accels, forces

def compute_unknowns(thetaPhi, thetaDelta, timeSeries, canon):

    solPhi = compute_coefficients(thetaPhi, timeSeries, canon)
    solDelta = compute_coefficients(thetaDelta, timeSeries, canon)
    M_id, C_id, K_id = identified_matrices([thetaPhi, thetaDelta], [solPhi,
        solDelta], canon)

    return M_id, C_id, K_id

def compute_coefficients(theta, timeSeries, canon):
    coordinates, rates, accels, forces = timeSeries
    M, C, K = canon
    massMat = [['m11', 'm12'],
               ['m21', 'm22']]
    dampMat = [['c11', 'c12'],
               ['c21', 'c22']]
    stifMat = [['k11', 'k12'],
               ['k21', 'k22']]

    coefs = [massMat[0] + dampMat[0] + stifMat[0],
             massMat[1] + dampMat[1] + stifMat[1]]

    if theta[0] in coefs[0]:
        eqIndex = 0
    elif theta[0] in coefs[1]:
        eqIndex = 1

    knownCoef = list(set(coefs[eqIndex]).difference(theta))

    A = np.zeros((coordinates.shape[1], len(theta)))
    B = forces[eqIndex]

    for i, par in enumerate(theta):
        if par in massMat[eqIndex]:
            col = accels[massMat[eqIndex].index(par)]
        elif par in dampMat[eqIndex]:
            col = rates[dampMat[eqIndex].index(par)]
        elif par in stifMat[eqIndex]:
            col = coordinates[stifMat[eqIndex].index(par)]
        A[:, i] = col

    for coef in knownCoef:
        if coef in massMat[eqIndex]:
            a = accels[massMat[eqIndex].index(coef)]
            B = B - M[eqIndex, massMat[eqIndex].index(coef)] * a
        elif coef in dampMat[eqIndex]:
            r = accels[dampMat[eqIndex].index(coef)]
            B = B - C[eqIndex, dampMat[eqIndex].index(coef)] * r
        elif coef in stifMat[eqIndex]:
            c = accels[stifMat[eqIndex].index(coef)]
            B = B - K[eqIndex, stifMat[eqIndex].index(coef)] * c

    x = np.linalg.lstsq(A, B)[0]

    coef = {k : v for k, v in zip(theta, x)}

    return coef

def identified_matrices(theta, sol, canon):
    M, C, K = canon
    massMat = [['m11', 'm12'],
               ['m21', 'm22']]
    dampMat = [['c11', 'c12'],
               ['c21', 'c22']]
    stifMat = [['k11', 'k12'],
               ['k21', 'k22']]

    coefs = [massMat[0] + dampMat[0] + stifMat[0],
             massMat[1] + dampMat[1] + stifMat[1]]

    massMat = np.array(massMat)
    dampMat = np.array(dampMat)
    stifMat = np.array(stifMat)

    M_id = np.zeros(massMat.shape)
    C_id = np.zeros(dampMat.shape)
    K_id = np.zeros(stifMat.shape)

    for i, row in enumerate(coefs):
        for j, par in enumerate(row):

            if par.startswith('m'):
                indices = np.nonzero(massMat == par)
                if par in theta[i]:
                    M_id[indices] = sol[i][par]
                else:
                    M_id[indices] = M[indices]

            elif par.startswith('c'):
                indices = np.nonzero(dampMat == par)
                if par in theta[i]:
                    C_id[indices] = sol[i][par]
                else:
                    C_id[indices] = C[indices]

            elif par.startswith('k'):
                indices = np.nonzero(stifMat == par)
                if par in theta[i]:
                    K_id[indices] = sol[i][par]
                else:
                    K_id[indices] = M[indices]

    return M_id, C_id, K_id
