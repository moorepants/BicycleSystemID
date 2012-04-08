import numpy as np
from scipy.io import loadmat
from uncertainties import ufloat
import bicycleparameters as bp
import bicycledataprocessor as bdp
from dtk.bicycle import benchmark_to_moore, pitch_from_roll_and_steer
from LateralForce import LinearLateralForce

def create_rst_table(tableData, roll, steer, fileName=None):
    """Returns a reStructuredText version of the table data.

    Parameters
    ----------
    tableData : list
        A list of rows for the table.
    roll : list
        The free parameters in the roll torque equation.
    steer : list
        The free parameters in the steer torque equation.
    fileName : string
        If a path to a file is given, the table will be written to that
        file.

    Returns
    -------
    rstTable : string
        reStructuredText version of the table.

    """

    latexMap = {'Mpp' : r'M_{\phi\phi}',
                'Mpd' : r'M_{\phi\delta}',
                'Mdp' : r'M_{\delta\phi}',
                'Mdd' : r'M_{\delta\delta}',
                'C1pp' : r'C_{1\phi\phi}',
                'C1pd' : r'C_{1\phi\delta}',
                'C1dp' : r'C_{1\delta\phi}',
                'C1dd' : r'C_{1\delta\delta}',
                'K0pp' : r'K_{0\phi\phi}',
                'K0pd' : r'K_{0\phi\delta}',
                'K0dp' : r'K_{0\delta\phi}',
                'K0dd' : r'K_{0\delta\delta}',
                'K2pp' : r'K_{2\phi\phi}',
                'K2pd' : r'K_{2\phi\delta}',
                'K2dp' : r'K_{2\delta\phi}',
                'K2dd' : r'K_{2\delta\delta}',
                'HpF' : r'H_{\phi F}',
                'HdF' : r'H_{\delta F}',
                }

    # top row
    head = [''] + [':math:`' + latexMap[p] + '`' for p in (roll + steer)]
    subhead = [['R', 'E'] + ['Value', ':math:`\sigma`', '% Difference'] * (len(roll) +
            len(steer))]

    allData = subhead + tableData

    # find the longest string in each column
    largest = [len(string) for string in tableData[0]]
    for row in allData:
        colSize = [len(string) for string in row]
        for i, pair in enumerate(zip(colSize, largest)):
            if pair[0] > pair[1]:
                largest[i] = pair[0]

    rstTable = ''

    rstTable += '+' + (sum(largest[:2]) + 5) * '=' + '+'

    for i in range(len(head[1:])):
        rstTable += (sum(largest[i * 3 + 2:i * 3 + 5]) + 8) * '=' + '+'
    rstTable += '\n'

    rstTable += '|' + (sum(largest[:2]) + 5) * ' ' + '|'
    for i, par in enumerate(head[1:]):
        x = sum(largest[i * 3 + 2:i * 3 + 5]) + 8 - len(par) - 1
        rstTable += ' ' + par + ' ' * x + '|'
    rstTable += '\n'

    for j, row in enumerate(allData):
        if j in [0, 1]:
            dash = '='
        else:
            dash = '-'

        line = ''
        for i in range(len(row)):
            line += '+' + dash * (largest[i] + 2)

        line += '+\n|'
        for i, item in enumerate(row):
            line += ' ' + item + ' ' * (largest[i] - len(item)) + ' |'
        line += '\n'
        rstTable += line

    for num in largest:
       rstTable += '+' + dash * (num + 2)
    rstTable += '+'

    if fileName is not None:
        f = open(fileName, 'w')
        f.write(rstTable)
        f.close()

    return rstTable

def table_data(roll, steer, mat, cov):
    """Returns a data set of string values for text output and formatting into
    a table.

    Parameters
    ----------
    roll : list
        The free parameters in the roll torque equation.
    steer : list
        The free parameters in the steer torque equation.
    mat : dictionary
        A dictionary of the resulting models.
    cov : dictionary
        A dictionary of the resulting covariance matrices of the free
        parameters.

    Returns
    -------
    data : list
        A list of rows for the table.

    """
    data = []
    allRiders = ['Charlie', 'Jason', 'Luke']
    # create the comparison models
    canon = load_benchmark_canon(allRiders)
    H = lateral_force_contribution(allRiders)

    def add_to_row(row, covar):

        def add_par_to_row(i, row, parameters):
            value = values[par]
            sigma = np.sqrt(covar[i].diagonal()[parameters.index(par)])
            uVal = ufloat((value, sigma))
            uStr = bp.tables.uround(uVal)
            valStr, sigStr = uStr.split('+/-')
            diff = (value - theory[par]) / theory[par]
            row += [valStr, sigStr, '{:.1%}'.format(diff)]

        for par in roll:
            add_par_to_row(0, row, roll)
        for par in steer:
            add_par_to_row(1, row, steer)

    for k, v in mat.items():
        if k == 'All':
            row = ['A', 'A']
            values = benchmark_canon_to_dict(*v)
            mean = mean_canon(allRiders, canon, H)
            theory = benchmark_canon_to_dict(*mean)
            covar = cov[k]
            add_to_row(row, covar)
            data.append(row)
        elif k == 'Horse Treadmill' or k == 'Pavillion Floor':
            row = ['A', k[0]]
            values = benchmark_canon_to_dict(*v)
            mean = mean_canon(allRiders, canon, H)
            theory = benchmark_canon_to_dict(*mean)
            covar = cov[k]
            add_to_row(row, covar)
            data.append(row)
        elif k in allRiders:
            M, C, K1, K2 = canon[k]
            theory = benchmark_canon_to_dict(M, C, K1, K2, H[k])
            for envName, envVal in v.items():
                if envName == 'All':
                    row = [k[0], 'A']
                    values = benchmark_canon_to_dict(*envVal)
                    covar = cov[k][envName]
                    add_to_row(row, covar)
                    data.append(row)
                else:
                    row = [k[0], envName[0]]
                    values = benchmark_canon_to_dict(*envVal)
                    covar = cov[k][envName]
                    add_to_row(row, covar)
                    data.append(row)
    return data

def benchmark_canonical_variance(A, B, xhat, resid):
    """Returns the variance in the ordinary least squares fit and the
    covariance matrix of the estimated parameters.

    Parameters
    ----------
    A : ndarray, shape(n,d)
        The left hand side matrix in Ax=B.
    B : ndarray, shape(n,)
        The right hand side vector in Ax=B.
    xhat : ndarray, shape(d)
        The best estimate of the free parameters.

    Returns
    -------
    var : float
        The variance of the fit.
    covar : ndarray, shape(d,d)
        The covariance of the parameters.

    """
    # I'm getting a memory error when trying to do this dot product even though
    # I have plenty of memory available.
    #Bhat = np.dot(A, xhat)
    # or this gives an error because I gave the wrong type for Bhat
    #Bhat = np.zeros_like(B)
    #np.dot(A, xhat, out=Bhat)
    #e = B - Bhat
    #del Bhat
    #var = np.dot(e.T, e) / (A.shape[0] - A.shape[1])

    # I am pretty sure that the residues from numpy.linalg.lstsq is the SSE
    # (the residual sum of squares). So my calculation above is not needed.

    var = resid / (A.shape[0] - A.shape[1])

    covar = var * np.linalg.inv(np.dot(A.T, A))

    return covar

def mean_arm(riders):
    """Returns a mean arm model for the given riders.

    Parameters
    ----------
    riders : list
        All or a subset of ['Charlie', 'Jason', 'Luke']

    Returns
    -------
    A : ndarray, shape(19, 19)
        The mean state matrix across the given riders.
    B : ndarray, shape(19, 4)
        The mean input matrix across the given riders.

    Notes
    -----
    This expects there to be a mat file in the current directory called
    `armsAB-<rider>.mat` for each rider (e.g. `armsAB-Charlie.mat`). This file
    should contain three things:

    speed : double, size(n)
        The constant speeds at which the state space was evaluated at.
    stateMatrices : double, size(n, 19, 19)
        The state matrices of the arm model for each speed.
    inputMatrices : double, size(n, 19, 4)
        The input matrices of the arm model for each speed.

    """

    data = {}
    for rider in riders:
        m = loadmat('armsAB-' + rider + '.mat', squeeze_me=True)
        data[rider] = {}
        data[rider]['speed'] = m['speed']
        data[rider]['A'] = m['stateMatrices']
        data[rider]['B'] = m['inputMatrices']

    x, y, z = data[riders[0]]['A'].shape
    As = np.zeros((len(riders), x, y, z))

    x, y, z = data[riders[0]]['B'].shape
    Bs = np.zeros((len(riders), x, y, z))

    for i, rider in enumerate(riders):
        As[i] = data[rider]['A']
        Bs[i] = data[rider]['B']

    return As.mean(axis=0), Bs.mean(axis=0), data[riders[0]]['speed']

def enforce_symmetry(runNums, trials, rollParams, steerParams, M, C1, K0, K2,
        H):
    """Returns the best estimate of the free parameters in the benchmark
    canonical matrices with respect to the given runs. The symmetry in the M
    and K0 matrices is enforce by finding the terms in the roll equation first
    and fixing the result in the steer equation.

    Parameters
    ----------
    runNums : list
        A list of run numbers to be used to formulate the least squares
        problem.
    trials : dictionary
        A dictionary of bicycledataprocess.Run objects. All in `runNums` must
        be present. These trials have an additional attribute `H` which holds
        the value for the bicycle model.
    rollParams : list
        The free parameters in the roll equation.
    steerParmams : list
        The free parameters in the steer equation.
    M, C1, K0, K2, H : ndarrays
        The the parameters not in `rollParams` and `steerParams` will be fixed
        to the values in these arrays.

    Returns
    -------
    idMatrices : tuple
        The identified model (M_id, C1_id, K0_id, K2_id, H_id) for the given
        set of runs.
    rollCoVar : ndarray
        The covariance matrix of the identified roll parameters.
    steerCoVar : ndarray
        The covariance matrix of the identified steer parameters.

    """

    rollAs, rollBs = lstsq_A_B(trials, rollParams)
    totRollA, totRollB = stack_A_B(rollAs, rollBs, runNums)
    rollSol, rollVar = np.linalg.lstsq(totRollA, totRollB)[:2]
    rollCoVar = benchmark_canonical_variance(totRollA, totRollB,
            rollSol, rollVar[0])

    overwrite = {'Mdp': rollSol[rollParams.index('Mpd')],
                 'K0dp': rollSol[rollParams.index('K0pd')]}
    steerAs, steerBs = lstsq_A_B(trials, steerParams, overwrite=overwrite)
    totSteerA, totSteerB = stack_A_B(steerAs, steerBs, runNums)
    steerSol, steerVar = np.linalg.lstsq(totSteerA, totSteerB)[:2]
    steerCoVar = benchmark_canonical_variance(totSteerA, totSteerB,
            steerSol, steerVar[0])

    M_mod = M.copy()
    M_mod[1, 0] = overwrite['Mdp']
    K0_mod = K0.copy()
    K0_mod[1, 0] = overwrite['K0dp']

    idMatrices = benchmark_identified_matrices((rollParams,
        steerParams), (rollSol, steerSol), M_mod, C1, K0_mod, K2, H)

    return idMatrices, rollCoVar, steerCoVar

def load_trials(runs, H):
    """Returns a dictionary of Run objects which were successfully computed and
    list of runs which had errors.

    Parameters
    ----------
    runs : list
        A list of run numbers to try to load.
    H : dictionary
        The H vector for the lateral force contributions to steer torque and
        roll torque for each rider.

    Returns
    -------
    trials : dictionary
        A dictionary of bicycledataprocess.Run objects. All in `runNums` must
        be present. These trials have an additional attribute `H` which holds
        the value for the bicycle model.
    errors : list
        A list of run numbers that had errors (primarily in the bump finding.


    """
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
            trial.H = H[trial.metadata['Rider']]
            trials[r] = trial

    return trials, errors

def mean_canon(riders, canon, H):
    """Returns the mean benchmark canoncial matrices across the given riders.

    Parameters
    ----------
    riders : list
        A list of riders to take the mean with respect to.
    canon : dictionary
        The M, C1, K0, and K2 matrices for each rider.
    H : dictionary
        The H vector for each rider.

    Returns
    -------
    M, C1, K0, K2, H : ndarrays
        The mean matrices across riders.

    """
    Ms = np.zeros((len(riders), 2, 2))
    C1s = np.zeros((len(riders), 2, 2))
    K0s = np.zeros((len(riders), 2, 2))
    K2s = np.zeros((len(riders), 2, 2))
    Hs = np.zeros((len(riders), 2, 1))

    for i, rider in enumerate(riders):
        M, C1, K0, K2 = canon[rider]
        h = H[rider]
        Ms[i] = M
        C1s[i] = C1
        K0s[i] = K0
        K2s[i] = K2
        Hs[i] = h.reshape(2, 1)

    return Ms.mean(0), C1s.mean(0), K0s.mean(0), K2s.mean(0), Hs.mean(0)

def load_benchmark_canon(riders):
    """Returns the benchmark canonical matrices for each rider.

    Parameters
    ----------
    riders : list
        A list of riders from ['Charlie', 'Jason', 'Luke'].

    Returns
    -------
    canon : dictionary
        The M, C1, K0, and K2 matrices for each rider on the instrumented
        bicycle.

    """
    path = '/media/Data/Documents/School/UC Davis/Bicycle Mechanics/BicycleParameters/data'
    canon = {}
    for rider in riders:
        if rider == 'Jason':
            bName = 'Rigid'
        else:
            bName= 'Rigidcl'
        bicycle = bp.Bicycle(bName, path)
        bicycle.add_rider(rider)
        canon[rider] = bicycle.canonical(nominal=True)

    return canon

def lateral_force_contribution(riders):
    """Returns the lateral force contribution vector.

    Parameters
    ----------
    riders : list
        A list of riders from ['Charlie', 'Jason', 'Luke'].

    Returns
    -------
    H : dictionary
        The lateral force contribution vector for the canonical form for each
        each rider.

    """
    H = {}
    for rider in riders:
        A, B, F = whipple_state_space(rider, 1.0)
        H[rider] = np.dot(np.linalg.inv(B[2:]), F[2:])

    return H

def lstsq_A_B(trials, params, overwrite=None):
    """Returns the known matrices for the least squares problems for multiple
    Runs.

    Parameters
    ----------
    trials : dictionary
        A dictionary of bicycledataprocessor.Run objects with the keyword being
        the run number.
    params : list
        This should be a list of either the free roll parameters or free steer
        parameters.
    overwrite : dictionary
        Specify this dictionary of matrix entries if you want to overwrite the
        any fixed values for the right hand side of the Ax=B.

    Returns
    -------
    eachA : dictionary
        The left hand side matrix, A, for each trial in Ax=B.
    eachB : dicationary
        The right hand side vector, B, for each trial in Ax=B.

    """

    eachA = {}
    eachB = {}

    for num, trial in trials.items():
        timeSeries = benchmark_time_series(trial, subtractMean=True)
        M, C1, K0, K2 = trial.bicycle.canonical(nominal=True)
        fixedValues = benchmark_canon_to_dict(M, C1, K0, K2, trial.H)

        if overwrite is not None:
            for k, v in overwrite.items():
                fixedValues[k] = v

        A, B = benchmark_lstsq_matrices(params, timeSeries, fixedValues)

        eachA[num] = A
        eachB[num] = B

    return eachA, eachB

def stack_A_B(As, Bs, runNums):
    """Returns the right and left hand side matrices for a set of runs.

    Parameters
    ----------
    As : dictionary
        A dictionary of the left hand side matrices for each run.
    Bs : dictionary
        A dictionary of the right hand side matrices for each run.
    runNums : list
        The run numbers to include in the matrices.

    Returns
    -------
    totA : ndarray
        The left hand side matrix for the set of runs.
    totB : ndarray
        The right hand side matric for the set of runs.

    """
    for num in runNums:
        A = As[num]
        B = Bs[num]

        try:
            totA = np.vstack((totA, A))
            totB = np.hstack((totB, B))
        except NameError:
            totA = A
            totB = B

    return totA, totB

def select_runs(riders, maneuvers, environments):
    """Returns a list of runs given a set of conditions.

    Parameters
    ----------
    riders : list
        All or a subset of ['Charlie', 'Jason', 'Luke'].
    maneuvers : list
        All or a subset of ['Balance', 'Balance With Disturbance', 'Track
        Straight Line', 'Track Straight Line With Disturbance'].
    environments : list
        All or a subset of ['Horse Treadmill', 'Pavillion Floor'].

    Returns
    -------
    runs : list
        List of run numbers for the given conditions.

    """

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

## the following deal with finding the entries to M, C and K

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
