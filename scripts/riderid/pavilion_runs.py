import numpy as np
import matplotlib.pyplot as plt
import bicycledataprocessor as bdp

data = bdp.DataSet()

# open the database for reading
data.open()

runTable = data.database.root.runTable

def query(row):

    man = row['Maneuver'].endswith('Disturbance')
    env = row['Environment'] == 'Pavillion Floor'
    cor = row['corrupt'] is not True and row['RunID'] not in range(326, 469)
    bik = row['Bicycle'] == 'Rigid Rider'

    if man and env and cor and bik:
        return True
    else:
        return False

good = [row['RunID'] for row in runTable if query(row)]

runs = []
for runNum in good:
    try:
        runs.append(bdp.Run(runNum, data.database, filterFreq=15.))
    except:
        pass

meanSpeeds = np.array([np.mean(r.taskSignals['ForwardSpeed']) for r in runs])

runIDs = np.array([int(r.metadata['RunID']) for r in runs])

riders = [r.metadata['Rider'] for r in runs]

lowSpeeds = meanSpeeds[meanSpeeds < 3.0]
medSpeeds = meanSpeeds[(3.0 < meanSpeeds) & (meanSpeeds < 4.5)]
highSpeeds = meanSpeeds[meanSpeeds > 4.5]

plt.boxplot(lowSpeeds, medSpeeds, highSpeeds)

# output mat files
for run in runs:
    run.export('mat', directory='../../data/pavilion')
