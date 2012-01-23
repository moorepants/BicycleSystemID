import numpy as np
import matplotlib.pyplot as plt

from bicycledataprocessor import main

# the path to the potentially necessary files
pathToBicycle = '/media/Data/Documents/School/UC Davis/Bicycle Mechanics/'
pathToDB = pathToBicycle + 'BicycleDataProcessor/InstrumentedBicycleData.h5'
pathToH5 = pathToBicycle + 'BicycleDAQ/data/h5'
pathToPar = pathToBicycle + 'BicycleParameters/data'
pathToCorrupt = pathToBicycle + 'BicycleDataProcessor/data-corruption.csv'

data = main.DataSet(fileName=pathToDB, pathToH5=pathToH5,
        pathToCorruption=pathToCorrupt)

# open the database for reading
data.open()

runTable = data.database.root.runTable

def query(row):

    man = row['Maneuver'] == 'Track Straight Line With Disturbance'
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
        runs.append(main.Run(runNum, data.database, pathToPar, filterFreq=30.))
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
    run.export('mat', directory='pavilion')
