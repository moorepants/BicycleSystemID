import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
from bicycledataprocessor import main as bdp

dataset = bdp.DataSet(fileName='../BicycleDataProcessor/InstrumentedBicycleData.h5')
dataset.open(mode='r')

pathToPar = '../BicycleParameters/data'

unfilteredRun = bdp.Run('00644', dataset.database, pathToPar)

dataset.close()

freqFigs = {}
for name, signal in unfilteredRun.taskSignals.items():
    f, a = signal.frequency()
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(f, a)
    ax.set_title(name)
    ax.set_xlabel('Frequency [hz]')
    ax.set_ylabel('Amplitude [' + signal.units + ']')
    ax.set_xlim((0.0, 50.0))
    freqFigs[name] = fig

# plot the frequency spectrum of the outputs

#filteredRun = bdp.Run('00644', dataset.database, pathToPar, filterFreq=30.0)
