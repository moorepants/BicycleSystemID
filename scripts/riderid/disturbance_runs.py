#!/usr/bin/env python

# This script builds a directory of processed runs intended for identification
# of the rider controller. Only disturbance runs are included.

import bicycledataprocessor as bdp

dataset = bdp.DataSet()
dataset.open()

table = dataset.database.root.runTable

runs = []
for row in table.iterrows():
    con = []
    con.append(row['Rider'] in ['Jason', 'Charlie', 'Luke'])
    con.append(row['Maneuver'] == 'Balance With Disturbance' or
               row['Maneuver'] == 'Track Straight Line With Disturbance')
    con.append(row['corrupt'] is not True)
    con.append(int(row['RunID']) > 100)
    if False not in con:
        runs.append(row['RunID'])

dataset.close()

timeShiftErrors = []
otherErrors = []
for r in runs:
    try:
        trial = bdp.Run(r, dataset, filterFreq=15.)
    except bdp.bdpexceptions.TimeShiftError:
        timeShiftErrors.append(r)
        dataset.close()
    except IndexError:
        otherErrors.append(r)
        dataset.close()
    else:
        trial.export('mat', directory='../../data/riderid/disturbance-runs')

print timeShiftErrors
print otherErrors
