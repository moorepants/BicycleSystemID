#!/usr/bin/env python

# This script builds a directory of processed runs intended for identification
# of the bicycle model.

import bicycledataprocessor as bdp

dataset = bdp.DataSet()
dataset.open()

table = dataset.database.root.runTable

runs = []
for row in table.iterrows():
    con = []
    con.append(row['Rider'] in ['Jason', 'Charlie', 'Luke'])
    con.append(row['Maneuver'] == 'Balance' or
               row['Maneuver'] == 'Track Straight Line' or
               row['Maneuver'] == 'Balance With Disturbance' or
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
        trial.verify_time_sync(show=False, saveDir='time-sync-plots')
    except bdp.bdpexceptions.TimeShiftError:
        timeShiftErrors.append(r)
        dataset.close()
    except:
        otherErrors.append(r)
        dataset.close()
    else:
        trial.export('mat')

f = open('time-shift-errors.txt', 'w')
for error in timeShiftErrors:
    f.write(str(error) + '\n')
f.close()

f = open('other-errors.txt', 'w')
for error in otherErrors:
    f.write(str(error) + '\n')
f.close()
