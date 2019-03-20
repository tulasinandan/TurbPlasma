import pandas as pd
import numpy as np
import Voyager as vy
import DataAnalysisRoutines as dar
from fbdr import *
from ffd import *

start = 1977
end = 1982
variables = ['Epoch2', 'scDistance', 'B1', 'B2', 'B3']
timevar = 'Epoch2'
SIGMA = 3.0
length = 30000000
FILL = np.NaN
ws = 100
BADDATAFILELIST = '/data/DATA-CODES/processed_data/BadIntervals_02012019/BadIntervals.txt'
SAVEFILEPATH = '/data/DATA-CODES/processed_data/Data_01082019/'
MAXVAL = 25
MINVAL = -25

# Create list of files to convert from CDF to dataframe
flst = vy.create_voyager_file_list_years(start,end)
d = dar.convert_cdfs_to_dataframe(flst, varlist=variables, time_var_str=timevar)
b1 = pd.Series(d['B1']).replace(to_replace = 999., value = FILL)
b2 = pd.Series(d['B2']).replace(to_replace = 999., value = FILL)
b3 = pd.Series(d['B3']).replace(to_replace = 999., value = FILL)
time = pd.Series(d['time'])
dist = pd.Series(d['scDistance']).replace(to_replace = 999., value = FILL)
dist[dist < 1] = np.NaN
timeidx = pd.date_range('1/1/2000',periods = len(dist), freq='48S')
newdist = pd.Series(d['scDistance'], index=timeidx)
newdist = newdist.resample('1.92S').asfreq().interpolate('linear',limit_direction='forward')
append = np.zeros(len(time)-len(newdist))
appended = pd.Series(np.append(newdist,append)).replace(to_replace=0, value=FILL)
appended.interpolate('linear', limit_direction = 'both')
data = pd.DataFrame(np.asarray([time, b1, b2, b3, appended]).transpose(), index=d['datetime'])
b1 = None
b2 = None
b3 = None
time = None
dist = None
d = None

# Making sure that there are no indices that were in the middle of files that went back in time
# If so, quickly remove them, + or minus 5, and go from there
t = data[0].reset_index(drop=True)
dt = t.shift(-1)-t
index = dt[dt<0].index
interestingtimes = np.asarray([]).astype(int)
for j in index:
        interestingtimes = np.append(interestingtimes,np.arange(j-5,j+5))
data.reset_index(drop=True, inplace=True)
data.drop(index=interestingtimes,inplace=True)
print('done with ridding data of negative dt')

# Fill the data gaps of the desired variables, original time values required
timefilled, distfilled = filldata(np.asarray(data[0]), np.asarray(data[4]), length, FILL)
timefilled, b1filled = filldata(np.asarray(data[0]), np.asarray(data[1]), length, FILL)
timefilled, b2filled = filldata(np.asarray(data[0]), np.asarray(data[2]), length, FILL)
timefilled, b3filled = filldata(np.asarray(data[0]), np.asarray(data[3]), length, FILL)
timefilled = np.trim_zeros(timefilled,'b')
distfilled = distfilled[:len(timefilled)]
b1filled = b1filled[:len(timefilled)]
b2filled = b2filled[:len(timefilled)]
b3filled = b3filled[:len(timefilled)]
print('Done Filling Time Gaps')

# Remove Data from provided file list
if BADDATAFILELIST != '':
	bd = np.loadtxt(BADDATAFILELIST).astype(int)
	for idx in bd:
		b1filled[idx[0]:idx[1]] = FILL 
		b2filled[idx[0]:idx[1]] = FILL
		b3filled[idx[0]:idx[1]] = FILL

### PYTHON VERSION, WILL RUN LONGER ###
for ws in [3,10,25]:
	b1filled = dar.drop_outliers(b1filled.copy(), MINVAL, MAXVAL, ws, SIGMA)
	b2filled = dar.drop_outliers(b2filled.copy(), MINVAL, MAXVAL, ws, SIGMA)
	b3filled = dar.drop_outliers(b3filled.copy(), MINVAL, MAXVAL, ws, SIGMA)
print('Done Cleaning')

# Convert from arrays to Series and Convert back to NaN
timefilled = pd.Series(np.asarray(timefilled));
distfilled = pd.Series(np.asarray(distfilled));
b1filled = pd.Series(b1filled); b1filled.replace(to_replace = FILL, value = np.NaN, inplace = True) 
b2filled = pd.Series(b2filled); b2filled.replace(to_replace = FILL, value = np.NaN, inplace = True) 
b3filled = pd.Series(b3filled); b3filled.replace(to_replace = FILL, value = np.NaN, inplace = True) 

# Save Cleaned and Filled data to filepath, The added part should change every run to avoid unwanted overwrites
timefilled.to_pickle(SAVEFILEPATH+'Vy1TimeFilledCleaned_'+str(start)+'_'+str(end)+'_03012019.pkl')
distfilled.to_pickle(SAVEFILEPATH+'Vy1DistanceFilledCleaned_'+str(start)+'_'+str(end)+'_03012019.pkl')
b1filled.to_pickle(SAVEFILEPATH+'Vy1BrFilledCleaned_'+str(start)+'_'+str(end)+'_03012019.pkl')
b2filled.to_pickle(SAVEFILEPATH+'Vy1BtFilledCleaned_'+str(start)+'_'+str(end)+'_03012019.pkl')
b3filled.to_pickle(SAVEFILEPATH+'Vy1BnFilledCleaned_'+str(start)+'_'+str(end)+'_03012019.pkl')
