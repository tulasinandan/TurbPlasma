import numpy as np
import pandas as pd
import DataAnalysisRoutines as dar

def create_voyager_file_list_years(start,end):
    import glob
# set directory for voyager files
    basedir='/archive/cdaweb.gsfc.nasa.gov/pub/data/voyager/voyager1/'
    filelist=[]
# sort all the files collected in the for each year in the range and append to file list
    for i in range(start,end+1):
        l=sorted(glob.iglob(basedir+'magnetic_fields_cdaweb/mag_2s/'+str(i)+'/*_v01.cdf'))
        for j in l:
            #print(j)
            filelist.append(j)
    return filelist

def create_voyager_file_list_for_velocity_years(start,end):
	import glob
	basedir = '/archive/cdaweb.gsfc.nasa.gov/pub/data/voyager/voyager1/coho1hr_magplasma/'
	filelist = []
	for i in range(start,end+1):
		l = sorted(glob.iglob(basedir+str(i)+'/*_v01.cdf'))
		for j in l:
			filelist.append(j)
	return filelist


def create_voyager_file_list_days(year, day):
    import glob
# set directory for voyager files
    basedir='/archive/cdaweb.gsfc.nasa.gov/pub/data/voyager/voyager1/'
    filelist=[]
# retrieve the file in directory by specified year, month and day (yyyy, mmdd)
    l=sorted(glob.iglob(basedir+'magnetic_fields_cdaweb/mag_2s/'+str(year) + '/voyager1_2s_mag_' + str(year) + day+'_v01.cdf'))
    for j in l:
        #print(j)
        filelist.append(j)
    return filelist

def voyager_kurtosis(series, lags, au, delta_t, scaling):
        """
        Routine to take the kurtosis of voyager data specifically of 1 component of the magnetic field
        Input:
              series: the series to calculate the kurtosis of
              lags: the array of lags with units of scaling
              au: the AU that the series is within the range of au to au+1
              delta_t: the time difference between each point in series
              scaling: either to be di or lambda to determine the lag as a constant multiple of di or a constant fraction of lambda
        Output:
              k: the Calculated Kurtosis of the series at the given lag with proper scaling
              lags: the lags in di or the proper di varied to keep the constant fraction of lambda
        """
        var = pd.read_pickle('/data/DATA-CODES/processed_data/Data/VY1DailyAverageVelocityDensityDistance.pkl')
        v = var.v
        n = var.n
        d = var.dist


        aus = np.loadtxt('/data/DATA-CODES/processed_data/Data/NewAuIndices_06222018.txt').astype(int)
        start = d[d > au].reset_index()['index'][0]
        end = d[d > au+1].reset_index()['index'][0]

        v_arr = v[start:end]
        n_arr = n[start:end]

        if v_arr.mean() == np.NaN:
                v_arr_mean = 450
        else:
                v_arr_mean = v_arr.mean()
        if n_arr.mean() == np.NaN:
                n_arr_mean = var.n[var.dist < 2].mean()/(au**2)
        else:
                n_arr_mean = n_arr.mean()
        series = series.copy()
        lags = np.asarray(lags)
        if scaling == 'di':
                ptlag = dar.diToDeltaT(lags, v_arr_mean, n_arr_mean)/delta_t
        elif scaling == 'lambda':
                au1 = d[d > 1].reset_index()['index'][0]
                au2 = d[d > 2].reset_index()['index'][0]
                ptlag = dar.lambdaToDeltaT(lags, v_arr_mean, n[au1:au2].mean(), au)/delta_t
        k, ptlag = dar.kurtosis(series, ptlag)
        return k, lags

