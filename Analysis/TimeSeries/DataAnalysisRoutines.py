# -*- coding: utf-8 -*-
"""
Created on Thu Jun 30 20:08:29 2016

@authors: Manuel, Tulasi
"""

import pandas as pd
import numpy as np


def logint(maxval):
	"""
        Routine to return logarithmically spaced integers as:
        1 2 3 4 5 6 7 8 9 10 20 30 ... 100 200 ... 1000 ... 10000 ....
	Input:
		maxval: The max value that the array will not surpass
	Output:
		lst: The list to be converted to an numpy array that is
			in logarithmic space
	"""
	pwr=0; lst=[0]
	while True:
		n=0
		for i in range(1,10):
			num=(n+i)*10**pwr
			lst.append(num)
			if num >= maxval:
				break;
		pwr += 1
		if num >=maxval:
			break;
	
	return np.asarray(lst)

def gen_log_space(limit, n):
	"""
	Routine to create an array of visually equidistant points in log space
	Input:
		limit: the max value the function will create its list to
		n: the number of points in the log space
	Output:
		result: A log space that are visually equidistant in a log space.
	"""
	result = [2]
	if n>1:  # just a check to avoid ZeroDivisionError
		ratio = (float(limit)/result[-1]) ** (1.0/(n-len(result)))
	while len(result)<n:
		next_value = result[-1]*ratio
		if next_value - result[-1] >= 1:
			# safe zone. next_value will be a different integer
	        	result.append(next_value)
		else:
	        	# problem! same integer. we need to find next_value by artificially incrementing previous value
	        	result.append(result[-1]+1)
	        	# recalculate the ratio so that the remaining values will scale correctly
	        	ratio = (float(limit)/result[-1]) ** (1.0/(n-len(result)))
			# round, re-adjust to 0 indexing (i.e. minus 1) and return np.uint64 array
	return np.array(map(lambda x: round(x)-1, result), dtype=np.int)

def drop_outliers(df, lower_limit, upper_limit, ws, std):
	"""
	Routine to drop outliers that survive the range limiter via standard deviation
	Input:
		df: The data that you want to clean
		lower_limit: The lower limit in the range limiter, so any value lower
			than this value results in NAN
		upper_limit: The upper limit in the range limiter, so any value greater
			than this value results in NAN
		ws: The window size to slide through df to remove outliers that survive the range limiter
		std: The multiplier of the standard deviation that if data points are outside of this result,
			they become NAN
	Output:
		df: The cleaned DataFrame
	"""
# make a copy to not overwrite the original
	df = df.copy()
# conditionally check all values in df to be out of specified range
	df[(df < lower_limit) | (df > upper_limit)] = np.NaN
# run through df in specified ranges set by plus_minus and then check all values 
# in range for data over specified factors of standard deviation
	for i in np.arange(ws, df.size-ws-1,1):
		IDX = np.s_[i-ws:i+ws]
		df[IDX][(np.abs(df[IDX]-df[IDX].mean())) >= (std*df[IDX].std())] = np.NaN
	return df


def PVI(ax,ay,az,lag):
	"""
	Routine to calculate the Partial Variance Increment between three time series, normalized
	Input:
		ax: series 1
		ay: series 2
		az: series 3
		lag: The number of points to shift each of the series ax, ay, and az
	Output:
		mag: normalized pvi
	"""
	ax = ax.copy()
	ay = ay.copy()
	az = az.copy()
# take the three factors of mag_field and find the derivatives
	dax=ax.shift(-lag)-ax
	day=ay.shift(-lag)-ay
	daz=az.shift(-lag)-az
# calc the magnitude of the derivatives of mag_field and return the pvi
	mag = dax.pow(2)+day.pow(2)+daz.pow(2)
	return (mag.div(mag.mean())).pow(.5)
    
def structure(ax, ay, az, ar_lags, ar_powers):
	"""
	Routine to compute the Structure coefficients of a certain series or number of series to different orders
	of structure at lags given in ar_lags
	Input:
		ax: series 1
		ay: series 2
		az: series 3
		ar_lags: The array consisting of lags, being the number of points to shift each of the series
		ar_powers: The array consisting of the Structure orders to perform on the series for all lags
	Output:
		df: The DataFrame containing the structure coefficients corresponding to ar_lags for each order in ar_powers 
	"""
# run through ar_lags and ar_powers so that for each power, run through each lag
	df = {}
	ax = ax.copy()
	ay = ay.copy()
	az = az.copy()

	for i in ar_powers:
		array = []
		for l in ar_lags:
			dax=ax.shift(-l)-ax
			day=ay.shift(-l)-ay
			daz=az.shift(-l)-az
			strct = (dax.pow(2)+day.pow(2)+daz.pow(2)).pow(0.5).pow(i).mean()
			array += [strct]
		df[str(i)] = array

	df = pd.Series(df)
	return df

def diToDeltaT(lags, v_mean, n_mean):
	"""
	Routine to convert a series of lags in di, to lags in seconds
	based on Taylor's hypothesis. 
	Input: 
		lags: lags in di to be converted to seconds
		v_mean: mean velocity of interval in km/s
		n_mean: mean density in cm^-3
	Output:
		lags: lags in seconds
	"""
	return (lags * 228 / v_mean / np.sqrt(n_mean))

def lambdaToDeltaT(lags, v_mean, n_const, au):
	"""
	Routine to convert a series of lags as a constant fraction of Lambda, to lags in seconds
	based on Taylor's Hypothesis
	Input:
		lags: lags as a constant fraction of Lambda to be converted to seconds
		v_mean: mean velocity of interval in km/s
		n_mean: mean density in cm^-3
		au: The AU range number containing the chosen interval
	Output:
		lags: lags in seconds
	"""
	di1 = 228 / np.sqrt(n_const)
	return (lags * di1 * np.sqrt(au) / v_mean)


def kurtosis(series, ptlag):
    """
    Routine to perform the Kurtosis on a series at a certain lag
    Input:
	series: the series to take the kurtosis of
	ptlag: the array of lags that represents the number of points to shift the series by
    Output:
	k: the kurtosis of the series
	ptlag: same array as input
    """
    k = []
    ptlag = np.asarray(ptlag).astype(int)
    series = series.copy()
    for i in ptlag:
        temp = (series.shift(-i) - series).copy()
        if (temp.pow(2).mean())**2 == 0 or (temp.pow(4).mean()) == 0:
            coeff = np.NaN
        else:
            coeff = temp.pow(4).mean()/(temp.pow(2).mean()**2)
        k += [coeff]
    k = pd.Series(k)

    return k, ptlag


def pdf(series, binsize):
    """
    Routine to determine the Probability Density Function of a certain Series 
    with a determined bin size
    Input:
    	series: The series to have the PDF performed on
    	binsize: The size of each bin in the PDF
    Output:
    	bins: the mean values of each bin in the pdf
    	pdf: the probability density of the series weighted to bins of size binsize.
    		The bigger the binsize, the lower number of bins, and vice versa
    """
# find rms, create empty array of arrays for bins
    series = series.copy()
#dropna and then sort the data from min to max, then reset the indices
    series.dropna(inplace=True)
#    rmsval = series.std()
#    series = series.div(rmsval)
# find the rms value of the series
    
    series.sort_values(inplace=True)
    series.reset_index(drop=True, inplace = True)
    length = len(series)/binsize
    pdf = np.zeros(length); bins=np.zeros(length)
# For each bin, take the size, divide by the max-min of that bin, then add to pdf
    acc = 0
    for i in range(binsize,len(series),binsize):
        temp = series[i-binsize:i]
        bins[acc] = temp.mean()
        pdf[acc]  = binsize/(temp.max()-temp.min())
        acc += 1
# return array of pdf
    return bins,pdf/len(series)


def waiting_time(series, time_var, ar_pvi):
	"""NOT FINISHED YET
	Routine to determine the time between two occurrences of the same PVI values
	Input: 
		series: The series containing the PVI data
		dt: The time between each point in data set that was used in calculation series
		ar_pvi: The PVI value(s) you want to determine the waiting times for
	Output:
		times: A series containing all time differences between each PVI occurrence
			for each PVI value provided in ar_pvi  
	"""
	series[0] = series.copy()
# to determine whether to look for pvi greater or less than desired pvi_lim
	if (constraint == 'greater'):
                cropped = series[0][series[0] > pvi_lim]
	else: 
                cropped = series[0][series[0] < pvi_lim]
# Calculate the time difference from one point to another       
	times = series[1].shift(-1) - series[1]
# return both the new cropped series and time intervals between each pvi
	return cropped, times

def get_dist_from_times_with_avg_speed(series, ar_time, ar_idx):
    """
    """
    avg_vels = []
    sec_times = []
    ar_dist = []
    for i in range(ar_idx) -1:
# temp for the data points bewteen the two occurences of the same pvi
        temp = series[ar_idx[i]:ar_idx[i+1]]
# find the avg velocity between the two occurences of the same pvi
        avg_vels += [temp.mean()]
# find the time in seconds between the two occurences of the same pvi
        sec_times += [(ar_idx[i+1]-ar_idx[i])*2]
# calculate the distance
        ar_dist += [avg_vels[i]*sec_times[i]]

    dictionary = {'vels':np.asarray(avg_vels), 'secs':np.asarray(sec_times), 'dist':np.asarray(ar_dist)}
    
    return pd.DataFrame(dictionary)
    
    
        
    
    
    
def correlation_coefficient(lagged_x, y, lag):
	"""
	Routine to determine the correlation coefficient between two data series
	with a certain lag shift on one of the series
	Input:
		lagged_x: The series that will be shifted by a certain lag
		y: The series that will not be affected by the lag shift
		lag: The number of points that lagged_x will be shifted by
	Output:
		top: the coefficient that will be normalized in correlation()
	"""
# Take copies of the time series
	lagged_x = lagged_x.copy()
	y = y.copy()
# take the shift of the spec'd series
	x = lagged_x.shift(-lag)
	ymean = np.zeros(len(y))
	xmean = np.zeros(len(x))
	ymean[0:] = y.mean()
	xmean[0:] = x.mean()
	ymean = pd.Series(ymean)
	xmean = pd.Series(xmean)
	top = ((x-xmean)*(y-ymean)).sum()
	return top

def correlation(lagged_x, y, ar_lags):
    """
    Routine to compute the Normalized Correlation between two arrays of points
    as one of the arrays is shifted by a series of lags
    Input:
    	lagged_x: The series that will be shifted by a certain lag
    	y: The Series that will not be affected by the lag shift
    	ar_lags: The array of lags that lagged_x will be shifted by
    Output:
    	corr: The series of normalized correlation coefficients each index
    		corresponding to the lag found at the same index from ar_lags
    """
    corr = np.zeros(np.size(ar_lags))
    acc = 0
    for i in range(np.size(ar_lags)):
        x = lagged_x.copy().shift(-ar_lags[i])
        ymean = np.zeros(len(y))
        xmean = np.zeros(len(x))
        ymean[0:] = y.mean()
        xmean[0:] = x.mean()
        ymean = pd.Series(ymean)
        xmean = pd.Series(xmean)
        print(np.sqrt((lagged_x-xmean).pow(2).sum() * (y-ymean).pow(2).sum()))
        corr[acc] = correlation_coefficient(lagged_x, y, ar_lags[i]) / np.sqrt((x-xmean).pow(2).sum() * (y-ymean).pow(2).sum())
        acc += 1
    return corr


def auto_correlation(series, ar_lags):
	"""
	Routine to take the correlation of a series with itself at different lag shifts
	Input:
		series: the series that will be correlated with itself
		ar_lags: The array of lags that series will be shifted by
	Output:
		corr: The series of normalized correlation coefficients, each index
			corresponding to the lag found at the same index from ar_lags
	"""
	return correlation(series, series, ar_lags)



def convert_cdfs_to_dataframe(filelist, varlist, time_var_str):
    """
    Routine to conver cdfs to a dictionary of arrays, keys are strings from varlist
    Input:
	filelist: the filename/paths to the cdf files you want to read in and convert to arrays
	varlist: the variable names you want to read in from the cdf files
	time_var_str: the string variable for the time you want from the file and also converted to datetime
    Output:
	dictionary: the dictionary, with keys that lead to their appropriate arrays from the cdf files
    """
    from spacepy import pycdf
    from delorean import Delorean
    
#create empty numpy arrays
    ll=len(varlist); varsdata=[np.zeros(1) for i in range(ll+1)]
# read data from cdf files and append the arrays.
    for i in filelist:
      #  print 'reading file '+i
        d = pycdf.CDF(i)
	# CHECK TIME STAMPS
	# FIND OUT WHAT BEGINNING INDEX TO USE
	# save last time stamp to compare with next file
        ctr = 0
        if i == filelist[0]:
            last = pycdf.VarCopy(d[time_var_str])[-1]
            
        else: 
            while (Delorean(pycdf.VarCopy(d[time_var_str])[ctr],timezone='UTC').epoch - Delorean(last,timezone='UTC').epoch < .01):
                ctr += 1
            last = pycdf.VarCopy(d[time_var_str])[-1]
        for j in varlist:
            idx=varlist.index(j)
            if j != 'scDistance':
                varsdata[idx] = np.append(varsdata[idx], pycdf.VarCopy(d[j])[ctr:])
            else:
                varsdata[idx] = np.append(varsdata[idx], pycdf.VarCopy(d[j])[int(ctr/25):])
    print('Done reading data')
#   For create an epoch array from time_var_str
#   (s)econds (s)ince (epoch) == ssepoch
    idxe = varlist.index(time_var_str); ldata=len(varsdata[0]); ssepoch=np.zeros(ldata)
    for i in range(1,ldata):
    	ssepoch[i] = Delorean(varsdata[idxe][i],timezone="UTC").epoch 
# drop the first zero before creating the data frame
    dictionary = {}; dictionary['time']=ssepoch[1:]
    for j in varlist:
        if j == time_var_str:
            dictionary['datetime']=varsdata[varlist.index(j)][1:]
        else:
            dictionary[j] = varsdata[varlist.index(j)][1:]
    return dictionary



def curve_fitting_poly(x, y, power):
    """
    Routine for polynomial fitting
    Input:
	x: Data along x axis
	y: Data along y axis
	power: the order of the polynomial you want the fit to model
    Output:
	x: same x as input
	y: same y as input
	y_new: the fitted values against 
	z: fitting parameters depending on inputs
    """
    z = np.polyfit(x, y, power)
    f = np.poly1d(z)

    y_new = f(x)

    return x, y, y_new, z


def lookup_indices_for_distance(time_val_arr, time):
	"""
	Routine to return the indices for when an AU begins/ends
	Input:
		time_val_arr: The Times at which AUs begin/end
		time: The Time series to look up the time values from time_val_arr
	Output:
		idx_aus: the index at which the times from time_val_arr occur in time
	"""
	idx_aus = [0]
	for val in time_val_arr:
		idx_aus += [np.asarray(time).index(val)]
	return idx_aus


