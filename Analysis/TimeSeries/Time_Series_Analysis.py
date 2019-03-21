#!/usr/bin/env python
import sys, os
sys.path.insert(0,os.environ['HOME']+'/AJGAR/TurbPlasma')
import random as rnd
import numpy as np
import numpy.fft as nf
from ..Simulations import AnalysisFunctions as af
import pandas as pd
#import seaborn as sns
import matplotlib.pyplot as plt
from .OLLibs import F90 as ftsa

pi = np.pi

def smooth_t(a,nums):
   b=a.copy()
   for i in range(nums):
      b[1::][:-1] = 0.25*b[:-2]+0.5*b[1::][:-1]+0.25*b[2:]
   return b

def gen_log_space(limit, n):
    result = [1]
    if n>1:  # just a check to avoid ZeroDivisionError
        ratio = (float(limit)/result[-1]) ** (1.0/(n-len(result)))
    while len(result)<n:
        next_value = result[-1]*ratio
        if next_value - result[-1] >= 1:
            # safe zone. next_value will be a different integer
            result.append(next_value)
        else:
            # problem! same integer. we need to find next_value 
            # by artificially incrementing previous value
            result.append(result[-1]+1)
            # recalculate the ratio so that the remaining values 
            # will scale correctly
            ratio = (float(limit)/result[-1]) ** (1.0/(n-len(result)))
    # round, re-adjust to 0 indexing (i.e. minus 1) and return np.uint64 array
    return np.array([round(x)-1 for x in result], dtype=np.int)
    #return np.array(list(map(lambda x: round(x)-1, result)), dtype=np.uint64)

# Correlation functions
def corr(x, y,lags,dt):
   return ftsa.correlation(x,y,lags,dt)
def autocorr(x,lags,dt):
   return ftsa.correlation(x,x,lags,dt)
def autocorrvec(x,y,z,lags,dt):
   r,rxx=autocorr(x,lags,dt)
   r,ryy=autocorr(y,lags,dt)
   r,rzz=autocorr(z,lags,dt)
   return r, (rxx+ryy+rzz)/3.

def lowpass_sharp(a,dt,freq):
   nt=len(a)
  #a=np.array(a)
  #an=np.append(np.append(a[1:nt/2][::-1],a),a[nt/2:-1][::-1])*af.windowff(2*nt-2,kind='tukey')
  #fa=nf.rfft(an)
  #fq=nf.rfftfreq(2*nt-2,d=dt)
  #idx=np.argmin(np.abs(fq-freq))
  #fa[idx:]=0j
  #b=nf.irfft(fa)[nt/2-1:-nt/2+1]
   fa=nf.rfft(a)
   fq=nf.rfftfreq(nt,d=dt)
   idx=np.argmin(np.abs(fq-freq))
   fa[idx:]=0j
   b=nf.irfft(fa)
   return b
##~   import numpy as np
##~   from scipy.signal import butter, lfilter, freqz
##~   import matplotlib.pyplot as plt
##~   
##~   
##~   def butter_lowpass(cutoff, fs, order=5):
##~       nyq = 0.5 * fs
##~       normal_cutoff = cutoff / nyq
##~       b, a = butter(order, normal_cutoff, btype='low', analog=False)
##~       return b, a
##~   
##~   def butter_lowpass_filter(data, cutoff, fs, order=5):
##~       b, a = butter_lowpass(cutoff, fs, order=order)
##~       y = lfilter(b, a, data)
##~       return y
##~   return butter_lowpass_filter(a,freq,1/dt)

def highpass_sharp(a,dt,freq):
   nt=len(a)
   fa=nf.rfft(a)
   fq=nf.rfftfreq(nt,d=dt)
   idx=np.argmin(np.abs(fq-freq))
   fa[:idx]=0j
   b=nf.irfft(fa)
   return b

def btspec(a_in,dt,downsample=10):
   r,dac=autocorr(a_in,list(range(int(len(a_in)/downsample))),dt)
   dac=np.append(dac[::-1][1:],dac[1:])
#  fda=nf.rfft(dac*af.windowff(len(dac),kind='tukey'))/len(dac)
   fda=nf.rfft(dac)/len(dac)
   fq=nf.rfftfreq(len(dac),d=dt)
   spec=0.5*abs(fda)**2 
   return fq,spec

def bmspec(a_in,dt,downsample=2.):
   """
      Code to compute power spectrum using the technique 
      from Bieber & Matthaeus, JGR 1993.
   """
   da=a_in[1:]-a_in[:-1]
   fq,spec=btspec(da,dt,downsample=downsample)
   filt=np.zeros_like(fq)
   filt[1:]=(np.pi*fq[1:]*dt)**2 / (4 * np.sin(np.pi*fq[1:]*dt)**4)
   return fq,spec*filt

def fspec(a_in,dt,winkind='tukey'):
   """
      Code to compute power spectrum using simple
      minded Fourier transform on windowed data.
   """
   a=a_in*af.windowff(len(a_in),kind=winkind)
   fa=nf.rfft(a)/len(a)
   fq=nf.rfftfreq(len(a),d=dt)
   spec=abs(fa)**2 
   return fq,spec

def sdk(series, lags, dt):
   k = []; tau=[]
   series = series.copy()
   for i in lags:
       temp = (series.shift(-i) - series).copy()
       coeff = temp.pow(4).mean()/(temp.pow(2).mean()**2)
       k   += [coeff]
       tau += [i*dt]
   k = pd.Series(k)
   tau = pd.Series(tau)
   return tau, k

def calc_pdf(series,weight=100,inc=None,Normalize=False):
# find rms, create empty array of arrays for bins
   if inc is not None:
      series = (series-series.shift(inc))
   else:
      series=series.copy()
#dropna and then sort the data from min to max, then reset the indices
   series.dropna(inplace=True)
   if Normalize:
      rmsval = series.std()
      series = series/rmsval
# find the rms value of the series
    
   series.sort_values(inplace=True)
   series.reset_index(drop=True, inplace = True)
   length = int(len(series)/weight)
   pdf = np.zeros(length); bins=np.zeros(length)
# For each bin, take the size, divide by the max-min of that bin, then add to pdf
   acc = 0
   for i in range(weight,len(series),weight):
      temp = series[i-weight:i]
      bins[acc] = temp.mean()
      pdf[acc]  = weight/(temp.max()-temp.min())
      acc += 1
# return array of pdfs
   return bins,pdf/len(series)

def pdf_inc(series,inc,weight=100):
   tmp=(series.shift(-inc)-series).copy()
   return calc_pdf(tmp,weight)
