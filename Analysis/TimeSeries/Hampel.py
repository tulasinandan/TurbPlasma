#!/usr/bin/env python

# Apply Hampel filter
# Then recover the time series

import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
import scipy.fftpack
import cmath as cm
import math
import sys
import os
import readline
import re

##   # Delete some points from the start and end
##   def cut_edge(xin, ncut):
##   
##      xout=np.zeros((len(xin[:,0])-int(2.0*ncut), len(xin[0,:])))
##   
##      ## Delete the first ncut and last ncut elements
##      for i in range(len(xin[0,:])):
##         if ncut == 0:
##            xout[:,i] = xin[:,i]
##         else:
##            xout[:,i] = xin[ncut:-ncut,i]
##   
##      return xout  

# Check NaNs in time series and replace with linear interpolation
def check_NaN(xin):

   xout = xin.copy()

   for i in range(len(xin[0,:])):
      for j in range(len(xin[:,0])):
         if math.isnan(xin[j,i])==True:
            print("NaN at", i, j)
            # Check if the element just before is a number       
            if math.isnan(xin[j-1,i])==False:
               istart=i; jstart=j-1
            # Check if the element just after is a numer
            if math.isnan(xin[j+1,i])==False:
               iend=i; jend=j+1
               x = np.array([xin[jstart,0], xin[jend,0]])
               fp = np.array([xin[jstart,i], xin[jend,i]])
               f = interpolate.interp1d(x,fp)
               xout[jstart+1:jend,i] = f(xin[jstart+1:jend,0])


   # Finally check if xout is NaN free
   for i in range(len(xout[0,:])):
      for j in range(len(xout[:,0])):
         if math.isnan(xout[j,i])==True:
            print("v array has NaN!!", i, j)

   return xout

def hampel(x,k, t0):
    '''adapted from hampel function in R package pracma
    x= 1-d numpy array of numbers to be filtered
    k= number of items in window/2 (# forward and backward wanted to capture in median filter)
    t0= number of standard deviations to use; 3 is default
    '''
    n = len(x)
    y = x.copy() #y is the corrected series
    L = 1.4826
    for i in range((k + 1),(n - k)):
        if np.isnan(x[(i - k):(i + k+1)]).all():
            continue
        x0 = np.nanmedian(x[(i - k):(i + k+1)])
        S0 = L * np.nanmedian(np.abs(x[(i - k):(i + k+1)] - x0))
        if (np.abs(x[i] - x0) > t0 * S0):
  	         y[i] = x0

    return y

def spectral_hampel_scalar(vx,dt,k,t0):

   '''adapted from hampel function in R package pracma
    x= 1-d numpy array of numbers to be filtered
    k= number of items in window/2 (# forward and backward wanted to capture in median filter)
    t0= number of standard deviations to use; 3 is default
   '''
   nt=len(vx)
   ## Do FFT
   f = np.linspace(0.0, 1.0/(2.0*dt), nt/2)	# freq axis
   df = f[1]-f[0]

   # Presuming v is vector, compute its spectrum
   ux = scipy.fftpack.fft(vx)
   usqd=abs(ux)**2
   
#  print "Energy in freq space, unprocessed data", (1.0/nt)*np.sum(usqd)
   
   
   ############## STAGE 2, HAMPEL FILTER ###############
   
   mx=np.abs(ux)
   ax=np.angle(ux)
   
   n = len(usqd)
   usqdf = usqd.copy() # usqdf is the corrected series
   uxf = mx.copy()

   L = 1.4826

   for i in range((k + 1),(n - k)):
      if np.isnan(usqd[(i - k):(i + k+1)]).all():
            continue
      usqd0 = np.nanmedian(usqd[(i - k):(i + k+1)])
      S0 = L*np.nanmedian(np.abs(usqd[(i - k):(i + k+1)] - usqd0))
      ux0 = np.nanmedian(mx[(i - k):(i + k+1)])
    # print usqd[i],usqd0
      if (np.abs(usqd[i] - usqd0) > t0 * S0):
  	      usqdf[i] = usqd0	
	      uxf[i] = ux0
   # Go back to real space
   vx2 = np.real(scipy.fftpack.ifft( uxf*(np.cos(ax)+cm.sqrt(-1)*np.sin(ax))))

   return vx2

def spectral_hampel_vector(vx,vy,vz,dt,k,t0):

   '''adapted from hampel function in R package pracma
    x= 1-d numpy array of numbers to be filtered
    k= number of items in window/2 (# forward and backward wanted to capture in median filter)
    t0= number of standard deviations to use; 3 is default
   '''
   nt=len(vx)
   ## Do FFT
   f = np.linspace(0.0, 1.0/(2.0*dt), nt/2)	# freq axis
   df = f[1]-f[0]

   # Presuming v is vector, compute its spectrum
   ux = scipy.fftpack.fft(vx)
   uy = scipy.fftpack.fft(vy)
   uz = scipy.fftpack.fft(vz)
   usqd=abs(ux)**2+abs(uy)**2+abs(uz)**2
   
#  print "Energy in freq space, unprocessed data", (1.0/nt)*np.sum(usqd)
   
   
   ############## STAGE 2, HAMPEL FILTER ###############
   
   mx=np.abs(ux)
   ax=np.angle(ux)
   my=np.abs(uy)
   ay=np.angle(uy)
   mz=np.abs(uz)
   az=np.angle(uz)
   
   n = len(usqd)
   usqdf = usqd.copy() # usqdf is the corrected series
   uxf = mx.copy()
   uyf = my.copy()
   uzf = mz.copy()

   L = 1.4826

   for i in range((k + 1),(n - k)):
      if np.isnan(usqd[(i - k):(i + k+1)]).all():
            continue
      usqd0 = np.nanmedian(usqd[(i - k):(i + k+1)])
      S0 = L*np.nanmedian(np.abs(usqd[(i - k):(i + k+1)] - usqd0))
      ux0 = np.nanmedian(mx[(i - k):(i + k+1)])
      uy0 = np.nanmedian(my[(i - k):(i + k+1)])
      uz0 = np.nanmedian(mz[(i - k):(i + k+1)])
      if (np.abs(usqd[i] - usqd0) > t0 * S0):
  	      usqdf[i] = usqd0	
	      uxf[i] = ux0
	      uyf[i] = uy0
	      uzf[i] = uz0
   # Go back to real space
   vx2 = np.real(scipy.fftpack.ifft( uxf*(np.cos(ax)+cm.sqrt(-1)*np.sin(ax))))
   vy2 = np.real(scipy.fftpack.ifft( uyf*(np.cos(ay)+cm.sqrt(-1)*np.sin(ay))))
   vz2 = np.real(scipy.fftpack.ifft( uzf*(np.cos(az)+cm.sqrt(-1)*np.sin(az))))
	    
   return vx2,vy2,vz2


