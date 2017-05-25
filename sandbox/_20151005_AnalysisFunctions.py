#!/usr/bin/env python
from pylab import *
from numpy.fft import *
import operator

##
## DEF TO CALCULATE THE SPECTRUM. IT GIVE PERP SPECTRUM PRESUMING
## THAT THE MEAN FIELD IS IN THE Z DIRECTION. NOT TOO GENERAL BUT 
## WE'LL GENERALIZE IT IF/WHEN NEED BE.
##
def spectrum(ar):
   if len(ar) == 0:
      print 'No array provided! Exiting!'
      return
   ar=ar-mean(ar)
   nx=shape(ar)[0];kx=fftshift(fftfreq(nx))*nx
   ny=shape(ar)[1];ky=fftshift(fftfreq(ny))*ny
   nz=shape(ar)[2];kz=fftshift(fftfreq(nz))*nz
   kp=np.zeros((nx,ny))
   fekp=np.zeros(min(nx/2,ny/2))
   for x in range(nx):
      for y in range(ny):
         kp[x,y]=sqrt(kx[x]**2+ky[y]**2)
   dk=abs(kp[nx/2,1]-kp[nx/2,0])
  
   far = fftshift(fftn(ar))/(nx*ny*nz); fftea=0.5*abs(far)**2
   ffteb=sum(fftea,axis=2)
#  for i in range(size(fekp)):
#     fekp[i]= np.sum(np.ma.MaskedArray(ffteb, ~((kp[nx/2,i+ny/2]-dk < kp) & (kp < kp[nx/2,i+ny/2]+dk))))
   for x in range(nx):
      for y in range(ny):
         i=np.round(kp[x,y])
         if i <nx/2:
            fekp[i]=fekp[i]+ffteb[x,y]

   return kp[nx/2,ny/2:],fekp

##
## SPECTRUM OF A VECTOR
##
##
def specvec(arx,ary,arz):
   kwave,ekx = spectrum(arx)
   kwave,eky = spectrum(ary)
   kwave,ekz = spectrum(arz)
   return kwave,ekx,eky,ekz,ekx+eky+ekz

def autocorrelation(ar,ax=0,dx=1.):
   nlen=shape(ar)[ax]/2
   r=np.zeros(nlen);corr=np.zeros(nlen)
   for i in range(nlen):
      ars=np.roll(ar,i,axis=ax)
      corr[i]=mean(ar*ars)
      r[i]=i*dx
   corr = corr/mean(ar**2)
   return r,corr

def autocvec(xx,yy,zz,ax=0,dx=1.):
   if ax == 0:
      nlen = shape(xx)[0]/2
   elif ax == 1:
      nlen = shape(xx)[1]/2
   elif ax == 2:
      nlen = shape(xx)[2]/2
   r=np.zeros(nlen);corr=np.zeros(nlen)
   tmp=np.zeros(shape(xx))
   for i in range(nlen):
      xxs=np.roll(xx,i,axis=ax)
      yys=np.roll(yy,i,axis=ax)
      zzs=np.roll(zz,i,axis=ax)
      corr[i]=mean(xx*xxs+yy*yys+zz*zzs)
      r[i]=i*dx
   corr = corr/mean(xx**2+yy**2+zz**2)
   return r,corr

def correlation(ar1,ar2,ax=0,dx=1.):
   nlen=shape(ar)[ax]/2
   r=np.zeros(nlen);corr=np.zeros(nlen)
   for i in range(nlen):
      ars=np.roll(ar2,i,axis=ax)
      corr[i]=mean(ar1*ars)
      r[i]=i*dx
   corr = corr/mean(ar1*ar2)
   return r,corr

def calc_pdf(ar,min=99999,max=99999,weight=100,inc=0,ax=0):
   if len(ar) == 0:
      print 'No array provided! Exiting!'
      return
   if min == 99999:
      min=ar.min()
   if max == 99999:
      max=ar.max()
# If PDF of increment, then increment the array
   if inc > 0:
      ar = ar - np.roll(ar,inc,axis=ax)
# Find the total length of data set
   arsize=reduce(operator.mul, shape(ar),1)
# Find the RMS value of data set and normalize to it.
   rmsval = sqrt(mean(ar**2))
   if rmsval != 0:
      ar = ar/rmsval
# Reshape the array to 1D & sort it.
   arr=np.reshape(ar,arsize)
   np.ndarray.sort(arr,kind='heapsort')
# Empty arrays for output
   bins=int(arsize/weight); pdf=np.zeros(bins); binvals=np.zeros(bins)
# Fill the bins 
   for i in range(bins):
      start=i*weight
      binvals[i] = mean(arr[start:start+weight])
      pdf[i] = weight/(arr[start:start+weight].max()-arr[start:start+weight].min())
   pdf = pdf/arsize
   return binvals,pdf

def normhist(ar,min=99999,max=99999,nbins=100,inc=0,ax=0):
   if len(ar) == 0:
      print 'No array provided! Exiting!'
      return
# If PDF of increment, then increment the array
   if inc > 0:
      ar = ar - np.roll(ar,inc,axis=ax)
# Find the total length of data set
   arsize=reduce(operator.mul, shape(ar),1)
# Find the RMS value of data set and normalize to it.
   rmsval = sqrt(mean(ar**2))
   if rmsval != 0:
      ar = ar/rmsval
# Reshape the array to 1D & sort it.
   arr=np.reshape(ar,arsize)
   np.ndarray.sort(arr,kind='heapsort')
# Empty arrays for output
   if min == 99999:
      min=ar.min()
   if max == 99999:
      max=ar.max()
   bins=np.linspace(min,max,nbins)
   dbin=bins[1]-bins[0]
   normhist=np.zeros(nbins)
# Rescale everything to have zero at min
   arr = arr-min
# Fill the bins 
   for i in range(len(arr)):
      j = floor(arr[i]/dbin)
      normhist[j] = normhist[j]+1
   normhist = normhist/(arsize*dbin)
   return bins,normhist

def kurtosis(ar):
   return mean(ar**4)/mean(ar**2)**2

def windowff(nslices,kind=""):
   if len(kind) == 0:
      kind = "blackmannharris"

   windowff=np.zeros(nslices)
   for i in range(nslices):
      tht = 2*pi*i/(nslices-1)
#  Hanning
      if kind.lower() == "hanning":
         windowff[i]=0.5*(1-cos(tht))
#  Blackman
      if kind.lower() == "blackmann":
         windowff[i]=0.42659-0.49656*cos(tht)+0.076849*cos(2*tht)
#  BlackmanHarris
      if kind.lower() == "blackmannharris": 
         windowff[i]=0.35875-0.48829*cos(tht)+0.14128*cos(2*tht)-0.01168*cos(3*tht)
#  Flat top Window
      if kind.lower() == "flattop": 
         windowff[i] = 1. - 1.93*cos(tht) + 1.29*cos(2*tht) - 0.388*cos(3*tht) + 0.028*cos(4*tht)
#  Nuttall
      if kind.lower() == "nuttall": 
         windowff[i] = 0.355768 - 0.487396*cos(tht) + 0.144232*cos(2*tht) - 0.012604*cos(3*tht)
#  Welch
      if kind.lower() == "welch": 
         windowff[i] = 1. - ((i - (nslices-1)/2.)/((nslices+1)/2.))**2
# Tukey
#     if kind.lower() == "tukey":
#        if i < int(0.1*(nslices-1)/2):
#           windowff[i] = 0.5*(1+cos(tht/0.1 - pi))
#        if int(0.1*(nslices-1)/2) < i and i < int((nslices-1)*0.95):
#           windowff[i] = 1
#        if i > int((nslices-1)*0.95):
#           windowff[i] = 0.5*(1+cos(tht/0.1 - 2*pi/0.1 +pi))

#
# THE FOLLOWING ROUTINE CALCULATES THE SCALE DEPENDENT
# KURTOSIS OF A GIVEN ARRAY IN THE GIVEN DIRECTION
#
def sdk(ar,ax=0):
   nx=shape(ar)[ax]
   ark=np.zeros(nx/2)
   ddx=np.arange(nx/2)
   for i in range(1,nx/2):
      tmp=np.roll(ar,i,axis=ax)      
      tmp=ar-tmp; rmsval=sqrt(mean(tmp**2))
      if rmsval != 0:
         tmp=tmp/rmsval 
      ark[i] = kurtosis(tmp)
   return ddx,ark
#
# THE FOLLOWING ROUTINE READS A GENERIC 3D ARRAY FROM A
# DIRECT ACCESS UNFORMATTED FILE AT A GIVEN SNAPSHOT
#
def readsu(ff,nx,ny,nz,timeslice):
   ff.seek(8*timeslice*nx*ny*nz)
   field = np.fromfile(ff,dtype='float64',count=nx*ny*nz)
   field = np.reshape(field,(nx,ny,nz),order='F')
   return field
#
# Savitzky Golay filter for smoothing signals.
# http://wiki.scipy.org/Cookbook/SavitzkyGolay
# http://en.wikipedia.org/wiki/Savitzky%E2%80%93Golay_filter
#
#
def savitzky_golay(y, window_size, order, deriv=0, rate=1):
    r"""Smooth (and optionally differentiate) data with a Savitzky-Golay filter.
    The Savitzky-Golay filter removes high frequency noise from data.
    It has the advantage of preserving the original shape and
    features of the signal better than other types of filtering
    approaches, such as moving averages techniques.
    Parameters
    ----------
    y : array_like, shape (N,)
        the values of the time history of the signal.
    window_size : int
        the length of the window. Must be an odd integer number.
    order : int
        the order of the polynomial used in the filtering.
        Must be less then `window_size` - 1.
    deriv: int
        the order of the derivative to compute (default = 0 means only smoothing)
    Returns
    -------
    ys : ndarray, shape (N)
        the smoothed signal (or it's n-th derivative).
    Notes
    -----
    The Savitzky-Golay is a type of low-pass filter, particularly
    suited for smoothing noisy data. The main idea behind this
    approach is to make for each point a least-square fit with a
    polynomial of high order over a odd-sized window centered at
    the point.
    Examples
    --------
    t = np.linspace(-4, 4, 500)
    y = np.exp( -t**2 ) + np.random.normal(0, 0.05, t.shape)
    ysg = savitzky_golay(y, window_size=31, order=4)
    import matplotlib.pyplot as plt
    plt.plot(t, y, label='Noisy signal')
    plt.plot(t, np.exp(-t**2), 'k', lw=1.5, label='Original signal')
    plt.plot(t, ysg, 'r', label='Filtered signal')
    plt.legend()
    plt.show()
    References
    ----------
    .. [1] A. Savitzky, M. J. E. Golay, Smoothing and Differentiation of
       Data by Simplified Least Squares Procedures. Analytical
       Chemistry, 1964, 36 (8), pp 1627-1639.
    .. [2] Numerical Recipes 3rd Edition: The Art of Scientific Computing
       W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery
       Cambridge University Press ISBN-13: 9780521880688
    """
    import numpy as np
    from math import factorial

    try:
        window_size = np.abs(np.int(window_size))
        order = np.abs(np.int(order))
    except ValueError, msg:
        raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    order_range = range(order+1)
    half_window = (window_size -1) // 2
    # precompute coefficients
    b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
    m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
    lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    return np.convolve( m[::-1], y, mode='valid')

#
#
# FUNCTION TO CALCULATE PHASE SPEEDS OF THE THREE BRANCHES OF TWO FLUID
# DISPERSION RELATION (STRINGER 1963)
#
def tfps(beta=0.6, ca=1., de2=0.000545, theta=0., kk=1.):
   """
      beta: Total plasma beta,
      ca: Alfven speed based on mean field, 
      de2: me/mi, 
      theta: Angle of propagation as fraction of pi/2), 
      kk: Wavenumber of interest in units of kdi
      
      Output is frequencies of the roots and the phase speeds w/k
      The roots are w[0]:Fast/Whistler, w[1]:Alfven/KAW, w[2]: Slow/Cyclotron
   """
#  theta=float(raw_input("Angle of propagation as fraction of pi/2: ")) 
#  kk=float(raw_input("What k value? "))
#  ca = float(raw_input("Alfven Speed: "))
#  beta=float(raw_input("Total plasma beta?: "))
#  de2= float(raw_input("me/mi : "))

   ct = cos(theta*pi/2)
   st = sin(theta*pi/2)
   tt = st/ct
   cs=sqrt(beta/2.)*ca
   di =1.
   caksq=cos(theta*pi/2.)**2*ca**2
   cmsq=ca**2 + cs**2
   
   pcs=np.zeros(4)
   D = 1 + kk**2*de2
   # Find out the root of the quadratic dispersion relation
   pcs[0] = 1.
   pcs[1] = -(ca**2/D + cs**2 + caksq*(1+kk**2*di**2/D)/D)*kk**2
   pcs[2] = caksq*kk**4*(ca**2/D + cs**2 + cs**2*(1+kk**2*di**2/D))/D
   pcs[3] = -(cs**2*kk**6*caksq**2)/D**2
   w2 = np.roots(pcs); w = sqrt(w2)
   speeds= w/kk
   return w,speeds

def tfst(beta=0.6, ca=1., de2=0.000545, theta=0., kmin=1e-2, kmax=10., npoints=200,wrt='n'):
   """
      beta: Total plasma beta,
      ca: Alfven speed based on mean field, 
      de2: me/mi, 
      theta: Angle of propagation as fraction of pi/2), 
      kkmin: Minimum Wavenumber of interest in units of kdi
      kkmax: Maximum wavenumber of interest in units of kdi

      Output is an array with 4 columns, k, w-fast, w-alf, w-slow
   """
   
   kmmn=np.log10(kmin)
   kmmx=np.log10(kmax)
   kk=np.logspace(kmmn,kmmx,npoints)
   warray=np.zeros((3,npoints))
   for i in range(0,npoints):
      f,s = tfps(beta, ca, de2, theta, kk[i])
      warray[:,i]=f
   plt.loglog(kk,warray[0,:], label='Fast/Magnetosonic')
   plt.loglog(kk,warray[1,:], label='Alfven/KAW')
   plt.loglog(kk,warray[2,:], label='Slow')
   plt.xlabel('$kd_i$')
   plt.ylabel('$\omega/\omega_{ci}$')
   plt.legend(loc='best',fancybox=True,framealpha=0.2)
   plt.title('Dispersion Relation for beta='+str(beta)+' and me/mi='+str(de2))
   plt.show()
   if wrt == 'y':
      ofile=open('disp'+str(theta*pi/2.)+'.dat','w')
      print>> ofile,'#  k', 'Acoustic', 'Alfven', 'Fast'
      for i in range(npoints):
         print>> ofile, kk[i],warray[2,i],warray[1,i],warray[0,i]
      ofile.close()
