#!/usr/bin/env python
import numpy.fft as nf
import numpy as np
from functools import reduce
pi = np.pi
eye = 0+1j 
import operator
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter as gf
from numpy import newaxis as nna
import os
from .OLLibs.F90 import faf

#==================================================

def create_kgrid(nx, ny, nz, lx=2*pi, ly=2*pi, lz=2*pi):
   """
   Create a 3D k grid for Fourier space calculations
   """

   print(lx, ly, lz)

   kx = nf.fftshift(nf.fftfreq(nx))*nx*2*pi/lx
   ky = nf.fftshift(nf.fftfreq(ny))*ny*2*pi/ly
   kz = nf.fftshift(nf.fftfreq(nz))*nz*2*pi/lz
  
   mg = np.meshgrid(kx,ky,kz)

   km = np.sqrt(np.sum((m**2 for m in mg)))

   return kx[:,nna,nna], ky[nna,:,nna], kz[nna,nna,:], km

#==================================================

def kfilter(ar, kf, lx=2*pi, ly=2*pi, lz=2*pi):
   """
   Function to filter a 3D array by zeroing out larger k values
   """

   while len(ar.shape) < 3:
      ar = ar.reshape(ar.shape + (1,))

   #  COMPUTE THE ARRAY SIZE 
   kx, ky, kz, km = create_kgrid(*ar.shape, lx=lx, ly=ly, lz=lz)

   #  FOURIER TRANSFORM THE ARRAY 
   far = nf.fftshift(nf.fftn(ar))

   #  SET VALUES ABOVE kf AS 0+0i
   far = (np.sign(km - kf) - 1.)/(-2.)*far

   #  BACK TRANSFORM TO REAL SPACE
   arf = np.real(nf.ifftn(nf.ifftshift(far)))
   return arf

#==================================================

def kdiv(arx,ary,arz,kf=None,lx=2*pi,ly=2*pi,lz=2*pi):
   """
   Function to compute divergence in Fourier space
   """

   #  COMPUTE THE ARRAY SIZE
   kx, ky, kz, km = create_kgrid(*np.shape(arx), lx=lx, ly=ly, lz=lz)

   #  FOURIER TRANSFORM THE ARRAY
   far = [nf.fftshift(nf.fftn(a)) for a in (arx, ary, arz)]

   # COMPUTE div=i*(kx*ax+ky*ay+kz*az)
   mg = np.meshgrid(kx,ky,kz)

   divarf = 1.j*reduce(operator.add, [a*b for a,b in zip(mg, far)])

   #  SET VALUES ABOVE kf AS 0+0i if kf non zero
   if kf is not None:
      divarf = (np.sign(km - kf) - 1.)/(-2.)*divarf

   divar = np.real(nf.ifftn(nf.ifftshift(divarf)))

   return divar

def vecpot(arx,ary,arz, kf, lx=2*np.pi, ly=2*np.pi, lz=2*np.pi):
   """
   Function to compute vector potential of a 3D array
   """
   nx,ny,nz=arx.shape

   #  COMPUTE THE ARRAY SIZE
   kx, ky, kz, km = create_kgrid(*arx.shape, lx=lx, ly=ly, lz=lz)
  #kx, ky, kz, k2 = kx[:,nna,nna],ky[nna,:,nna],kz[nna,nna,:], km**2
   k2=km**2
   k2[nx/2,ny/2,nz/2]=1.

   #  FOURIER TRANSFORM THE ARRAY
   farx = nf.fftshift(nf.fftn(arx))
   fary = nf.fftshift(nf.fftn(ary))
   farz = nf.fftshift(nf.fftn(arz))

   #  SET VALUES ABOVE kf AS 0+0i
   farx = (np.sign(km - kf) - 1.)/(-2.)*farx
   fary = (np.sign(km - kf) - 1.)/(-2.)*fary
   farz = (np.sign(km - kf) - 1.)/(-2.)*farz

   #  FIND THE CORRESPONDING VECTOR POTENTIAL A = -ik x B /k^2
   axf = -eye*(ky*farz-kz*fary)/k2
   ayf = -eye*(kz*farx-kx*farz)/k2
   azf = -eye*(kx*fary-ky*farx)/k2

   #  BACK TRANSFORM TO REAL SPACE
   ax  = np.real(nf.ifftn(nf.ifftshift(axf)))
   ay  = np.real(nf.ifftn(nf.ifftshift(ayf)))
   az  = np.real(nf.ifftn(nf.ifftshift(azf)))
   return ax,ay,az

def kurl(arx,ary,arz, kf, lx=2*np.pi, ly=2*np.pi, lz=2*np.pi):
   kx, ky, kz, km = create_kgrid(*arx.shape, lx=lx, ly=ly, lz=lz)
   
   #  FOURIER TRANSFORM THE ARRAY
   farx = nf.fftshift(nf.fftn(arx))
   fary = nf.fftshift(nf.fftn(ary))
   farz = nf.fftshift(nf.fftn(arz))

   #  SET VALUES ABOVE kf AS 0+0i
   farx = (np.sign(km - kf) - 1.)/(-2.)*farx
   fary = (np.sign(km - kf) - 1.)/(-2.)*fary
   farz = (np.sign(km - kf) - 1.)/(-2.)*farz

   #  COMPUTE VORTICITY
   axf = eye*(ky*farz-kz*fary)
   ayf = eye*(kz*farx-kx*farz)
   azf = eye*(kx*fary-ky*farx)

   #  BACK TRANSFORM TO REAL SPACE
   wx  = np.real(nf.ifftn(nf.ifftshift(axf)))
   wy  = np.real(nf.ifftn(nf.ifftshift(ayf)))
   wz  = np.real(nf.ifftn(nf.ifftshift(azf)))
   return wx,wy,wz

def calc_psi_2d(bx,by,dx,dy):
   import numpy as np
# Calculating Psi
   psi = np.zeros(np.shape(bx))
   for k in range(1,np.shape(bx)[0]):
      psi[k,0] = psi[k-1,0] + by[k,0]
   for k in range(0,np.shape(bx)[0]):
      for l in range(1,np.shape(bx)[1]):
         psi[k,l] = psi[k,l-1] - bx[k,l]
   psi=psi*dx; psi=psi-np.mean(psi)
   return psi


##
## DEF TO CALCULATE THE PERP SPECTRUM. 
##
def PerpSpectrum(ar,sumax=2,lenx=2*pi,leny=2*pi,lenz=2*pi):
   """
      PerpSpectrum(ar,sumax=2,lenx=2*pi,leny=2*pi,lenz=2*pi)
      ar -> Array to compute the spectrum of
      sumax -> Axis of magnetic field direction. Right now only x,y,z = 0,1,2
      lenx,leny,lenz -> System size in x,y,z directions to take into 
                        account the anisotropy of system if any

      RETURNS:
      kk -> Wavenumber array
      fekp -> Spectrum of the array
   """
   if len(ar) == 0:
      print('No array provided! Exiting!')
      return
   ar=ar-np.mean(ar)
   nx=np.shape(ar)[0];kx=nf.fftshift(nf.fftfreq(nx))*nx*(2*pi/lenx)
   ny=np.shape(ar)[1];ky=nf.fftshift(nf.fftfreq(ny))*ny*(2*pi/leny)
   nz=np.shape(ar)[2];kz=nf.fftshift(nf.fftfreq(nz))*nz*(2*pi/lenz)
  
   far = nf.fftshift(nf.fftn(ar))/(nx*ny*nz); fftea=0.5*np.abs(far)**2
   ffteb=np.sum(fftea,axis=sumax)
## DEFINE A TEMPORARY XY PLANE
   if sumax==0:
      nnx=ny; nny=nz
      kkx=ky; kky=kz
   elif sumax==1:
      nnx=nx; nny=nz
      kkx=kx; kky=kz
   elif sumax==2:
      nnx=nx; nny=ny
      kkx=kx; kky=ky
## DEFINE THE KPERP VALUES AND CORRESPONDING SPECTRUM
   fekp=np.zeros(min(nnx/2,nny/2))
   kp=np.zeros((nnx,nny))
   for x in range(nnx):
      for y in range(nny):
         kp[x,y]=np.sqrt(kkx[x]**2+kky[y]**2)
   if nnx == 1:
      dk=np.abs(kp[0,1]-kp[0,0])
      kk=kp[0,nny/2:]
   elif nny == 1:
      dk=np.abs(kp[1,0]-kp[0,0])
      kk=kp[nnx/2:,0]
   else:
      dk=np.abs(kp[1,0]-kp[0,0])
      kk=kp[nnx/2,nny/2:]

   for i in range(len(fekp)):
      fekp[i]= np.sum(np.ma.MaskedArray(ffteb, ~((kp[nx/2,i+ny/2]-dk < kp) & (kp < kp[nx/2,i+ny/2]+dk))))
   #for x in range(nnx):
   #   for y in range(nny):
   #      i=np.round(kp[x,y]/dk)
   #      if i < len(fekp):
   #         fekp[i]=fekp[i]+ffteb[x,y]

   return kk,fekp/dk

##
## FORTRAN PERPENDICULAR SPECTRUM OF A VECTOR
##
def fperpspecvec(arx,ary,arz,sumax=2):
   """
      fperpspecvec(arx,ary,arz,sumax=2)
      Computes the spectrum of a vector assuming Lx=Ly=Lz=2*pi
      Uses Fortran perpspectrum function
   """
   global nk
   nx=np.shape(arx)[0]; ny=np.shape(arx)[1]; nz=np.shape(arx)[2]
   ## DEFINE A TEMPORARY XY PLANE
   if sumax==0:
      nk=np.min([ny,nz])/2
   elif sumax==1:
      nk=np.min([nx,nz])/2
   else: #sumax==2
      nk=np.min([nx,ny])/2
   #  print sumax, nx, ny, nz, nk+1
   kwave,ekx = faf.perpspec(arx,sumax,nk)
   kwave,eky = faf.perpspec(ary,sumax,nk)
   kwave,ekz = faf.perpspec(arz,sumax,nk)
   return kwave,ekx,eky,ekz,ekx+eky+ekz
##
## FORTRAN PERPENDICULAR SPECTRUM OF A VECTOR
## IN DOMAIN WITH VARIABLE LX, LY, LZ
##
def flperpspecvec(arx,ary,arz,sumax=2,lx=2*pi,ly=2*pi,lz=2*pi):
   """
      flperpspecvec(arx,ary,arz,sumax=2,lx=2*pi,ly=2*pi,lz=2*pi)
      Computes the spectrum of a vector in a domain with variable lengths
      Uses Fortran lperpspec function
   """
   global nk
   nx=np.shape(arx)[0]; ny=np.shape(arx)[1]; nz=np.shape(arx)[2]
   ## DEFINE A TEMPORARY XY PLANE
   if sumax==0:
      nk=np.min([ny,nz])/2
   elif sumax==1:
      nk=np.min([nx,nz])/2
   else: #sumax==2
      nk=np.min([nx,ny])/2
   #  print sumax, nx, ny, nz, nk+1
   kwave,ekx = faf.lperpspec(arx,sumax,nk,lx,ly,lz)
   kwave,eky = faf.lperpspec(ary,sumax,nk,lx,ly,lz)
   kwave,ekz = faf.lperpspec(arz,sumax,nk,lx,ly,lz)
   return kwave,ekx,eky,ekz,ekx+eky+ekz
##
## PERPENDICULAR SPECTRUM OF A VECTOR
##
def PerpSpecVec(arx,ary,arz,sumx=2,lx=2*pi,ly=2*pi,lz=2*pi):
   """
      PerpSpecVec(arx,ary,arz,sumx=2,lx=2*pi,ly=2*pi,lz=2*pi)
      arx,ary,arz -> Components of a vector to compute the spectrum of
      sumax -> Axis of magnetic field direction. Right now only x,y,z = 0,1,2
      lenx,leny,lenz -> System size in x,y,z directions to take into 
                        account the anisotropy of system if any

      RETURNS:
      kk -> Wavenumber array
      ekx,eky,ekz,ekx+eky+ekz -> Spectrum of the array
   """
   kwave,ekx = PerpSpectrum(arx,sumax=sumx,lenx=lx,leny=ly,lenz=lz)
   kwave,eky = PerpSpectrum(ary,sumax=sumx,lenx=lx,leny=ly,lenz=lz)
   kwave,ekz = PerpSpectrum(arz,sumax=sumx,lenx=lx,leny=ly,lenz=lz)
   return kwave,ekx,eky,ekz,ekx+eky+ekz

##
## DEF TO CALCULATE 2D SPECTRUM. 
##
def Spec2D(ar,ax=2,lenx=2*pi,leny=2*pi,lenz=2*pi):
   """
      Spec2D(ar,ax=2,lenx=2*pi,leny=2*pi,lenz=2*pi)
      
      2D spectrum of ar perpendicular to axis ax
      
   """
   if len(ar) == 0:
      print('No array provided! Exiting!')
      return
   ar=ar-np.mean(ar)
   nx=np.shape(ar)[0];kx=nf.fftshift(nf.fftfreq(nx))*nx*(2*pi/lenx)
   ny=np.shape(ar)[1];ky=nf.fftshift(nf.fftfreq(ny))*ny*(2*pi/leny)
   nz=np.shape(ar)[2];kz=nf.fftshift(nf.fftfreq(nz))*nz*(2*pi/lenz)

   if ax==0:
      k1=ky; k2=kz
   elif ax==1:
      k1=kx; k2=kz
   elif ax==2:
      k1=kx; k2=ky
  
   far = nf.fftshift(nf.fftn(ar))/(nx*ny*nz); fftea=0.5*np.abs(far)**2
   ffteb=fftea.sum(axis=ax)
   return k1,k2,ffteb

##
## 2D SPECTRUM OF A VECTOR
##
def SpecVec2D(arx,ary,arz,axx=2,lx=2*pi,ly=2*pi,lz=2*pi):
   """
      SpecVec2D(arx,ary,arz,axx=2,lx=2*pi,ly=2*pi,lz=2*pi)

      2D spectrum of a vector perpendicular to axis ax
   """
   k1,k2,ekx = Spec2D(arx,ax=axx,lenx=lx,leny=ly,lenz=lz)
   k1,k2,eky = Spec2D(ary,ax=axx,lenx=lx,leny=ly,lenz=lz)
   k1,k2,ekz = Spec2D(arz,ax=axx,lenx=lx,leny=ly,lenz=lz)
   return k1,k2,ekx,eky,ekz,ekx+eky+ekz

##
## DEF TO CALCULATE REDUCED SPECTRUM. 
##
def ReducedSpec(ar,ax=2,lenx=2*pi,leny=2*pi,lenz=2*pi):
   """
      ReducedSpec(ar,ax=2,lenx=2*pi,leny=2*pi,lenz=2*pi)
      
      Reduced spectrum of ar along axis ax
      
   """
   if len(ar) == 0:
      print('No array provided! Exiting!')
      return
   ar=ar-np.mean(ar)
   nx=np.shape(ar)[0];kx=nf.fftshift(nf.fftfreq(nx))*nx*(2*pi/lenx)
   ny=np.shape(ar)[1];ky=nf.fftshift(nf.fftfreq(ny))*ny*(2*pi/leny)
   nz=np.shape(ar)[2];kz=nf.fftshift(nf.fftfreq(nz))*nz*(2*pi/lenz)

   if ax==0:
      # Both are one because after the first sum, the new array as dim=2
      # First sum along y, then along z
      ax1=1; ax2=1
      kk=kx; nn=nx
   elif ax==1:
      # First sum along x (0), then along z (1 for the new configuration)
      ax1=0; ax2=1
      kk=ky; nn=ny
   elif ax==2:
      # First sum along x (0), then along y (0 for the new configuration)
      ax1=0; ax2=0
      kk=kz; nn=nz
  
   far = nf.fftshift(nf.fftn(ar))/(nx*ny*nz); fftea=0.5*np.abs(far)**2
   ffteb=fftea.sum(axis=ax1).sum(axis=ax2)
   dk = kk[1]-kk[0]
   return kk[nn/2:],ffteb[nn/2:]/dk

##
## REDUCED SPECTRUM OF A VECTOR
##
def ReducedSpecVec(arx,ary,arz,axx=2,lx=2*pi,ly=2*pi,lz=2*pi):
   """
      ReducedSpecVec(arx,ary,arz,axx=2,lx=2*pi,ly=2*pi,lz=2*pi)

      Reduced spectrum of a vector along axis ax
   """
   kwave,ekx = ReducedSpec(arx,ax=axx,lenx=lx,leny=ly,lenz=lz)
   kwave,eky = ReducedSpec(ary,ax=axx,lenx=lx,leny=ly,lenz=lz)
   kwave,ekz = ReducedSpec(arz,ax=axx,lenx=lx,leny=ly,lenz=lz)
   return kwave,ekx,eky,ekz,ekx+eky+ekz

def autocorrelation(ar,ax=0,dx=1.):
   nlen=np.shape(ar)[ax]/2
   r=np.zeros(nlen);corr=np.zeros(nlen)
   for i in range(nlen):
      ars=np.roll(ar,i,axis=ax)
      corr[i]=np.mean(ar*ars)
      r[i]=i*dx
   corr = corr/np.mean(ar**2)
   return r,corr

def autocvec(xx,yy,zz,ax=0,dx=1.):
   if ax == 0:
      nlen = np.shape(xx)[0]/2
   elif ax == 1:
      nlen = np.shape(xx)[1]/2
   elif ax == 2:
      nlen = np.shape(xx)[2]/2
   r=np.zeros(nlen);corr=np.zeros(nlen)
   tmp=np.zeros(np.shape(xx))
   for i in range(nlen):
      xxs=np.roll(xx,i,axis=ax)
      yys=np.roll(yy,i,axis=ax)
      zzs=np.roll(zz,i,axis=ax)
      corr[i]=np.mean(xx*xxs+yy*yys+zz*zzs)
      r[i]=i*dx
   corr = corr/np.mean(xx**2+yy**2+zz**2)
   return r,corr

def fcorr(ar1, ar2, ax=0,step=1, dx=1.):
   if ax == 0:
      rlen = np.shape(ar1)[0]/(2*step)
   elif ax == 1:
      rlen = np.shape(ar1)[1]/(2*step)
   elif ax == 2:
      rlen = np.shape(ar1)[2]/(2*step)
   r,cor=faf.correlation(ar1,ar2,ax,dx,step,rlen)
   return r,cor

def correlation(ar1,ar2,ax=0,dx=1.):
   nlen=np.shape(ar1)[ax]/2
   r=np.zeros(nlen);corr=np.zeros(nlen)
   for i in range(nlen):
      ars=np.roll(ar2,i,axis=ax)
      corr[i]=np.mean(ar1*ars)
      r[i]=i*dx
   corr = corr/np.mean(ar1*ar2)
   return r,corr

def calc_pdf(ar,min=99999,max=99999,weight=100,inc=0,ax=0,Normalized=False):
   if len(ar) == 0:
      print('No array provided! Exiting!')
      return
   if min == 99999:
      min=ar.min()
   if max == 99999:
      max=ar.max()
   # If PDF of increment, then increment the array
   if inc > 0:
      ar = ar - np.roll(ar,inc,axis=ax)
   # Find the total length of data set
   arsize=reduce(operator.mul, np.shape(ar),1)
   # Find the RMS value of data set and normalize to it.
   if Normalized == True:
      rmsval = np.sqrt(np.mean((ar-np.mean(ar))**2))
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
      binvals[i] = np.mean(arr[start:start+weight])
      pdf[i] = weight/(arr[start:start+weight].max()-arr[start:start+weight].min())
   pdf = pdf/arsize
   return binvals,pdf

def normhist(ar,min=99999,max=99999,nbins=100,inc=0,ax=0):
   if len(ar) == 0:
      print('No array provided! Exiting!')
      return
   # If PDF of increment, then increment the array
   if inc > 0:
      ar = ar - np.roll(ar,inc,axis=ax)
   # Find the total length of data set
   arsize=reduce(operator.mul, np.shape(ar),1)
   # Find the RMS value of data set and normalize to it.
   rmsval = np.sqrt(np.mean(ar**2))
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
      j = np.floor(arr[i]/dbin)
      normhist[j] = normhist[j]+1
   normhist = normhist/(arsize*dbin)
   return bins,normhist

def kurtosis(ar):
   return np.mean(ar**4)/np.mean(ar**2)**2

def windowff(nslices,kind=None):
   """
      Returns a window function to enforce quasi-preiodicity
      on an aperiodic signal. Window options are:
      BlackmanHarris, Hanning, Blackman, FlatTop, Welch, Tukey
   """
   import scipy.signal.windows as w
   if kind is None: 
      return np.ones(nslices)
   elif kind.lower() == 'tukey':
      return w.get_window((kind.lower(),0.1),nslices)
   else:
      return w.get_window(kind.lower(),nslices)

#  windowff=np.zeros(nslices)
   
#  for i in range(nslices):
#     tht = 2*pi*i/(nslices-1)
#     #  No Window
#     if kind is None:
#        windowff[i]=1.
#     #  Hanning
#     elif kind.lower() == "hanning":
#        windowff[i]=0.5*(1-np.cos(tht))
#     #  Blackman
#     elif kind.lower() == "blackmann":
#        windowff[i]=0.42659-0.49656*np.cos(tht)+0.076849*np.cos(2*tht)
#     #  BlackmanHarris
#     elif kind.lower() == "blackmannharris": 
#        windowff[i]=0.35875-0.48829*np.cos(tht)+0.14128*np.cos(2*tht)-0.01168*np.cos(3*tht)
#     #  Flat top Window
#     elif kind.lower() == "flattop": 
#        windowff[i] = 1. - 1.93*np.cos(tht) + 1.29*np.cos(2*tht) - 0.388*np.cos(3*tht) + 0.028*np.cos(4*tht)
#     #  Nuttall
#     elif kind.lower() == "nuttall": 
#        windowff[i] = 0.355768 - 0.487396*np.cos(tht) + 0.144232*np.cos(2*tht) - 0.012604*np.cos(3*tht)
#     #  Welch
#     elif kind.lower() == "welch": 
#        windowff[i] = 1. - ((i - (nslices-1)/2.)/((nslices+1)/2.))**2
#     # Tukey
#     elif kind.lower() == "tukey":
#        if i <= int(0.1*(nslices-1)/2):
#           windowff[i] = 0.5*(1+np.cos(tht/0.1 - pi))
#        if int(0.1*(nslices-1)/2) < i and i < int((nslices-1)*0.95):
#           windowff[i] = 1
#        if i >= int((nslices-1)*0.95):
#           windowff[i] = 0.5*(1+np.cos(tht/0.1 - 2*pi/0.1 +pi))
#  return windowff

#
# THE FOLLOWING ROUTINE CALCULATES THE SCALE DEPENDENT
# KURTOSIS OF A GIVEN ARRAY IN THE GIVEN DIRECTION
#
#def fsdk(ar,ax=0):
#   from CFLIB.lib.faf import kurtosis as kurt
#   nx=np.shape(ar)[ax]
#   ark=np.zeros(nx/2)
#   ddx=np.arange(nx/2)
#   for i in range(1,nx/2):
#      tmp=np.roll(ar,i,axis=ax)      
#      tmp=ar-tmp; rmsval=np.sqrt(np.mean(tmp**2))
#      if rmsval != 0:
#         tmp=tmp/rmsval 
#      ark[i] = kurt(tmp)
#   return ddx,ark
#def sdk(ar,ax=0):
#   nx=np.shape(ar)[ax]
#   ark=np.zeros(nx/2)
#   ddx=np.arange(nx/2)
#   for i in range(1,nx/2):
#      tmp=np.roll(ar,i,axis=ax)      
#      tmp=ar-tmp; rmsval=np.sqrt(np.mean(tmp**2))
#      if rmsval != 0:
#         tmp=tmp/rmsval 
#      ark[i] = kurtosis(tmp)
#   return ddx,ark
#
def sdk(ar,bs=1,fs=2,step=1,dx=1,ax=0,fort='y'):
   nx=np.shape(ar)[ax]
   ark=np.zeros((fs-bs+1)/step)
   ddx=np.zeros((fs-bs+1)/step)
   for i in range(bs,fs,step):
      idx = (i-bs)/step
      tmp=np.roll(ar,i,axis=ax)      
      tmp=ar-tmp; rmsval=np.sqrt(np.mean(tmp**2))
      if rmsval != 0:
         tmp=tmp/rmsval 
      if fort == 'y':
         ark[idx] = faf.kurtosis(tmp)
      else:
         ark[idx] = kurtosis(tmp)
      ddx[idx] = i
   return ddx*dx,ark
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
    except ValueError as msg:
        raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    order_range = list(range(order+1))
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
      theta: Angle of propagation in degrees
      kk: Wavenumber of interest in units of kdi
      
      Output is frequencies of the roots and the phase speeds w/k
      The roots are w[0]:Fast/Whistler, w[1]:Alfven/KAW, w[2]: Slow/Cyclotron
   """
   tht = theta*pi/180.
   ct = np.cos(tht)
   st = np.sin(tht)
   tt = st/ct
   cs=np.sqrt(beta/2.)*ca
   di =1.
   caksq=np.cos(tht)**2*ca**2
   cmsq=ca**2 + cs**2
   
   pcs=np.zeros(4)
   D = 1 + kk**2*de2
   # Find out the root of the quadratic dispersion relation
   pcs[0] = 1.
   pcs[1] = -(ca**2/D + cs**2 + caksq*(1+kk**2*di**2/D)/D)*kk**2
   pcs[2] = caksq*kk**4*(ca**2/D + cs**2 + cs**2*(1+kk**2*di**2/D))/D
   pcs[3] = -(cs**2*kk**6*caksq**2)/D**2
   w2 = np.roots(pcs); w = np.sqrt(w2)
   speeds= w/kk
   return w,speeds

def tfst(beta=0.6, ca=1., de2=0.000545, theta=0., kmin=1e-2, kmax=10., npoints=200,wrt='n'):
   """
      beta: Total plasma beta,
      ca: Alfven speed based on mean field, 
      de2: me/mi, 
      theta: Angle of propagation in degrees
      kkmin: Minimum Wavenumber of interest in units of kdi
      kkmax: Maximum wavenumber of interest in units of kdi

      Output is an array with 4 columns, k, w-fast, w-alf, w-slow
   """
   import matplotlib.pyplot as plt
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
      ofile=open('disp'+str(theta)+'.dat','w')
      print('#  k', 'Acoustic', 'Alfven', 'Fast', file=ofile)
      for i in range(npoints):
         print(kk[i],warray[2,i],warray[1,i],warray[0,i], file=ofile)
      ofile.close()

def tfev(beta=0.6, ca=1., de2=0.000545, theta=0., k=1.,aa=0.1):
   """
      beta: Total plasma beta,
      ca: Alfven speed based on mean field, 
      de2: me/mi, 
      theta: Angle of propagation in degrees
      k: wavenumber of interest in units of kdi

      Output: Prints the two fluid eigenvector to the screen
   """
   import numpy as np
   def amp(beta,de2,k,w,theta,aa):
      th=theta*pi/180.
      bb=1-w**2*(1+de2*k**2)/(k**2*np.cos(th)**2)
      sk='sin('+str(round(k,3))+'x)'
      ck='cos('+str(round(k,3))+'x)'
      def st(a):
         return str(round(a,3))
      return 'bx=0.'+\
            '\tby = '+st(round(2*aa,3))+ck+\
            '\tbz = '+st(-np.cos(th)*bb*2*aa/w)+sk+\
            '\nux = '+st(aa*k*bb*np.sin(2*th)/(w**2-beta*k**2))+sk+\
            '\tuy = '+st(-2*aa*k*np.cos(th)/w)+ck+\
            '\tuz = '+st(2*aa*k*bb*np.cos(th)**2/w**2)+sk+\
            '\t n = '+st((k*np.cos(th)/w)*aa*k*bb*np.sin(2*th)/(w**2-beta*k**2))+sk
   f,s=tfps(beta,ca,de2,theta,k)
#     The roots are w[0]:Fast/Whistler, w[1]:Alfven/KAW, w[2]: Slow/Cyclotron
   print('Fast/ Whistler')
   print(amp(beta,de2,k,f[0],theta,aa))
   print('##################')
   print()
   print('Alfven/ KAW')
   print(amp(beta,de2,k,f[1],theta,aa))
   print('##################')
   print()
   print('Slow/ Cyclotron')
   print(amp(beta,de2,k,f[2],theta,aa))
   

##
## FUNCTION TO FIT A POWERLAW BETWEEN TWO SERIES (X,Y) BETWEEN
## TWO POINTS XI, XF
##
def fitpowerlaw(ax,ay,xi,xf):
   idxi=np.argmin(abs(ax-xi))
   idxf=np.argmin(abs(ax-xf))
   xx=np.linspace(xi,xf,100)
   z=np.polyfit(np.log(ax[idxi:idxf]),np.log(ay[idxi:idxf]),1); 
   p=np.poly1d(z);
   pwrl=np.exp(p(np.log(xx)))
   return z,xx,pwrl
##
## FUNCTION TO FIT A POLYNOMIAL OF ORDER N BETWEEN TWO SERIES 
## (X,Y) BETWEEN TWO POINTS XI, XF
##
def fitpoly(ax,ay,n,xi,xf,numpts=100):
   idxi=np.argmin(abs(ax-xi))
   idxf=np.argmin(abs(ax-xf))
   xx=np.linspace(xi,xf,numpts)
   z=np.polyfit(ax[idxi:idxf],ay[idxi:idxf],n); 
   p=np.poly1d(z);
   poly=p(xx)
   return z,xx,poly
##
## Plot a powerlaw
##

def pltpwrl(x0,y0,xi=1,xf=10,alpha=-1.66667,ax=None,**kwargs):
   """
      Plots a power-law with exponent alpha between the 
      xrange (xi,xf) such that it passes through x0,y0
   """
   import numpy as np
   import matplotlib.pyplot as plt
   x=np.linspace(xi,xf,50)
   if ax is None:
      ax=plt.gca()
   ax.plot(x,(y0*x0**-alpha)*x**alpha,**kwargs)

##
## First derivative for periodic arrays
##
def pderiv(ar,dx=1.,ax=0,order=2,smth=None):
   """
      pderiv gives the first partial derivative
      of a periodic array along a given axis.

      Inputs:
         ar - The input array
         dx - Grid spacing, defaults to 1.
         ax - Axis along which to take the derivative
         order - Order of accuracy, (1,2) defaults to 2

      Output:
         dar - The derivative array
   """
   if smth is not None:
      ar = gf(ar,sigma=smth)
   if order == 1:
      dar = (np.roll(ar,-1,axis=ax)-ar)/dx
   elif order == 2:
      dar = (np.roll(ar,-1,axis=ax)-np.roll(ar,1,axis=ax))/(2*dx)
   
   return dar 

##
## Second derivative for periodic arrays
##
def pdderiv(ar,dx=1.,ax=0,order=4,smth=None):
   """
      pderiv gives the double partial derivative
      of a periodic array along a given axis.

      Inputs:
         ar - The input array
         dx - Grid spacing, defaults to 1.
         ax - Axis along which to take the derivative
         order - Order of accuracy, (2,4) defaults to 2

      Output:
         dar - The derivative array
   """
   if smth is not None:
      ar=gf(ar,sigma=smth)
   if order == 2:
      dar = (np.roll(ar,-1,axis=ax) - 2*ar + np.roll(ar,1,axis=ax))/dx**2
   elif order == 4:
      dar = (-np.roll(ar,-2,axis=ax) + 16*np.roll(ar,-1,axis=ax) - 30*ar + 16*np.roll(ar,1,axis=ax)-np.roll(ar,2,axis=ax))/(12*dx**2)

   return dar

##
## Divergence of a periodic array
## 
def pdiv(arx,ary,arz,dx=1,dy=1,dz=1,smth=None):
   return pderiv(arx,dx=dx,ax=0,smth=smth)+pderiv(ary,dx=dy,ax=1,smth=smth)+pderiv(arz,dx=dz,ax=2,smth=smth)

##
## Curl of a periodic array
##
def pcurl(arx,ary,arz,dx=1.,dy=1.,dz=1.,smth=None):
   return pderiv(arz,dx=dy,ax=1,smth=smth)-pderiv(ary,dx=dz,ax=2,smth=smth),\
          pderiv(arx,dx=dz,ax=2,smth=smth)-pderiv(arz,dx=dx,ax=0,smth=smth),\
          pderiv(ary,dx=dx,ax=0,smth=smth)-pderiv(arx,dx=dy,ax=1,smth=smth)

def pcurlx(ary,arz,dx=1.,dy=1.,dz=1.,smth=None):
   return pderiv(arz,dx=dy,ax=1,smth=smth)-pderiv(ary,dx=dz,ax=2,smth=smth)

def pcurly(arz,arx,dx=1.,dy=1.,dz=1.,smth=None):
   return pderiv(arx,dx=dz,ax=2,smth=smth)-pderiv(arz,dx=dx,ax=0,smth=smth)

def pcurlz(arx,ary,dx=1.,dy=1.,dz=1.,smth=None):
   return pderiv(ary,dx=dx,ax=0,smth=smth)-pderiv(arx,dx=dy,ax=1,smth=smth)

##
## Gradient of a periodic scalar array
##
def pgrad(ar,dx=1,dy=1,dz=1,smth=None):
   return pderiv(ar,dx=dx,ax=0,smth=smth),\
          pderiv(ar,dx=dy,ax=1,smth=smth),\
          pderiv(ar,dx=dz,ax=2,smth=smth)

##
## Conditional PDFs on signed variable
##
def scpdf(d,a=None,b=None,mag=1.,sav='n'):
   import matplotlib.pyplot as plt
   ms = str(mag)+'$\sigma$'
   rmsa = np.sqrt(np.mean(d[a]**2))
   rmsb = np.sqrt(np.mean(d[b]**2)); mag=mag*rmsb; 
   bn,pn=calc_pdf(d[a][np.where( d[b] < -mag)])
   bs,ps=calc_pdf(d[a][np.where((d[b] < mag) & (d[b] > -mag))])
   bp,pp=calc_pdf(d[a][np.where( d[b] > mag)])
   plt.clf()
   plt.plot(bn*rmsa,pn,label=b+' < -'+ms)
   plt.plot(bs*rmsa,ps,label='|'+b+'| < '+ms)
   plt.plot(bp*rmsa,pp,label=b+' >  '+ms)
   plt.xlabel(a)
   plt.yscale('log')
   plt.legend(loc='best')
   if sav == 'y':
      plt.savefig(a+'-conditioned-on-'+b+'-signed.png',dpi=100)
   return {'bn':bn,'pn':pn,'bs':bs,'ps':ps,'bp':bp,'pp':pp, 'rmsbn':rmsa}

##
## Conditional PDFs on unsigned variable
##
def ucpdf(d,a=None,b=None,mag=1.,sav='n'):
   import matplotlib.pyplot as plt
   ms  = str(mag)+'$\sigma$'
   ms2 = str(2*mag)+'$\sigma$'
   rmsa = np.sqrt(np.mean(d[a]**2))
   rmsb = np.sqrt(np.mean(d[b]**2)); mag=mag*rmsb; 
   bn,pn=calc_pdf(d[a][np.where( abs(d[b]) < mag)])
   bs,ps=calc_pdf(d[a][np.where((abs(d[b]) > mag) & (abs(d[b]) < 2*mag))])
   bp,pp=calc_pdf(d[a][np.where( abs(d[b]) > 2*mag)])
   plt.clf()
   plt.plot(bn*rmsa,pn,label='|'+b+'| < '+ms)
   plt.plot(bs*rmsa,ps,label=ms+'< |'+b+'| < '+ms2)
   plt.plot(bp*rmsa,pp,label='|'+b+'| > '+ms2)
   plt.xlabel(a)
   plt.yscale('log')
   plt.legend(loc='best')
   if sav == 'y':
      plt.savefig(a+'-conditioned-on-'+b+'-usigned.png',dpi=100)
   return {'bn':bn,'pn':pn,'bs':bs,'ps':ps,'bp':bp,'pp':pp, 'rmsbn':rmsa}

