from numpy.ctypeslib import ndpointer
import matplotlib.pyplot as plt
import numpy as np
import ctypes

def correlation(ar1,ar2,ax=0,dx=1.):
   nlen=np.shape(ar1)[ax]/2
   r=np.zeros(nlen);corr=np.zeros(nlen)
   for i in range(nlen):
      ars=np.roll(ar2,i,axis=ax)
      corr[i]=np.mean(ar1*ars)
      r[i]=i*dx
   corr = corr/np.mean(ar1*ar2)
   return r,corr

def c_correlation(ar1,ar2,ax=0,dx=1.):
   lib = ctypes.cdll.LoadLibrary('/home/tulasi/P3D-PLASMA-PIC/p3dpy/helloworld.so')
   func = lib.c_correlation
   func.restype = None
   func.argtypes = [ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),   #ar1
                    ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),   #ar2
                    ctypes.c_double,              #dx
                    ctypes.c_int,                 #nlen
                    ctypes.c_int,                 #nx
                    ctypes.c_int,                 #ny
                    ctypes.c_int,                 #nz
                    ctypes.c_int,                 #ax
                    ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),   #r
                    ndpointer(ctypes.c_double, flags="C_CONTIGUOUS")]   #corr
# nlen finds the length of the array in the specified direction
   nlen=np.shape(ar2)[ax]/2; 
   nx=np.shape(ar1)[0]; 
   ny=np.shape(ar1)[1]; 
   nz=np.shape(ar1)[2]
   r=np.zeros(nlen);corr=np.zeros(nlen)
   func(np.ascontiguousarray(ar1),np.ascontiguousarray(ar2),dx,nlen,nx,ny,nz,ax,r,corr)
#  func(ar1,ar2,dx,nlen,nx,ny,nz,ax,r,corr)
   return r,corr

lib = ctypes.cdll.LoadLibrary('/home/tulasi/P3D-PLASMA-PIC/p3dpy/helloworld.so')
cmean=lib.mean
cmean.restype = ctypes.c_double
cmean.argtypes = [ndpointer(ctypes.c_double), ctypes.c_int, ctypes.c_int, ctypes.c_int]

ckurt=lib.kurtosis
ckurt.restype = ctypes.c_double
ckurt.argtypes = [ndpointer(ctypes.c_double), ctypes.c_int, ctypes.c_int, ctypes.c_int]

from subs import *
from scipy.ndimage import gaussian_filter as gf
rdir,rc=create_p3do_object()
rc.vars2load(['bx','pxx','jtotz'])
rc.loadslice(30)
px=gf(rc.pxx,sigma=2)
jz=gf(rc.jtotz,sigma=2)
bx=gf(rc.bx,sigma=2)
ar1=px
ar2=px
r1,corr1=correlation(ar1,ar2)
r2,corr2=c_correlation(ar1,ar2)
print np.mean(px)
print cmean(px, rc.nx, rc.ny, rc.nz)
plt.plot(r1,corr1)
plt.plot(r2,corr2/corr2[0],'x')
plt.show()

def csdk(ar,ax=0):
   csdk=lib.c_sdk
   csdk.restype = None
   csdk.argtypes = [ndpointer(ctypes.c_double), ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_int,ndpointer(ctypes.c_double),ndpointer(ctypes.c_double)]
   nx=np.shape(ar)[0]; ny=np.shape(ar)[1]; nz=np.shape(ar)[2]
   ddx=np.zeros(np.shape(ar)[ax]/2)
   ark=np.zeros(np.shape(ar)[ax]/2)
   csdk(ar,ax,nx,ny,nz,ddx,ark)
   return ddx,ark

pmt=lib.printmat
pmt.restype=None
pmt.argtypes = [ndpointer(ctypes.c_double), ctypes.c_int, ctypes.c_int]
a=np.array([[0.,0.,1.],[1.,4.,6.],[3.,9.,0.]])
pmt(a,3,3)
#   for j in range(64):
#      bx[i,j,0]=np.sin(x[i])*np.cos(y[j])+1.203215
#r1,corr1=correlation(bx,bx)
#r2,corr2=c_correlation(bx,bx)
