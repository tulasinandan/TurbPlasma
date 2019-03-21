#!/usr/bin/env python
from subs import create_object, ask_for_steps
import AnalysisFunctions as af
import numpy.fft as nf
import numpy as np

rc=create_object()
rc.vars2load(['bx','by','bz','jix','jiy','jiz','jex','jey','jez','ni','ne'])
t  = rc.dtmovie*np.arange(rc.numslices)
tt = rc.dtmovie*np.arange(rc.numslices)*2*np.pi/rc.lx
bs,fs,step=ask_for_steps(rc.numslices)

ii=0+1j; kmin=2*np.pi/rc.lx;kminsq=kmin**2

kx,ky,kz,km=af.create_kgrid(rc.nx,rc.ny,rc.nz,lx=rc.lx,ly=rc.ly,lz=rc.lz)
ksq=km**2
ax=np.zeros((rc.nx,rc.ny,rc.nz),dtype=complex)
ay=np.zeros((rc.nx,rc.ny,rc.nz),dtype=complex)
az=np.zeros((rc.nx,rc.ny,rc.nz),dtype=complex)

wx=np.zeros((rc.nx,rc.ny,rc.nz),dtype=complex)
wy=np.zeros((rc.nx,rc.ny,rc.nz),dtype=complex)
wz=np.zeros((rc.nx,rc.ny,rc.nz),dtype=complex)

outf=open('AW2.'+rc.dirname+'.dat','w')
print('t[it],\t tt[it],\t eb,\t az2,\t ev,\t wz2,\t eb+ev-az2-wz2', file=outf)
for it in range(bs,fs,step):
   print('time slice',it)
   rc.loadslice(it)
   # FOURIER TRANSFORM THE MAGNETIC FIELD FLUCTUATIONS
   bx=rc.bx-np.mean(rc.bx); fbx=nf.fftshift(nf.fftn(bx))/(rc.nx*rc.ny*rc.nz)
   by=rc.by-np.mean(rc.by); fby=nf.fftshift(nf.fftn(by))/(rc.nx*rc.ny*rc.nz)
   bz=rc.bz-np.mean(rc.bz); fbz=nf.fftshift(nf.fftn(bz))/(rc.nx*rc.ny*rc.nz)
   # COMPUTE THE CENTER OF MASS VELOCITY
   rc.jix=rc.jix/rc.ni; rc.jiy=rc.jiy/rc.ni;rc.jiz=rc.jiz/rc.ni
   rc.jex=rc.jex/rc.ne; rc.jey=rc.jey/rc.ne;rc.jez=rc.jez/rc.ne
   vcx = (rc.jix+rc.jex*rc.m_e)/(1+rc.m_e)
   vcy = (rc.jiy+rc.jey*rc.m_e)/(1+rc.m_e)
   vcz = (rc.jiz+rc.jez*rc.m_e)/(1+rc.m_e)
   # FOURIER TRANSFORM THE CENTER OF MASS VELOCITY
   vcx=vcx-np.mean(vcx); fvx=nf.fftshift(nf.fftn(vcx))/(rc.nx*rc.ny*rc.nz)
   vcy=vcy-np.mean(vcy); fvy=nf.fftshift(nf.fftn(vcy))/(rc.nx*rc.ny*rc.nz)
   vcz=vcz-np.mean(vcz); fvz=nf.fftshift(nf.fftn(vcz))/(rc.nx*rc.ny*rc.nz)
   # COMPUTE THE VORTICITY AND VECTOR POTENTIAL
   for i in range(rc.nx):
      for j in range(rc.ny):
         for k in range(rc.nz):
            if ksq[i,j,k] != 0.:
               # VECTOR POTENTIAL
               ax[i,j,k] = ii*(ky[j]*fbz[i,j,k]-kz[k]*fby[i,j,k])/ksq[i,j,k]
               ay[i,j,k] = ii*(kz[k]*fbx[i,j,k]-kx[i]*fbz[i,j,k])/ksq[i,j,k]
               az[i,j,k] = ii*(kx[i]*fby[i,j,k]-ky[j]*fbx[i,j,k])/ksq[i,j,k]
               # VORTICITY
               wx[i,j,k] = ii*(ky[j]*fvz[i,j,k]-kz[k]*fvy[i,j,k])/ksq[i,j,k]
               wy[i,j,k] = ii*(kz[k]*fvx[i,j,k]-kx[i]*fvz[i,j,k])/ksq[i,j,k]
               wz[i,j,k] = ii*(kx[i]*fvy[i,j,k]-ky[j]*fvx[i,j,k])/ksq[i,j,k]

   asq = abs(ax)**2 + abs(ay)**2 + abs(az)**2; asq=np.sum(asq,axis=2)
   wsq = abs(wx)**2 + abs(wy)**2 + abs(wz)**2; wsq=np.sum(wsq,axis=2)

   eb=np.sum(abs(fbx**2+fby**2+fbz**2))*0.5
   ev=np.sum(abs(fvx**2+fvy**2+fvz**2))*0.5
   az2=kminsq*np.sum(asq)*0.5
   wz2=np.sum(wsq)*0.5
   
   print(t[it], tt[it], eb, az2, ev, wz2, eb+ev-az2-wz2, file=outf)

outf.close()
rc.fin()
