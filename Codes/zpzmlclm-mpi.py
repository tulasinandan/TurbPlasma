#!/usr/bin/env python
#
#  Calculating the Z+ & Z- values based on center of mass 
#  velocity instead of proton velocity.
#
import os
import sys
sys.path.insert(0,os.environ['HOME']+'/P3D-PLASMA-PIC/p3dpy/')
from subs import create_object
import numpy as np
from mpi4py import MPI
import AnalysisFunctions as af

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()
status = MPI.Status()

rc=create_object()
rc.vars2load(['bx','by','bz','jix','jiy','jiz','ni','jex','jey','jez','ne'])
comm.Barrier()
nt=rc.numslices-rc.numslices%size
ie=1/np.e
if rank == 0:
   print('NUMBER OF SLICES ',nt)
   t_start=MPI.Wtime()
   tt =np.zeros(nt)
   ezp=np.zeros(nt) 
   ezm=np.zeros(nt)
   lcp=np.zeros(nt)
   lcm=np.zeros(nt)

#for it in range(bs,fs,step):
for it in range(rank,nt,comm.size):
   rc.loadslice(it)
 
   rho=rc.ni+rc.ne*rc.m_e
   rc.jix=rc.jix/rc.ni; rc.jiy=rc.jiy/rc.ni;rc.jiz=rc.jiz/rc.ni
   rc.jex=rc.jex/rc.ne; rc.jey=rc.jey/rc.ne;rc.jez=rc.jez/rc.ne
# 
   vcx = (rc.jix+rc.jex*rc.m_e)/(1+rc.m_e)
   vcy = (rc.jiy+rc.jey*rc.m_e)/(1+rc.m_e)
   vcz = (rc.jiz+rc.jez*rc.m_e)/(1+rc.m_e)
   rc.bx=(rc.bx-np.mean(rc.bx))/np.sqrt(rho)
   rc.by=(rc.by-np.mean(rc.by))/np.sqrt(rho)
   rc.bz=(rc.bz-np.mean(rc.bz))/np.sqrt(rho)
   vcx=vcx-np.mean(vcx);vcy=vcy-np.mean(vcy);vcz=vcz-np.mean(vcz)
#
#
   zpx=vcx+rc.bx; zpy=vcy+rc.by; zpz=vcz+rc.bz
#  
   lezp=0.5*np.mean(rho*(np.abs(zpx)**2+np.abs(zpy)**2+np.abs(zpz)**2))
   zmx=vcx-rc.bx; zmy=vcy-rc.by; zmz=vcz-rc.bz
# 
   lezm=0.5*np.mean(rho*(np.abs(zmx)**2+np.abs(zmy)**2+np.abs(zmz)**2))
## #Correlation for z+
  #kwave,zkx,zky,zkz,zk = af.PerpSpecVec(zpx,zpy,zpz)
   kwave,zkx,zky,zkz,zk = af.fperpspecvec(zpx,zpy,zpz,2,rc.nx/2)
   llcp = (np.pi/sum(zk[1:]))*sum(zk[1:]/kwave[1:])
## #Correlation for z-
  #kwave,zkx,zky,zkz,zk = af.PerpSpecVec(zmx,zmy,zmz)
   kwave,zkx,zky,zkz,zk = af.fperpspecvec(zmx,zmy,zmz,2,rc.nx/2)
   llcm = (np.pi/sum(zk[1:]))*sum(zk[1:]/kwave[1:])
# 
   ltt = round(it*rc.dtmovie,4)

## PRINT TO STDOUT FOR RECAPTURE IF THE RUN CRASHES
   print(ltt,ltt*2*np.pi/rc.lx,lezp,lezm,llcp,llcm)

   if rank == 0:
      tt[ it] = ltt
      ezp[it] = lezp
      ezm[it] = lezm
      lcp[it] = llcp
      lcm[it] = llcm
   
   if rank > 0:
      snddata=[it,ltt,lezp,lezm,llcp,llcm]
      comm.send(snddata, dest=0, tag=13)
   else:
#  if rank == 0:
      for src in range(1,comm.size):
         rcvdata=comm.recv(source=src,tag=13,status=status)
         tt[ rcvdata[0]] = rcvdata[1]
         ezp[rcvdata[0]] = rcvdata[2]
         ezm[rcvdata[0]] = rcvdata[3]
         lcp[rcvdata[0]] = rcvdata[4]
         lcm[rcvdata[0]] = rcvdata[5]
   print('Done time slice ',it,' on proc ',rank)
comm.Barrier()

if rank == 0:
   print('Done Computing. Writing the file now')
   outf=open('zpzmlc.'+rc.dirname+'.dat','w')
   print('#','t','tt','ezp','ezm','lcp','lcm', file=outf)
   for i in range(nt):
      print(tt[i],tt[i]*2*np.pi/rc.lx,ezp[i],ezm[i],lcp[i],lcm[i], file=outf)
   outf.close()
   t_fin = MPI.Wtime()-t_start
   print('Total time taken %0.3f'%t_fin)

rc.fin()

