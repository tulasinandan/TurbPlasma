#!/usr/bin/env python
import os
import sys
sys.path.insert(0,os.environ['HOME']+'/P3D-PLASMA-PIC/p3dpy/')
import numpy as np
from mpi4py import MPI
from subs import create_object
from AnalysisFunctions import fcorr 

#
# MPI INITIALIZATION
#
comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()
status = MPI.Status()

#
# CREATE P3D OBJECT & OPEN OUTPUT FILE
#
rc=create_object()
rc.vars2load(['bx','by','bz'])
ie=1/np.e

kind = 'int'

#
# DECIDE THE NUMBER OF SLICES
#
extslc=rc.numslices%size
nt=(rc.numslices-extslc)/size
if rank==0:
   bs=0; fs=nt+extslc
else:
   bs=extslc+rank*nt; fs=bs+nt
#
# CREATE OUTPUT ARRAYS
#
if rank==0:
   lxc=np.zeros(rc.numslices);lyz=np.zeros(rc.numslices);lc=np.zeros(rc.numslices);tt=np.zeros(rc.numslices)
else:
   snddata=np.zeros((4,nt))
#
# MAIN LOOP
#
if rank == 0:
   t_start=MPI.Wtime()
comm.Barrier()
for i in range(bs,fs):
   rc.loadslice(i)
   rx,bxcor=fcorr(rc.bx,rc.bx,ax=0,dx=rc.dx)
   ry,bycor=fcorr(rc.by,rc.by,ax=1,dx=rc.dy)

   snddata[0,i]=rc.time
   if kind == "ie":
      snddata[1,i]=rx[abs(bxcor-ie).argmin()]
      snddata[2,1]=ry[abs(bycor-ie).argmin()]
   elif kind == "int":
      snddata[1,i]=np.sum(bxcor)*rc.dx
      snddata[2,1]=np.sum(bycor)*rc.dy
   snddata[3,i] = 0.5*(snddata[1,i]+snddata[2,i])

   print('t,lxc,lyc,lc\t',i,snddata[0,i],snddata[1,i],snddata[2,i],snddata[3,i])

#
# PROC 0 COLLECTS DATA AND WRITES THE FILE
#
if rank > 0:
#  snddata=[tt,lxc,lyc,lc]
   comm.send(snddata, dest=0, tag=13)
else:
   for src in range(1,comm.size):
      tbs=extslc+src*nt; tfs=extslc+(src+1)*nt
      rcvdata=comm.recv(source=src,tag=13,status=status)
      for j in range(tbs,tfs):
         tt[j] =rcvdata[0,j-tbs]
         lxc[j]=rcvdata[1,j-tbs]
         lyc[j]=rcvdata[2,j-tbs]
         lc[j] =rcvdata[3,j-tbs]
#
comm.Barrier()
#
if rank == 0:
   print('Done Computing. Writing the file now')
   outf=open('acorl.'+rc.dirname+'.dat','w')
   print('#','t,\t tt,\t lxc,\t lyc,\t lc', file=outf)
   for i in range(rc.numslices):
      print(tt[i],tt[i]*2*np.pi/rc.lx,lxc[i],lyc[i],lc[i], file=outf)
   outf.close()
   t_fin = MPI.Wtime()-t_start
   print('Total time taken %0.3f'%t_fin)
#
