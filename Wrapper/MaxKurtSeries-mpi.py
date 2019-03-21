#!/usr/bin/env python
import os
import sys
sys.path.insert(0,os.environ['HOME']+'/P3D-PLASMA-PIC/p3dpy/')
import numpy as np
from mpi4py import MPI
from subs import create_object
import AnalysisFunctions as af 

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
rc.vars2load(['bx','by'])
rc.loadenergies()

#
# DECIDE THE NUMBER OF SLICES
#
extslc=rc.numslices%size
nt=(rc.numslices-extslc)/size
if rank==0:
   bs=0; fs=nt+extslc
else:
   bs=extslc+rank*nt; fs=bs+nt

print(rank, bs, fs, nt)
print() 
print()
#
# CREATE OUTPUT ARRAYS
#
if rank==0:
   sclmx=np.zeros(rc.numslices); krtmx=np.zeros(rc.numslices); tsris=np.zeros(rc.numslices)
else:
   snddata=np.zeros((3,nt))

#
# MAIN LOOP
#
bl=0.2; fl=5.; stepl=2
if fl > 0.5*rc.lx:
   fl=0.5*rc.lx
if bl < rc.dx:
   bl=rc.dx
bsl=int(bl/rc.dx); fsl=int(fl/rc.dx)

if rank == 0:
   t_start=MPI.Wtime()
comm.Barrier()
for i in range(bs,fs):
   rc.loadslice(i)
   idx = (bs-i)
   #Assign the appropriate nonlinear time to the timeseries
   if rank == 0:
      tsris[i] = rc.ta[np.argmin(np.abs(rc.t-rc.time))]
   else:
      snddata[0,idx] = rc.ta[np.argmin(np.abs(rc.t-rc.time))]

   #Compute the scale dependent kurtosis between bl,fl
   rx,sdx=af.sdk(rc.bx,bs=bsl,fs=fsl,step=stepl,dx=rc.dx,ax=0)
   ry,sdy=af.sdk(rc.by,bs=bsl,fs=fsl,step=stepl,dx=rc.dx,ax=1)

   #Find the location of maximum value of kurtosis between sdx,sdy
   tsdk=np.array([sdx,sdy])
   rr =np.array([rx,ry])
   xx=np.unravel_index(tsdk.argmax(),tsdk.shape)

   if rank == 0:
      sclmx[i]=rr[xx]
      krtmx[i]=tsdk[xx]
      print('rank,t,sclmx,mxkurt\t',rank,tsris[i],sclmx[i],krtmx[i])
   else:
      snddata[1,idx]=rr[xx]
      snddata[2,idx]=tsdk[xx]
      print('rank,t,sclmx,mxkurt\t',rank,snddata[0,idx],snddata[1,idx],snddata[2,idx])


#
# PROC 0 COLLECTS DATA AND WRITES THE FILE
#
if rank > 0:
   comm.send(snddata, dest=0, tag=13)
else:
   for src in range(1,comm.size):
      tbs=extslc+src*nt; tfs=extslc+(src+1)*nt
      rcvdata=comm.recv(source=src,tag=13,status=status)
      for j in range(tbs,tfs):
         tsris[j] =rcvdata[0,j-tbs]
         sclmx[j]=rcvdata[1,j-tbs]
         krtmx[j]=rcvdata[2,j-tbs]
#
comm.Barrier()
#
if rank == 0:
   print('Done Computing. Writing the file now')
   outf=open('mxkurtdi.'+rc.dirname+'.dat','w')
   print('#','t,\t tt,\t lxc,\t lyc,\t lc', file=outf)
   for i in range(rc.numslices):
      print(tsris[i],sclmx[i],krtmx[i], file=outf)
   outf.close()
   t_fin = MPI.Wtime()-t_start
   print('Total time taken %0.3f'%t_fin)
#
