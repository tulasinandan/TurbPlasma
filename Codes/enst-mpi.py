#!/usr/bin/env python
#
#
import os
import sys
import numpy as np
from mpi4py import MPI
sys.path.insert(0,os.environ['HOME']+'/AJGAR/TurbPlasma/')
from p3d import p3d
#from subs import create_object

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()
status = MPI.Status()
#
# INITIALIZATION
rc=p3d('.','000')
rc.vars2load(['ensti','enste','pali','pale','jx','jy','jz'])
nm=rc.lambdae/rc.dx
nt=rc.numslices-rc.numslices%size
#
# END INITIALIZATION

# RANK 0 (MASTER CORE) SETS UP SOME STUFF FOR ITSELF
if rank == 0:
   t_start=MPI.Wtime()
   tt =np.zeros(nt)
   j2 =np.zeros(nt)
   wi2=np.zeros(nt)
   we2=np.zeros(nt)
   Pi =np.zeros(nt)
   Pe =np.zeros(nt)

comm.Barrier()

for it in range(rank,nt,comm.size):
   rc.loadslice(it,smth=nm)
   ltt =rc.time 
   lj2 =np.mean(rc.jx**2+rc.jy**2+rc.jz**2)
   lwi2=np.mean(rc.omix**2+rc.omiy**2+rc.omiz**2)
   lwe2=np.mean(rc.omex**2+rc.omey**2+rc.omez**2)
   lPi =np.mean(rc.pali)
   lPe =np.mean(rc.pale)
## PRINT TO STDOUT FOR RECAPTURE IF THE RUN CRASHES
   print(ltt,ltt*2*np.pi/rc.lx,lj2,lwi2,lwe2,lPi,lPe)

   if rank == 0:
      tt[ it] = ltt
      j2[ it] = lj2 
      wi2[it] = lwi2
      we2[it] = lwe2
      Pi[ it] = lPi 
      Pe[ it] = lPe 
   
   if rank > 0:
      snddata=[it,ltt,lj2,lwi2,lwe2,lPi,lPe]
      comm.send(snddata, dest=0, tag=13)
   else:
#  if rank == 0:
      for src in range(1,comm.size):
         rcvdata=comm.recv(source=src,tag=13,status=status)
         tt[ rcvdata[0]] = rcvdata[1]
         j2[ rcvdata[0]] = rcvdata[2]
         wi2[rcvdata[0]] = rcvdata[3]
         we2[rcvdata[0]] = rcvdata[4]
         Pi[ rcvdata[0]] = rcvdata[5]
         Pe[ rcvdata[0]] = rcvdata[6]
   print('Done time slice ',it,' on proc ',rank)
comm.Barrier()

if rank == 0:
   ofl=open('enst.'+rc.dirname+'.dat','w')
   print('# Time\t tnl\t j^2\t wi^2\t we^2\t Pi\t Pe', file=ofl)
   for i in range(nt):
      print(tt[i], tt[i]*2*np.pi/rc.lx, j2[i], wi2[i], we2[i], Pi[i], Pe[i], file=ofl)
   t_fin = MPI.Wtime()-t_start
   print('Total time taken %0.3f'%t_fin)

