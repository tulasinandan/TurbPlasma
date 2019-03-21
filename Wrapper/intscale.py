#!/usr/bin/env python
import numpy as np
from subs import create_object,ask_for_steps
from AnalysisFunctions import fcorr 
rc=create_object()
rc.vars2load(['bx','by'])
ofile=open('intscale.'+rc.dirname+'.dat','w')
ie=1/np.e
bs,fs,step=ask_for_steps(rc.numslices)
for i in range(bs,fs,step):
   print(i)
   rc.loadslice(i)
   rx,bxcor=fcorr(rc.bx,rc.bx,ax=0,dx=rc.dx)
#  lxc=rx[abs(bxcor-ie).argmin()]
   lxint=np.sum(bxcor)*rc.dx
   ry,bycor=fcorr(rc.by,rc.by,ax=1,dx=rc.dy)
#  lyc=ry[abs(bycor-ie).argmin()]
   lyint=np.sum(bycor)*rc.dy

   print(rc.time,lxint,lyint,0.5*(lxint+lyint), file=ofile)
ofile.close()
