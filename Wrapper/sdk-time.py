#!/usr/bin/env python
import AnalysisFunctions as af
import numpy as np
from subs import create_object,ask_for_steps

print(" NEED TO SIT DOWN AND FIX IT SUCH THAT ")
print(" IT CALCULATES THE INCREMENT ARRAY AND ")
print(" KURTOSIS ACCORDING TO THE DIRECTION   ")
print(" THAT YOU PROVIDE.                     ")
rc=create_object()
rc.vars2load(['bx','by','bz'])
bs,fs,step=ask_for_steps(rc.numslices)
nt=(fs-bs)/step
bxk=np.zeros(rc.nx/2)
byk=np.zeros(rc.nx/2)
bzk=np.zeros(rc.nx/2)
xk = np.arange(max([rc.nx,rc.ny,rc.nz])/2)*rc.dx

for it in range(bs,fs,step):
   print('time slice',it)
   rc.loadslice(it)
   bxk[:rc.ny/2] = af.sdk(rc.bx,ax=1)[1]
   byk[:rc.nx/2] = af.sdk(rc.by,ax=0)[1]
   bzk[:rc.nx/2] = af.sdk(rc.bz,ax=0)[1]
   outf=open('sdk-'+rc.dirname+'-%03d.dat'%it,'w')
   for i in range(rc.nx/2):
      print(xk[i],bxk[i],byk[i],bzk[i], file=outf) 
   outf.close()

rc.fin()
   
