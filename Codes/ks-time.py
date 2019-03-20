#!/usr/bin/env python
import AnalysisFunctions as af
import numpy as np
from subs import create_object, ask_for_steps

rc=create_object()
rc.vars2load(['bx','by','bz'])
bs,fs,step=ask_for_steps(rc.numslices)
nt=(fs-bs)/step

for it in range(bs,fs,step):
   print('time slice',it)
   rc.loadslice(it)
   kk,ekbx,ekby,ekbz,ekb=af.fperpspecvec(rc.bx,rc.by,rc.bz)
   outf=open('perpspec-'+rc.dirname+'-%03d.dat'%it,'w')
   for i in range(len(kk)):
      print(kk[i],ekbx[i],ekby[i],ekbz[i],ekb[i], file=outf)
   outf.close()

rc.fin()
   
