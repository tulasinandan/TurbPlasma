#!/usr/bin/env python
import AnalysisFunctions as af
import numpy as np
from subs import create_object,ask_for_steps

rc=create_object()
rc.vars2load(['bx','by','bz'])
bs,fs,step=ask_for_steps(rc.numslices)
nt=(fs-bs)/step

for it in range(bs,fs,step):
   print('time slice',it)
   rc.loadslice(it)
   kx,xkbx,xkby,xkbz,xkb=af.ReducedSpecVec(rc.bx,rc.by,rc.bz,axx=0,lx=rc.lx)
   ky,ykbx,ykby,ykbz,ykb=af.ReducedSpecVec(rc.bx,rc.by,rc.bz,axx=1,ly=rc.ly)
   outf=open('spectrum/xspec-'+rc.dirname+'-%03d.dat'%it,'w')
   for i in range(len(kx)):
      print(kx[i],xkbx[i],xkby[i],xkbz[i],xkb[i], file=outf)
   outf.close()
   outf=open('spectrum/yspec-'+rc.dirname+'-%03d.dat'%it,'w')
   for i in range(len(ky)):
      print(ky[i],ykbx[i],ykby[i],ykbz[i],ykb[i], file=outf)
   outf.close()

rc.fin()
   
