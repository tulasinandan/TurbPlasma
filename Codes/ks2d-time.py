#!/usr/bin/env python
import AnalysisFunctions as af
import numpy as np
from subs import create_object,ask_for_steps,imss
import matplotlib.pyplot as plt

rc=create_object()
rc.vars2load(['bx','by','bz'])
rcd=rc.__dict__
bs,fs,step=ask_for_steps(rc.numslices)
nt=(fs-bs)/step

cut=input(" What part of domain? Enter nothing for whole domain")
if cut is None:
   cut=[0,rc.lx,0.,rc.ly]
else:
   cut=cut.split()
   cut[0]=float(cut[0]); cut[1]=float(cut[1])
   cut[2]=float(cut[2]); cut[3]=float(cut[3])

for it in range(bs,fs,step):
   print('time slice',it)
   rc.loadslice(it)
   k1,k2,ekx,eky,ekz,ekb=af.SpecVec2D(rc.bx,rc.by,rc.bz,axx=2,lx=rc.lx,ly=rc.ly,lz=rc.lz)
   if it == bs:
      vmax=np.log10(ekb).max(); vmin=vmax-6
   clf();imss(rcd,np.log10(ekb),cmap=cm.gist_rainbow_r,cut=cut,cbar=1,interpolation='none', vmin=vmin, vmax=vmax,lblsz='large')
   plt.savefig('2dspec-'+rc.dirname+'-'+str(rc.time)+'.png')

rc.fin()
