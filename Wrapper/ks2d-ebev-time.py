#!/usr/bin/env python
import AnalysisFunctions as af
import numpy as np
from subs import create_object,ask_for_steps,imss
import matplotlib.pyplot as plt

rc=create_object()
rc.vars2load(['bx','by','bz','jix','jiy','jiz','ni'])
rcd=rc.__dict__
bs,fs,step=ask_for_steps(rc.numslices)
nt=(fs-bs)/step

cut=input(" What part of domain? Enter nothing for whole domain")
if cut=='':
   cut=None
else:
   cut=cut.split()
   cut[0]=float(cut[0]); cut[1]=float(cut[1])
   cut[2]=float(cut[2]); cut[3]=float(cut[3])

for it in range(bs,fs,step):
   print('time slice',it)
   rc.loadslice(it)
   rc.computevars(['vi'])
   k1,k2,ekx,eky,ekz,ekb=af.SpecVec2D(rc.bx,rc.by,rc.bz,axx=2,lx=rc.lx,ly=rc.ly,lz=rc.lz)
   k1,k2,ekx,eky,ekz,ekv=af.SpecVec2D(rc.vix,rc.viy,rc.viz,axx=2,lx=rc.lx,ly=rc.ly,lz=rc.lz)
   if it == bs:
      kcd={'xx':k1,'yy':k2}
      vmax=np.log10(ekb+ekv).max(); vmin=vmax-6
   plt.clf();ax=plt.gca(); imss(kcd,np.log10(ekb+ekv),cmap=plt.cm.pink_r,\
         ax=ax,cut=cut,cbar=1,interpolation='none', vmin=vmin, vmax=vmax)
   ax.set_title('$t\omega_{ci}$ = '+str(round(rc.time,3)))
   ax.set_xlabel('$k_x d_i$')
   ax.set_ylabel('$k_y d_i$')
   plt.savefig('2dspec-'+rc.dirname+'-%04d.png'%it,dpi=300,bbox_inches='tight')

rc.fin()
