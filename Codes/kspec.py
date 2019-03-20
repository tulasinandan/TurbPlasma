#!/usr/bin/env python

def kspec(rc):
   import AnalysisFunctions as af
   import numpy as np
   from subs import ask_for_steps

   rc.vars2load(['bx','by','bz'])
   bs,fs,step = ask_for_steps(rc.numslices)

   ebx=np.zeros(rc.nx/2+1)
   eby=np.zeros(rc.nx/2+1)
   ebz=np.zeros(rc.nx/2+1)
   for it in range(bs,fs,step):
      print('time slice',it)
      rc.loadslice(it)
      kk,ekx,eky,ekz,ekb=af.fperpspecvec(rc.bx,rc.by,rc.bz)
      ebx=ebx+ekx ; eby=eby+eky ; ebz=ebz+ekz
   
   kk=kk*2*np.pi/rc.lx; dk=kk[2]-kk[1]
   ebz=ebz*float(step)/(float(fs-bs)*dk)
   eby=eby*float(step)/(float(fs-bs)*dk)
   ebx=ebx*float(step)/(float(fs-bs)*dk)
   outf=open('kspec.'+rc.dirname+'.dat','w')
   for i in range(rc.nx/2):
      print(kk[i],ebx[i],eby[i],ebz[i], file=outf) 
   outf.close()
   
if __name__=="__main__":
   from subs import create_object
   rc=create_object()
   kspec(rc)
   rc.fin()
   
