#!/usr/bin/env python
def lccor(rc,bs=0,fs=1,step=1,kind='int'):
   import numpy as np
   from AnalysisFunctions import fcorr 
   ie=1/np.e
   rc.vars2load(['bx','by','bz'])
   tt  = np.zeros((fs-bs)/step)
   lxc = np.zeros((fs-bs)/step)
   lyc = np.zeros((fs-bs)/step)
   lc  = np.zeros((fs-bs)/step)
   for i in range(bs,fs,step):
      print(i); idx = (i-bs)/step
      rc.loadslice(i); 
      tt[idx] = rc.time
      rx,bxcor=fcorr(rc.bx,rc.bx,ax=0,dx=rc.dx)
      ry,bycor=fcorr(rc.by,rc.by,ax=1,dx=rc.dy)
      if kind == "ie":
         lxc[idx]=rx[abs(bxcor-ie).argmin()]
         lyc[idx]=ry[abs(bycor-ie).argmin()]
      elif kind == "int":
         lxc[idx]=np.sum(bxcor)*rc.dx
         lyc[idx]=np.sum(bycor)*rc.dy
      lc[idx] = 0.5*(lxc[idx]+lyc[idx])
      print(tt[idx],lxc[idx],lyc[idx],lc[idx])
   return tt,lxc,lyc,lc

if __name__ == '__main__':
   from subs import create_object, ask_for_steps
   rc=create_object()
   rc.vars2load(['bx','by','bz'])
   bs,fs,steps=ask_for_steps(rc.numslices)
   kind=input("Integral scale ('int') or 1/e scale ('ie')? ")
   tt,lxc,lyc,lc = lccor(rc,bs,fs,step,kind=kind)
   
   ofile=open('acorl.'+rc.dirname+'.dat','w')
   for i in range(len(lxc)):
      print(tt[i],lxc[i],lyc[i],lc[i], file=ofile)
   ofile.close()
