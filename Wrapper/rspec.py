#!/usr/bin/env python
from pylab import *
import AnalysisFunctions as af
from subs import create_object, ask_for_steps

rc = create_object()
rc.vars2load(['bx','by','bz'])

bs,fs,step = ask_for_steps(rc.numslices)

xbx=np.zeros(rc.nx/2)
xby=np.zeros(rc.nx/2)
xbz=np.zeros(rc.nx/2)
ybx=np.zeros(rc.ny/2)
yby=np.zeros(rc.ny/2)
ybz=np.zeros(rc.ny/2)
for it in range(bs,fs,step):
   print('time slice',it)
   rc.loadslice(it)

   kx,xkx,xky,xkz,xkb = af.ReducedSpecVec(rc.bx,rc.by,rc.bz,lenx=rc.lx,ax=0)
   xbx=xbx+xkx ; xby=xby+xky ; xbz=xbz+xkz

   ky,ykx,yky,ykz,ykb=af.ReducedSpecVec(rc.bx,rc.by,rc.bz,leny=rc.ly,ax=1)
   ybx=ybx+ykx ; yby=yby+yky ; ybz=ybz+ykz

xbz=xbz*float(step)/float(fs-bs)
xby=xby*float(step)/float(fs-bs)
xbx=xbx*float(step)/float(fs-bs)

ybz=ybz*float(step)/float(fs-bs)
yby=yby*float(step)/float(fs-bs)
ybx=ybx*float(step)/float(fs-bs)
outf=open('rxspec.'+rc.dirname+'.dat','w')
for i in range(len(kx)):
   print(kx[i],xbx[i],xby[i],xbz[i],xbx[i]+xby[i]+xbz[i], file=outf) 
outf.close()
outf=open('ryspec.'+rc.dirname+'.dat','w')
for i in range(len(ky)):
   print(ky[i],ybx[i],yby[i],ybz[i],ybx[i]+yby[i]+ybz[i], file=outf) 
outf.close()

rc.fin()
   
