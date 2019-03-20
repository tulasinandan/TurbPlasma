#!/usr/bin/env python
from pylab import *
import AnalysisFunctions as af
from subs import create_object, ask_for_steps

rc = create_object()
rc.vars2load(['bx','by','bz'])

bs,fs,step=ask_for_steps(rc.numslices)
nt=(fs-bs)/step
bxk=np.zeros(rc.nx/2)#;xk=np.zeros(rc.nx/2)
byk=np.zeros(rc.nx/2)
bzk=np.zeros(rc.nx/2)
xk = np.arange(rc.nx/2)*rc.dx

for it in range(bs,fs,step):
   print('time slice',it)
   rc.loadslice(it)
   bxk = bxk + af.fsdk(rc.bx)[1]
   byk = byk + af.fsdk(rc.by)[1]
   bzk = bzk + af.fsdk(rc.bz)[1]

bxk=bxk/nt; byk=byk/nt; bzk=bzk/nt

outf=open('sdkurtosis.'+rc.dirname+'.dat','w')
for i in range(rc.nx/2):
   print(xk[i],bxk[i],byk[i],bzk[i], file=outf) 
outf.close()

rc.fin()
   
