#!/usr/bin/env python
from pylab import *
import AnalysisFunctions as af
from subs import create_object,ask_for_steps

rc = create_object()

rc.vars2load(['jz'])
bs,fs,step=ask_for_steps(rc.numslices)
nt=(fs-bs)/step + 1
jzt=np.zeros((rc.nx,rc.ny,rc.nz,nt))

jzb=0; jzpdf=0
for it in range(bs,fs,step):
   print('time slice',it)
   rc.loadslice(it)
   tmp1,tmp2 = af.normhist(rc.jz,min=-20.,max=20.)
   jzb=jzb+tmp1
   jzpdf=jzpdf+tmp2
jzb=jzb/nt; jzpdf=jzpdf/nt

outf=open('jznh.'+rc.dirname+'.dat','w')
for i in range(len(jzb)):
   print(jzb[i],jzpdf[i], file=outf)
outf.close()

rc.fin()
   
