#!/usr/bin/env python
from pylab import *
import AnalysisFunctions as af
from subs import create_object, ask_for_steps

rc = create_object()

rc.vars2load(['jz'])
bs,fs,step=ask_for_steps(rc.numslices)
nt=(fs-bs)/step + 1
jzt=np.zeros((rc.nx,rc.ny,rc.nz,nt))

print 'What weight function to use?'
wght=int(raw_input())

for it in range(bs,fs,step):
   print 'time slice',it
   rc.loadslice(it)
   i = (it-bs)/step
   jzt[:,:,:,i]=rc.jz

jzb,jzpdf = af.calc_pdf(jzt,weight=wght)

outf=open('jzpdf.'+rc.dirname+'.dat','w')
for i in range(len(jzb)):
   print >> outf,jzb[i],jzpdf[i]
outf.close()

rc.fin()
   
