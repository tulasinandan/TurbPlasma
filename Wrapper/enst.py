#!/usr/bin/env python
from subs import create_object
import numpy as np
rc=create_object()
rc.vars2load(['ensti','enste','pali','pale','jx','jy','jz'])
ofl=open('enst.'+rc.dirname+'.dat','w')
nm=rc.lambdae/rc.dx
print('# Time\t j^2\t wi^2\t we^2\t Pi\t Pe', file=ofl)
for it in range(rc.numslices):
   print(it)
   rc.loadslice(it,smth=nm)
   print(rc.time, np.mean(rc.jx**2+rc.jy**2+rc.jz**2),\
                np.mean(rc.omix**2+rc.omiy**2+rc.omiz**2),\
                np.mean(rc.omex**2+rc.omey**2+rc.omez**2),\
                np.mean(rc.pali**2), np.mean(rc.pale**2), file=ofl)

ofl.close()
