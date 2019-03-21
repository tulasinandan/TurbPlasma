#!/usr/bin/env python
#1)run 2)B0 3)betai 4)betae 5)Z0 6)L 7)runtime 8)lambda_c 9)T_p(0) 10)T_p(final) 11)T_e(0) 12)T_e(final) 13)Qi/Qe 14)tci/tnl 15) deltaE/deltaT

def msparams(rc,it,ft):
   import numpy as np
   from Codes.lccor import lccor
   rc.loadenergies()
   idxi = np.argmin(np.abs(rc.ta-it)); idxf=np.argmin(np.abs(rc.ta-ft))
   idxm = np.int(rc.t[np.argmin(np.abs(rc.ta-0.5*(it+ft)))]/rc.dtmovie)
   nstep= len(rc.ta)
   tt,lxc,lyc,lc=lccor(rc,bs=idxm,fs=idxm+1,step=1)
   b0=np.sqrt(rc.b0x**2+rc.b0y**2+rc.b0z**2)
  #Z0=np.sqrt(0.5*rc.edz[0])
   Z0=2.*np.sqrt(rc.edz[0])
   va=b0/np.sqrt(rc.n_0)
  #tnlotci=(va/Z0)*(lc/(2*np.pi))**(1./3.)
   tnlotci=(va/Z0)*(lc)**(1./3.)
   dedt=(rc.edz[idxf]-rc.edz[idxi])/(ft-it)
   rundetails=np.array([\
   b0, rc.betai, rc.betae, Z0, rc.lx,\
   rc.ta.max(), lc, rc.eip[idxi], rc.eip[idxf], rc.eep[idxi], rc.eep[idxf],\
   (rc.eip[idxf]-rc.eip[idxi])/(rc.eep[idxf]-rc.eep[idxi]), 1/tnlotci,dedt
   ])
   return rundetails

if __name__=="__main__":
   from subs import create_object
   rc=create_object()
   msp=msparams(rc,2.,5.)
   print(msp)
