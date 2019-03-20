#!/usr/bin/env python
def maxkurtseries(rc,it=0.,ft=.1,stept=1,bl=0.2,fl=5.,stepl=1):
   import numpy as np
   import AnalysisFunctions as af

   rc.vars2load(['bx','by'])
   rc.loadenergies()

   bst=np.int(rc.t[np.argmin(np.abs(rc.ta-it))]/rc.dtmovie)
   fst=np.int(rc.t[np.argmin(np.abs(rc.ta-ft))]/rc.dtmovie)
   if fl > 0.5*rc.lx:
      fl=0.5*rc.lx
   if bl < rc.dx:
      bl=rc.dx
   bsl=int(bl/rc.dx); fsl=int(fl/rc.dx)

   sclmx=np.zeros((fst-bst)/stept)
   krtmx=np.zeros((fst-bst)/stept)
   tsris=np.zeros((fst-bst)/stept)

   for it in range(bst,fst,stept):
      print(it)
      #Find the appropriate index for movie
      idx = (it-bst)/stept

      #Load the corresponding time slice
      rc.loadslice(it)

      #Assign the appropriate nonlinear time to the timeseries
      tsris[idx] = rc.ta[np.argmin(np.abs(rc.t-rc.time))]

      #Compute the scale dependent kurtosis between bl,fl
      rx,sdx=af.sdk(rc.bx,bs=bsl,fs=fsl,step=stepl,dx=rc.dx,ax=0)
      ry,sdy=af.sdk(rc.by,bs=bsl,fs=fsl,step=stepl,dx=rc.dx,ax=1)

      #Find the location of maximum value of kurtosis between sdx,sdy
#     xx=np.argmin(np.abs(np.array([sdx.max(),sdy.max()])))
      tsdk=np.array([sdx,sdy])
      rr =np.array([rx,ry])
      xx=np.unravel_index(tsdk.argmax(),tsdk.shape)

      sclmx[idx]=rr[xx]
      krtmx[idx]=tsdk[xx]

#     if xx[1] ==0:
#        sclmx[idx] = rx[xx[0]]
#        krtmx[idx] = sdx.max()
#     else:
#        sclmx[idx] = ry[xx[0]]
#        krtmx[idx] = sdy.max()
   return tsris,sclmx,krtmx,np.mean(krtmx)

if __name__ == "__main__":
   from subs import create_object,ask_for_steps
   rc=create_object()
   rc.vars2load(['bx','by'])

   tim,scl,krtmx,avgkrt=maxkurtseries(rc,it=0.,ft=5.,stept=2)
   ofl=open('maxkurtdi.'+rc.dirname+'.dat','w')
   for i in range(len(scl)):
      print(tim,scl[i],krtmx[i], file=ofl)
   ofl.close()
   print(avgkrt)
