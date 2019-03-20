#!/usr/bin/env python
#
#  Calculating the Z+ & Z- values based on center of mass 
#  velocity instead of proton velocity.
#

def zpzmlclm(rc):
   import AnalysisFunctions as af
   import numpy as np
   from subs import ask_for_steps
   rc.vars2load(['bx','by','bz','jix','jiy','jiz','ni','jex','jey','jez','ne'])
   bs,fs,step=ask_for_steps(rc.numslices)
   #bs=0; fs=rc.numslices; step=1
   nt=(fs-bs)/step+1
   ie=1/np.e
   ezp=np.zeros(nt); ezm=np.zeros(nt)
   tt=np.zeros(nt); tnldi=np.zeros(nt)
   lcp=np.zeros(nt); lcm=np.zeros(nt)
   
   #time_series_output=raw_input("Output zp, zm etc as a time series? ")
   time_series_output = ''
   if time_series_output == 'y':
      zpkxt = np.zeros((rc.nx/2,nt))
      zpkyt = np.zeros((rc.nx/2,nt))
      zpkzt = np.zeros((rc.nx/2,nt))
      
      zmkxt = np.zeros((rc.nx/2,nt))
      zmkyt = np.zeros((rc.nx/2,nt))
      zmkzt = np.zeros((rc.nx/2,nt))
   
   for it in range(bs,fs,step):
      print('time slice',it)
      rc.loadslice(it)
      idx=(it-bs)/step
      rho=rc.ni+rc.ne*rc.m_e
      rc.jix=rc.jix/rc.ni; rc.jiy=rc.jiy/rc.ni;rc.jiz=rc.jiz/rc.ni
      rc.jex=rc.jex/rc.ne; rc.jey=rc.jey/rc.ne;rc.jez=rc.jez/rc.ne
   # 
      vcx = (rc.jix+rc.jex*rc.m_e)/(1+rc.m_e)
      vcy = (rc.jiy+rc.jey*rc.m_e)/(1+rc.m_e)
      vcz = (rc.jiz+rc.jez*rc.m_e)/(1+rc.m_e)
      rc.bx=(rc.bx-np.mean(rc.bx))/np.sqrt(rho)
      rc.by=(rc.by-np.mean(rc.by))/np.sqrt(rho)
      rc.bz=(rc.bz-np.mean(rc.bz))/np.sqrt(rho)
   #
      vcx=vcx-np.mean(vcx);vcy=vcy-np.mean(vcy);vcz=vcz-np.mean(vcz)
   #
      zpx=vcx+rc.bx; zpy=vcy+rc.by; zpz=vcz+rc.bz
   #  
      ezp[idx]=0.5*np.mean(rho*(np.abs(zpx)**2+np.abs(zpy)**2+np.abs(zpz)**2))
      zmx=vcx-rc.bx; zmy=vcy-rc.by; zmz=vcz-rc.bz
   # 
      ezm[idx]=0.5*np.mean(rho*(np.abs(zmx)**2+np.abs(zmy)**2+np.abs(zmz)**2))
   ## #Correlation for z+
     #kwave,zkx,zky,zkz,zk = af.PerpSpecVec(zpx,zpy,zpz)
      kwave,zkx,zky,zkz,zk = af.fperpspecvec(zpx,zpy,zpz)
      if time_series_output == 'y':
         zpkt[:,idx]=zk
      lcp[idx] = (2*np.pi/sum(zk[1:]))*sum(zk[1:]/kwave[1:])
   ## #Correlation for z-
     #kwave,zkx,zky,zkz,zk = af.PerpSpecVec(zmx,zmy,zmz)
      kwave,zkx,zky,zkz,zk = af.fperpspecvec(zmx,zmy,zmz)
      if time_series_output == 'y':
         zmkt[:,idx]=zk
      lcm[idx] = (2*np.pi/sum(zk[1:]))*sum(zk[1:]/kwave[1:])
   # 
      tt[idx]=round(it*rc.dtmovie,4)
      # wci*tnl(di) = Va/Z (L/di)^(1/3)
      # Matthaeus EA, ApJ 2014
      tnldi[idx]=(rc.B0*((lcp[idx]+lcm[idx])/2.)**(1./3.))/(2*np.sqrt(ezp[idx]+ezm[idx]))

   
   outf=open('zpzmlc.'+rc.dirname+'.dat','w')
   print('#','t','tt','ezp','ezm','lcp','lcm','tnldi', file=outf)
   for i in range(nt):
      print(tt[i],tt[i]*2*np.pi/rc.lx,ezp[i],ezm[i],lcp[i],lcm[i],tnldi[i], file=outf)
   outf.close()
   
   rc.fin()
   
   if time_series_output == 'y':
      zpkxt.tofile('zpkxt.'+rc.dirname+'.dat')
      zpkyt.tofile('zpkyt.'+rc.dirname+'.dat')
      zpkzt.tofile('zpkzt.'+rc.dirname+'.dat')
      zmkxt.tofile('zmkxt.'+rc.dirname+'.dat')
      zmkyt.tofile('zmkyt.'+rc.dirname+'.dat')
      zmkzt.tofile('zmkzt.'+rc.dirname+'.dat')
if __name__ == '__main__':
   from subs import create_object
   rc=create_object()
   zpzmlclm(rc)
