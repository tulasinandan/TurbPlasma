#!/usr/bin/env python
#
#  Calculating the Z+ & Z- values based on center of mass 
#  velocity instead of proton velocity.
#
from pylab import *
from AnalysisFunctions import *
from subs import create_object

rc = create_object()

rc.vars2load(['bx','by','bz','jix','jiy','jiz','ni','jex','jey','jez','ne'])

nt=rc.numslices; bs=0; fs=rc.numslices; step=1
ie=1/e
ezp=np.zeros(nt); ezm=np.zeros(nt)
tt=np.zeros(nt)
lcp=np.zeros(nt); lcm=np.zeros(nt)

#zpkxt = np.zeros((rc.nx/2,nt))
#zpkyt = np.zeros((rc.nx/2,nt))
#zpkzt = np.zeros((rc.nx/2,nt))
#
#zmkxt = np.zeros((rc.nx/2,nt))
#zmkyt = np.zeros((rc.nx/2,nt))
#zmkzt = np.zeros((rc.nx/2,nt))

for it in range(bs,fs,step):
   print('time slice',it)
   rc.loadslice(it)
   rho=rc.ni+rc.ne*rc.m_e
   rc.jix=rc.jix/rc.ni; rc.jiy=rc.jiy/rc.ni;rc.jiz=rc.jiz/rc.ni
   rc.jex=rc.jex/rc.ne; rc.jey=rc.jey/rc.ne;rc.jez=rc.jez/rc.ne
# 
   vcx = (rc.jix+rc.jex*rc.m_e)/(1+rc.m_e)
   vcy = (rc.jiy+rc.jey*rc.m_e)/(1+rc.m_e)
   vcz = (rc.jiz+rc.jez*rc.m_e)/(1+rc.m_e)
   rc.bx=(rc.bx-mean(rc.bx))/sqrt(rho)
   rc.by=(rc.by-mean(rc.by))/sqrt(rho)
   rc.bz=(rc.bz-mean(rc.bz))/sqrt(rho)
#
   vcx=vcx-mean(vcx);vcy=vcy-mean(vcy);vcz=vcz-mean(vcz)
#
   zpx=vcx+rc.bx; zpy=vcy+rc.by; zpz=vcz+rc.bz
#  zpkxt[:,it] = spectrum(zpx)[1]
#  zpkyt[:,it] = spectrum(zpy)[1]
#  zpkzt[:,it] = spectrum(zpz)[1]
   ezp[(it-bs)/step]=0.5*mean(rho*(abs(zpx)**2+abs(zpy)**2+abs(zpz)**2))
   zmx=vcx-rc.bx; zmy=vcy-rc.by; zmz=vcz-rc.bz
#  zmkxt[:,it] = spectrum(zmx)[1]
#  zmkyt[:,it] = spectrum(zmy)[1]
#  zmkzt[:,it] = spectrum(zmz)[1]
   ezm[(it-bs)/step]=0.5*mean(rho*(abs(zmx)**2+abs(zmy)**2+abs(zmz)**2))
#  clf()
## #Correlation for z+
   kwave,zkx,zky,zkz,zk = specvec(zpx,zpy,zpz)
   lcp[(it-bs)/step] = (pi/sum(zk[1:]))*sum(zk[1:]/kwave[1:])
#  r,l=autocvec(zpx,zpy,zpz,ax=0,dx=rc.dx); lcp[(it-bs)/step]=r[abs(l-ie).argmin()]
#  if make_plot == 'y':
#     subplot(2,1,1)
#     title('t = '+str(it*rc.movieout_full*2*pi/rc.lx))
#     plot(r,l,label='$z^+$')
#     ylabel('Correlation')
#     legend(loc=3)
## #Correlation for z-
   kwave,zkx,zky,zkz,zk = specvec(zmx,zmy,zmz)
   lcm[(it-bs)/step] = (pi/sum(zk[1:]))*sum(zk[1:]/kwave[1:])
#  r,l=autocvec(zmx,zmy,zmz,ax=0,dx=rc.dx); lcm[(it-bs)/step]=r[abs(l-ie).argmin()]
#  if make_plot == 'y':
#     subplot(2,1,2)
#     plot(r,l,label='$z^-$')
#     ylabel('Correlation')
#     legend(loc=3)
#     savefig('lplm%03d.png' %it)
   tt[(it-bs)/step]=round(it*rc.movieout_full,4)

outf=open('zpzmlc.'+rc.dirname+'.dat','w')
print('#','t','tt','ezp','ezm','lcp','lcm', file=outf)
for i in range(nt):
   print(tt[i],tt[i]*2*pi/rc.lx,ezp[i],ezm[i],lcp[i],lcm[i], file=outf)
outf.close()

rc.fin()

#zpkxt.tofile('zpkxt.'+rc.dirname+'.dat')
#zpkyt.tofile('zpkyt.'+rc.dirname+'.dat')
#zpkzt.tofile('zpkzt.'+rc.dirname+'.dat')
#zmkxt.tofile('zmkxt.'+rc.dirname+'.dat')
#zmkyt.tofile('zmkyt.'+rc.dirname+'.dat')
#zmkzt.tofile('zmkzt.'+rc.dirname+'.dat')
