#!/usr/bin/env python
from pylab import *
from numpy.fft import *
from scipy.signal import detrend
from subs import create_object

rc=create_object()

ebdir=rc.rundir+'/staging/EBSeries'
arrs=['t','ex','ey','ez','bx','by','bz']
data=np.loadtxt(ebdir+'/EBseries.0000'); 
#idx=abs(data[:,0]*2*pi/rc.lx-5.0).argmin()
#ntimes=len(data[idx:,0])
ntimes=len(data[:,0])
for j in range(1,len(arrs)):
   exec('e'+arrs[j]+'=np.zeros(ntimes)')
for i in range(0,rc.nprocs,20):
   print('Processor ',str(i))
   file=ebdir+'/EBseries.%04d' % (i)
   data=np.loadtxt(file)
#  data=data[idx:,:]
#  wff=windowff(len(data[:,0]),kind="blackmannharris")
   t=data[:,0]
   for j in range(1,len(arrs)):
#     exec(arrs[j]+'=data[:,'+str(j)+']*wff')
#     exec(arrs[j]+'=ddetrend(data[:,'+str(j)+'],degree=int(float(ntimes)/40))')
      exec(arrs[j]+'=detrend(data[:,'+str(j)+'])')
      exec(arrs[j]+'f=fft('+arrs[j]+')/ntimes')
      exec('e'+arrs[j]+'t=0.5*abs('+arrs[j]+'f)**2')
      exec('e'+arrs[j]+'=e'+arrs[j]+'+e'+arrs[j]+'t')
for i in range(1,len(arrs)):
   exec('e'+arrs[j]+'=e'+arrs[j]+'/rc.nprocs')
ffq=fftfreq(ntimes,d=rc.dt)*2*pi
outf=open('Eulerian.'+rc.dirname+'.dat','w')
for i in range(ntimes):
   print(ffq[i],eex[i],eey[i],eez[i],ebx[i],eby[i],ebz[i], file=outf)
outf.close()

plot(ffq,eex+eey+eez, label='E')
plot(ffq,ebx+eby+ebz, label='B')
#xlim([ffq[2]-ffq[1],ffq.max()])
xlim([2.5e-2,3e2])
ylim([1e-7,5e2])
xscale('log')
yscale('log')
xlabel('$\omega$')
ylabel('$E(\omega)$')
legend(loc=3)
savefig('Eulerian-spec-'+rc.dirname+'.png')
