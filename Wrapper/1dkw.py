from pylab import *
from subs import create_object,ask_for_steps

rc = create_object()

rc.vars2load(['bx','by','bz','n'])

recalc=input('Recalculate spectra? ')
if (recalc == 'y'):
   print(''); print('')
   bs,fs,step=ask_for_steps(rc.numslices)
   nslices=fs-bs
   windowff=np.zeros(nslices)
   for i in range(nslices):
      tht = 2*pi*i/(nslices-1)
   #  Hanning
   #  windowff[i]=0.5*(1-cos(tht))
   
   #  Blackman
   #  windowff[i]=0.42659-0.49656*cos(tht)+0.076849*cos(2*tht)
   
   #  BlackmanHarris
      windowff[i]=0.35875-0.48829*cos(tht)+0.14128*cos(2*tht)-0.01168*cos(3*tht)
   
   #  Flat top Window
   #  windowff[i] = 1. - 1.93*cos(tht) + 1.29*cos(2*tht) \
   #  - 0.388*cos(3*tht) + 0.028*cos(4*tht)
   
   #  Nuttall
   #  windowff[i] = 0.355768 - 0.487396*cos(tht) + \
   #                0.144232*cos(2*tht) - 0.012604*cos(3*tht)
   
   #  Welch
   #  windowff[i] = 1. - ((i - (nslices-1)/2.)/((nslices+1)/2.))**2
   
   
   bxt=np.zeros((nslices,rc.nx))
   byt=np.zeros((nslices,rc.nx))
   bzt=np.zeros((nslices,rc.nx))
   nt =np.zeros((nslices,rc.nx))
   
   for i in range(bs,fs,step):
         rc.loadslice(i)
         bxt[(i-bs)/step,:]=(rc.bx-mean(rc.bx))[:,0,0]*windowff[(i-bs)/step]
         byt[(i-bs)/step,:]=(rc.by-mean(rc.by))[:,0,0]*windowff[(i-bs)/step]
         bzt[(i-bs)/step,:]=(rc.bz-mean(rc.bz))[:,0,0]*windowff[(i-bs)/step]
         nt[(i-bs)/step,:] =(rc.n -mean(rc.n ))[:,0,0]*windowff[(i-bs)/step]
   
   
   tmp=fftn(nt)/sqrt(nslices*rc.nx)
   nft=fftshift(tmp)
   enk=0.5*abs(nft)**2
   
   tmp=fftn(bxt)/sqrt(nslices*rc.nx)
   bxft=fftshift(tmp)
   tmp=fftn(byt)/sqrt(nslices*rc.nx)
   byft=fftshift(tmp)
   tmp=fftn(bzt)/sqrt(nslices*rc.nx)
   bzft=fftshift(tmp)
   ebk=0.5*(abs(bxft)**2+abs(byft)**2+abs(bzft)**2)
   
   btime=round(bs*rc.movieout_full,2)
   ftime=round(fs*rc.movieout_full,2)
   k=fftshift(fftfreq(rc.nx)*2*pi*rc.nx/rc.lx)
   o=fftshift(fftfreq(nslices)*2*pi*nslices/(ftime-btime))

#dispersion=np.loadtxt('disp.dat')
#kd=dispersion[:,0]
#a=dispersion[:,2]
#f=dispersion[:,3]

### Voitenko's KAW relation
#kaw=np.zeros(len(k))
#ct=cos(0.9889032*pi/2); st=sin(0.9889032*pi/2)
#from scipy.special import i0 as i0
#for i in range(len(k)):
#   kkk=k[i]
#   kaw[i]=kkk*ct*kkk*st*sqrt(1/(1-i0(kkk**2*st**2)*exp(-kkk**2*st**2))+10)

figure(1)
clf()
#contourf(k,o,ebk,256)
tmp=np.clip(ebk,ebk.max()*1e-4,ebk.max())/ebk.max()
contourf(k,o,np.log10(tmp),32,cmap=cm.gist_ncar_r)
xlim([-7.,7.])
ylim([-10.,10.])
#xscale('log')
#yscale('log')
colorbar()
#plot(kd,a,'-o', kd,f,'-o')
ylabel('$\omega/\omega_{ci}$')
xlabel('$kd_i$')
title('$k-\omega$ plot for log($E_B$) between t=('+str(btime)+','+str(ftime)+')')
savefig('bkw.png')

figure(2)
clf()
#contourf(k,o,enk,256)
tmp=np.clip(enk,enk.max()*1e-4,enk.max())/enk.max()
contourf(k,o,np.log10(tmp),32,cmap=cm.gist_ncar_r)
xlim([-20.,20.])
ylim([-5.,5.])
#yscale('log')
#xscale('log')
#plot(kd,a,'-o', kd,f,'-o')
colorbar()
ylabel('$\omega/\omega_{ci}$')
xlabel('$kd_i$')
title('$k-\omega$ plot for log($E_n$) between t=('+str(btime)+','+str(ftime)+')')
savefig('nkw.png')

#show()

