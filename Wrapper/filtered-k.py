#!/usr/bin/env python
from pylab import *
from numpy.fft import *
from subs import create_object,ask_for_steps

rc = create_object()

rc.vars2load(['bx','by','bz'])
bs,fs,step=ask_for_steps(rc.numslices)
nt=(fs-bs)/step
bxk=np.zeros(rc.nx/2)#;xk=np.zeros(rc.nx/2)
byk=np.zeros(rc.nx/2)
bzk=np.zeros(rc.nx/2)
xk = np.arange(rc.nx/2)*rc.dx

kx=fftshift(fftfreq(rc.nx))*rc.nx
ky=fftshift(fftfreq(rc.ny))*rc.ny
kz=fftshift(fftfreq(rc.nz))*rc.nz
kp=np.zeros((rc.nx,rc.ny))
for x in range(rc.nx):
   for y in range(rc.ny):
      kp[x,y]=sqrt(kx[x]**2+ky[y]**2)
dk=abs(kp[rc.nx/2,1]-kp[rc.nx/2,0])

kde=int(round(rc.lx/(sqrt(rc.m_e)*2*pi)))

for it in range(bs,fs,step):
   print('time slice',it)
   rc.loadslice(it)
   fbx = fftshift(fftn(rc.bx[:,:,0]))
   fby = fftshift(fftn(rc.by[:,:,0]))
   fbz = fftshift(fftn(rc.bz[:,:,0]))
   for x in range(rc.nx):
      for y in range(rc.ny):
         i=np.round(kp[x,y])
         if i > 5*kde:
            fbx[x,y] = complex(0,0)
            fby[x,y] = complex(0,0)
            fbz[x,y] = complex(0,0)
   bxf = real(ifftn(ifftshift(fbx)))
   byf = real(ifftn(ifftshift(fby)))
   bzf = real(ifftn(ifftshift(fbz)))
   bxk = bxk + sdk(bxf)[1]
   byk = byk + sdk(byf)[1]
   bzk = bzk + sdk(bzf)[1]

bxk=bxk/nt; byk=byk/nt; bzk=bzk/nt

outf=open('filtered-sdk.'+rc.dirname+'.dat','w')
for i in range(rc.nx/2):
   print(xk[i],bxk[i],byk[i],bzk[i], file=outf) 
outf.close()
rc.fin()
