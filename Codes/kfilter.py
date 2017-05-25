def kfilter(ar,kf):
   nx=shape(ar)[0];kx=fftshift(fftfreq(nx))*nx
   ny=shape(ar)[1];ky=fftshift(fftfreq(ny))*ny
   nz=shape(ar)[2];kz=fftshift(fftfreq(nz))*nz
   km=np.zeros((nx,ny,nz))
   for x in range(nx):
      for y in range(ny):
         for z in range(nz):
            km[x,y,z]=sqrt(kx[x]**2+ky[y]**2+kz[z]**2)
   
   fbx = fftshift(fftn(ar))
   for x in range(nx):
      for y in range(ny):
         for z in range(nz):
            i=np.round(kp[x,y,z])
            if i > kf:
               fbx[x,y,z] = complex(0,0)
   bxf = real(ifftn(ifftshift(fbx)))
   return bxf
