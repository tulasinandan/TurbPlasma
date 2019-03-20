#!/usr/bin/env python
import numpy as np
from os.path import basename, realpath, exists
   
###-------------------------------------------------------------------------------
###
### THESE ROUTINES BELONG TO A MORE GENERAL LIBRARY BUT COPYING THEM HERE FOR
### COMPACTNESS FOR NOW. I HAVE THEM IN A LIBRARY CALLED AnalysisFunctions.
### I WILL REMOVE THEM WHEN I SHARE WITH YOU THE COMPLETE PACKAGE.
###                                     TULASI - May 24, 2017
###
###-------------------------------------------------------------------------------
##
## Second derivative for periodic arrays
##
def pdderiv(ar,dx=1.,ax=0,order=4,smth=None):
   """
      pderiv gives the double partial derivative
      of a periodic array along a given axis.

      Inputs:
         ar - The input array
         dx - Grid spacing, defaults to 1.
         ax - Axis along which to take the derivative
         order - Order of accuracy, (2,4) defaults to 2

      Output:
         dar - The derivative array
   """
   if smth is not None:
      ar=gf(ar,sigma=smth)
   if order == 2:
      dar = (np.roll(ar,-1,axis=ax) - 2*ar + np.roll(ar,1,axis=ax))/dx**2
   elif order == 4:
      dar = (-np.roll(ar,-2,axis=ax) + 16*np.roll(ar,-1,axis=ax) - 30*ar + 16*np.roll(ar,1,axis=ax)-np.roll(ar,2,axis=ax))/(12*dx**2)

   return dar

##
## Curl of a periodic array
##
def pcurl(arx,ary,arz,dx=1.,dy=1.,dz=1.,smth=None):
   return pderiv(arz,dx=dy,ax=1,smth=smth)-pderiv(ary,dx=dz,ax=2,smth=smth),\
          pderiv(arx,dx=dz,ax=2,smth=smth)-pderiv(arz,dx=dx,ax=0,smth=smth),\
          pderiv(ary,dx=dx,ax=0,smth=smth)-pderiv(arx,dx=dy,ax=1,smth=smth)

def pcurlx(ary,arz,dx=1.,dy=1.,dz=1.,smth=None):
   return pderiv(arz,dx=dy,ax=1,smth=smth)-pderiv(ary,dx=dz,ax=2,smth=smth)

def pcurly(arz,arx,dx=1.,dy=1.,dz=1.,smth=None):
   return pderiv(arx,dx=dz,ax=2,smth=smth)-pderiv(arz,dx=dx,ax=0,smth=smth)

def pcurlz(arx,ary,dx=1.,dy=1.,dz=1.,smth=None):
   return pderiv(ary,dx=dx,ax=0,smth=smth)-pderiv(arx,dx=dy,ax=1,smth=smth)

### ----------------------------------------------------------------------------------------
###
###  END OUTSIDE ROUTINES. WILL BE REMOVED IN THE COMPLETE PACKAGE SHARED LATER
###
### ----------------------------------------------------------------------------------------

class spc(object):
   """SPectral Code Reader: First Cut
	First Cut - 2017/01/13
   """

   def __init__(self,shelldirname=None):
      # If no rundir specified
      if shelldirname is None: 
         shelldirname = input('Please enter the rundir: ') 
      self.rundir = realpath(shelldirname)
      self.dirname= basename(self.rundir)
      self.primitives=['ax','ay','az','vx','vy','vz']
      self.derived={'bx':['ay','az'],'by':['ax','az'],'bz':['ax','ay'],\
        'ex':['by','bz','vy','vz'],'ey':['bx','bz','vx','vz'],\
        'ez':['bx','by','vx','vy']}
      self.allvars=self.primitives+['bx','by','bz','ex','ey','ez']
      self._readinit()

   def _readinit(self):
      f=open(self.rundir+'/mhd3-00-000.dat','rb')
      dum0=np.fromfile(f,dtype='int32',count=1)
      ih=np.fromfile(f,dtype='i4',count=9)
      f.close()
      self.nxc=ih[0]+1;self.nyc=ih[1];self.nzp=ih[2]
      self.nx=ih[3]; self.ny=ih[4]; self.nz=ih[5]
      self.lbox=2*np.pi*ih[6]; self.nprocs=ih[8]
      self._readglobs()
      self.xx=np.linspace(0,self.lbox,self.nx)
      self.yy=np.linspace(0,self.lbox,self.ny)
      self.zz=np.linspace(0,self.lbox,self.nz)
      self.dx=self.lbox/self.nx
      self.dy=self.lbox/self.ny
      self.dz=self.lbox/self.nz
      
   def _readglobs(self):
      fl=open(self.rundir+'/globs.dat','r')
      self.comment=fl.readline().strip(' \t\r\n')
      tmp=fl.readline().split(); self.ntglobs=np.int(tmp[0]); self.numglobs=np.int(tmp[1])
      globs=[]
      for i in fl.readlines():
         globs+=i.split()
      globs=np.asarray([float(i) for i in globs]).reshape((self.ntglobs,self.numglobs))
      fl.close()
      half_pi=np.pi/2.
      self.time      = globs[:,0]
      self.Ev        = globs[:,1]
      self.Eb        = globs[:,2]
      self.asqd      = globs[:,3]
      self.jsqd      = globs[:,4]
      self.enst      = globs[:,5]
      self.Hm        = globs[:,6]
      self.Hc        = globs[:,7]
      self.Hk        = globs[:,8]
      self.da0       = globs[:,9]
      self.emfx      = globs[:,10]
      self.emfy      = globs[:,11]
      self.emfz      = globs[:,12]
      self.Ev_2D     = globs[:,13]
      self.Eb_2D     = globs[:,14]
      self.asqd_2D   = globs[:,15]
      self.jsqd_2D   = globs[:,16]
      self.enst_2D   = globs[:,17]
      self.Hc_2D     = globs[:,18]
      self.jom       = globs[:,19]
      self.Hj        = globs[:,20]
      self.jom_2D    = globs[:,21]
      self.Lv_2D     = globs[:,22]/self.Ev_2D*half_pi
      self.Lb_2D     = globs[:,23]/self.Eb_2D*half_pi
      self.LHc_2D    = globs[:,24]*half_pi
      self.Lv        = globs[:,25]/self.Ev*half_pi
      self.Lb        = globs[:,26]/self.Eb*half_pi
      self.LHc       = globs[:,27]*half_pi
      self.SigVom    = globs[:,28]
      self.SigJB     = globs[:,29]
      self.SigVB     = globs[:,30]
      self.max_om123 = globs[:,31:37]
      self.max_j123  = globs[:,37:]

   def vars2load(self,v2lu):
      """
         Based on user input, find the dependencies and create the v2l
         array.
      """
      
      if len(v2lu) == 1:
         if v2lu[0] == 'min':
               self.vars2l=self.primitives
         else:
               self.vars2l=v2lu
      else:
         self.vars2l = v2lu

      v2l=v2lu[:]
      
      while any([x in self.derived for x in v2l]):
         toload=v2l[:]
         v2l=[]
         while len(toload) > 0:
            current=toload.pop()
            if current in self.primitives and current not in v2l:
              v2l.append(current)
            elif current in self.derived:
              for i in self.derived[current]:
                if i not in v2l:
                  v2l.append(i)
            elif current not in self.primitives+list(self.derived.keys()):
              print(current+' not implemented. Implement it!')
      for i in v2lu:
        if i in self.derived: 
          for j in self.derived[i]:
            if j not in v2l: v2l.append(j)
      for i in v2lu:
        if i in self.derived and i not in v2l: v2l.append(i)
      self.vars2l=v2l
#     self.vars2l=sorted(v2l,key=self.allvars.index)
      
      for i in self.primitives:
         self.__dict__[i+'c']=np.zeros((self.nxc,self.ny,self.nz),dtype='complex')
         self.__dict__[i]    =np.zeros((self.nx, self.ny, self.nz))
      for i in self.vars2l:
         if i in self.derived:
            self.__dict__[i]=np.zeros((self.nx,self.ny,self.nz))

   def loadslice(self,timeslice):
#     from scipy.io import FortranFile
      for i in range(self.nprocs):
         f=open('mhd3-%02d-%03d.dat'%(timeslice,i),'rb')      
         d={}
         dum0=np.fromfile(f,dtype='int32',count=1)
         ih=np.fromfile(f,dtype='i4',count=9)
         dum1=np.fromfile(f,dtype='int32',count=1)
         ncount=(ih[0]+1)*ih[1]*ih[2]
         
         dum2=np.fromfile(f,dtype='int32',count=1)
         d['time']=np.fromfile(f,dtype='float',count=1)
         d['ax']  =np.fromfile(f,dtype='complex',count=ncount).reshape((ih[0]+1,ih[1],ih[2]),order='F')
         d['ay']  =np.fromfile(f,dtype='complex',count=ncount).reshape((ih[0]+1,ih[1],ih[2]),order='F')
         d['az']  =np.fromfile(f,dtype='complex',count=ncount).reshape((ih[0]+1,ih[1],ih[2]),order='F')
         dum2=np.fromfile(f,dtype='int32',count=1)
         
         dum2=np.fromfile(f,dtype='int32',count=1)
         d['time']=np.fromfile(f,dtype='float',count=1)
         d['vx']  =np.fromfile(f,dtype='complex',count=ncount).reshape((ih[0]+1,ih[1],ih[2]),order='F')
         d['vy']  =np.fromfile(f,dtype='complex',count=ncount).reshape((ih[0]+1,ih[1],ih[2]),order='F')
         d['vz']  =np.fromfile(f,dtype='complex',count=ncount).reshape((ih[0]+1,ih[1],ih[2]),order='F')
         dum2=np.fromfile(f,dtype='int32',count=1)
         f.close()

         for j in self.primitives:
            self.__dict__[j+'c'][:,:,i*ih[2]:(i+1)*ih[2]]=d[j]
      print('Loaded primitives in Fourier Space')

      for i in self.primitives:
         self.__dict__[i] = np.fft.irfftn(self.__dict__[i+'c'].T)
      print('Transformed primitives to real space')
      
      for i in self.vars2l:
         if i in self.derived:
            self.__dict__[i] = self._derivedv(i)
            

   def _derivedv(self,varname):
     #import AnalysisFunctions as af
      if varname == 'bx':
        #return af.pcurlx(self.ay,self.az)
         return pcurlx(self.ay,self.az)
      if varname == 'by':
        #return af.pcurly(self.az,self.ax)
         return pcurly(self.az,self.ax)
      if varname == 'bz':
        #return af.pcurlz(self.ax,self.ay)
         return pcurlz(self.ax,self.ay)
      if varname == 'ex':
         return self.vy*self.bz - self.vz*self.by
      if varname == 'ey':
         return self.vz*self.bx - self.vx*self.bz
      if varname == 'ez':
         return self.vx*self.by - self.vy*self.bx 
