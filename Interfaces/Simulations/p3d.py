#!/usr/bin/env python
###########################################################################
##
## The class to handle P3D simulation data in python. An object is created 
## that has the run parameters attached with it. The basic usage of the 
##  class can be done in the following ways.
##
## Example 1:
## > from p3d import p3d
## > r1=p3d('OTV')
## > r1.vars2load(['bx','by','bz'])
## > for i in range(20):
## >     r1.loadslice(i)
## >     DO WHATEVER!!!
## > r1.fin()
##
##
## Example 2:
## > from p3d import p3d
## > r1=p3d('OTV')
## > bxf=open('bx','rb')
## > bx=r1.readslice(bxf,10)
## > DO WHATEVER!!!
## > r1.fin()
##
## Methods in the class:
## 1) self.__init__(): Initialize the run object and load run parameters
## 2) self.print_params(): Print the run parameters
## 3) self.vars2load(): Define what variables to load, create variables
##    and open files for reading data.
## 4) self.readslice(): Read a time slice from an open file.
## 5) self.loadslice(): Load a snapshot of time for all the defined variables.
## 6) 
##
##                                        Tulasi Nandan Parashar
##                                        2014/08/23
## Hat tip to ColbyCH for teaching me the basics of classes in python.
##
###########################################################################

import numpy as np
from os.path import basename, realpath, exists
import TurbPlasma.Analysis.Simulations as af
from scipy.ndimage import gaussian_filter as gf
from TurbPlasma.Utilities.subs import calc_dep
   
class p3d(object):
   """p3d object:
         Tulasi Parashar's version to read data from
         P3D github version.
         Created on 08/22/2014
         Last modified on 09/09/2015
         Last modified on 10/13/2016
   """
   def __init__(self,shelldirname=None,filenum=None):
      # If no rundir specified
      if shelldirname is None: 
         shelldirname = input('Please enter the rundir: ') 
      self.rundir = realpath(shelldirname)
      self.dirname= basename(self.rundir)

      # If filenum not given
      if filenum is None:
         self.filenum=input("Please enter the file number to load (e.g. 000): ")
      else:
         self.filenum=filenum

      # Check where the paramfile is
      if exists(self.rundir+'/param_'+self.dirname):
         self.paramfile=self.rundir+'/param_'+self.dirname
      elif exists(self.rundir+'/staging/param_'+self.dirname):
         self.paramfile=self.rundir+'/staging/param_'+self.dirname
      elif exists(self.rundir+'/paramfile'):
         self.paramfile=self.rundir+'/paramfile'
      else:
         raise ValueError('Paramfile not found in '+self.dirname)
      # load parameters
      self.params=loadparams(self.paramfile)
      for i in list(self.params.keys()):
         self.__dict__[i]=self.params[i]

      self.primitives=['bx','by','bz','ex','ey','ez','jix','jiy','jiz','jex','jey'\
      ,'jez','jx','jy','jz','ni','pixx','piyy','pizz','pixy','pixz','piyz'\
      ,'pexx','pexy','pexz','peyy','peyz','pezz','ne','rho']
      self.derived={\
      'tex':['pexx','ne'],'tey':['peyy','ne'],'tez':['pezz','ne'],\
      'te':['tex','tey','tez'],\
      'tix':['pixx','ni'],'tiy':['piyy','ni'],'tiz':['pizz','ni'],\
      'ti':['tix','tiy','tiz'],\
      'vix':['jix','ni'],'viy':['jiy','ni'],'viz':['jiz','ni'],\
      'vi':['vix','viy','viz'],\
      'vex':['jex','ne'],'vey':['jey','ne'],'vez':['jez','ne'],\
      've':['vex','vey','vez'],\
      'omix':['viy','viz'],'omiy':['viz','vix'],'omiz':['vix','viy'],\
      'omi':['omix','omiy','omiz'],'ensti':['omi'],'pali':['omi'],\
      'omex':['vey','vez'],'omey':['vex','vez'],'omez':['vex','vey'],\
      'ome':['omex','omey','omez'],'enste':['ome'],'pale':['ome'],\
      'omx':['cmy','cmz'],'omy':['cmz','cmx'],'omz':['cmx','cmy'],\
      'om':['omx','omy','omz'], 'enst':['om'], 'pal':['om'],\
      'dui':['vix','viy','viz'],'due':['vex','vey','vez'],\
      'den':['ni','ne'],\
      'cmx':['vix','vex'],'cmy':['viy','vey'],'cmz':['viz','vez'],\
      'zpx':['bx','cmx','den'],'zpy':['by','cmy','den'],'zpz':['bz','cmz','den'],\
      'zmx':['bx','cmx','den'],'zmy':['by','cmy','den'],'zmz':['bz','cmz','den'],\
      'zpzm':['zpx','zpy','zpz','zmx','zmy','zmz']\
      }

      self.allvars=self.primitives+list(self.derived.keys())


####
#### Method to print the parameters associated with the run
####
   def print_params(self):
      """
         A quick method to print the parameters and variables attached with
         the p3d run object.
      """
      for i in list(self.params.keys()):
         print(i,' = ',self.params[i])
####
#### Method to define the variables to load, create variables and
#### open corresponding files.
#### 
   def vars2load(self,v2lu):
      """
         Define the variables to load, define corresponding numpy arrays &
         open the files
      """
      if len(v2lu) == 1:
         if v2lu[0] == 'min':
            self.vars2l=['bx','by','jix','jiy','jz','ni']
         elif v2lu[0] == 'prim':
            self.vars2l=self.primitives
         elif v2lu[0] == 'all':
            self.vars2l=self.allvars
         else:
            # Compute all the dependencies and set variables 2 load 
            self.vars2l=calc_dep(v2lu,self.primitives,self.derived)
      else:
         # Compute all the dependencies and set variables 2 load 
         self.vars2l=calc_dep(v2lu,self.primitives,self.derived)

      # If byte or double byte data, open the log file.
      if self.data_type in ("b", "bb"):
         if exists(self.rundir+'/log'):
            print(self.rundir+'/log')
            self.logfile=open(self.rundir+"/log","r")
         else:
            print(self.rundir+'/staging/movie.log.'+self.filenum)
            self.logfile=open(self.rundir+"/staging/movie.log."+self.filenum,"r")
         self.szl=np.size(self.logvars)
         self.alllogvals=np.loadtxt(self.logfile)

      # Create arrays and open files
      for i in self.vars2l:
         self.__dict__[i]=np.array((self.nx,self.ny,self.nz))
         if i in self.primitives:
            if exists(self.rundir+'/'+i):
               self.__dict__[i+'f']=open(self.rundir+'/'+i,"rb")
            else:
               self.__dict__[i+'f']=open(self.rundir+'/staging/movie.'+i+\
                                    '.'+str(self.filenum),"rb")
      # If 'b' or 'bb' data type, load minmax for each loaded variable
      if self.data_type in ('b', 'bb'):
         for i in self.vars2l:
            if i in self.primitives:
               self.__dict__[i+'minmax']=self.alllogvals[self.logvars.index(i)\
                                         ::len(self.logvars),:]
      # Find out the number of slices for the open file
      self.__dict__[self.vars2l[0]+'f'].seek(0,2)
      filesize=self.__dict__[self.vars2l[0]+'f'].tell()
      numbersize=2**['b','bb','f','d'].index(self.data_type)
      self.numslices=filesize/(self.nx*self.ny*self.nz*numbersize)
####
#### Method to read a particular time slice from an open file.
#### 
   def readslice(self,f,timeslice,v=""):
      """
         This method reads a particular slice of time from a given file. The
         explicit inputs are file name, time slice and data type. It is used 
         as:
            output=self.readslice(filename,time)
      """
      ### Byte data #########################################################
      if self.data_type == 'b':
         f.seek(timeslice*self.nx*self.ny*self.nz)
         field = np.fromfile(f,dtype='uint8',count=self.nx*self.ny*self.nz)
         field = np.reshape(field,(self.nx,self.ny,self.nz),order='F')
         minmax=self.__dict__[v+'minmax'][timeslice]
         field = minmax[0]+(minmax[1]-minmax[0])*field/255.
      ### Double byte data #################################################
      elif self.data_type == 'bb':
         f.seek(2*timeslice*self.nx*self.ny*self.nz)
         field = np.fromfile(f,dtype='int16',count=self.nx*self.ny*self.nz)
         field = np.reshape(field,(self.nx,self.ny,self.nz),order='F')
         minmax=self.__dict__[v+'minmax'][timeslice]
         field = minmax[0]+(minmax[1]-minmax[0])*(field+32678.)/65535.
      ### Four byte (single precision) data ################################
      elif self.data_type == 'f':
         f.seek(4*timeslice*self.nx*self.ny*self.nz)
         field = np.fromfile(f,dtype='float32',count=self.nx*self.ny*self.nz)
         field = np.reshape(field,(self.nx,self.ny,self.nz),order='F').astype('float64')
      ### Double precision data ############################################
      elif self.data_type == 'd':
         f.seek(8*timeslice*self.nx*self.ny*self.nz)
         field = np.fromfile(f,dtype='float64',count=self.nx*self.ny*self.nz)
         field = np.reshape(field,(self.nx,self.ny,self.nz),order='F')
      return field
####
#### Method to load time slices for the loaded variables.
####
   def loadslice(self,it,smth=None):
      """
         Load the variables initialized by self.vars2load()
      """
      for i in self.vars2l:
         if i in self.primitives:
            self.__dict__[i]=self.readslice(self.__dict__[i+'f'],it,i)
      for i in self.vars2l:
         if i in self.derived:
           #self.__dict__[i] = self._derivedv(i)
            self._derivedv(i)

      self.mmd={}
      for i in self.vars2l:
         if smth is not None:
            self.__dict__[i]=gf(self.__dict__[i],sigma=smth)
         self.mmd[i]=[self.__dict__[i].min(),self.__dict__[i].max()]
      self.time = it*self.dtmovie

####
#### Method to add attributes to the object
####
   def addattr(self,key,val):
      for i in key:
         print('Adding '+i) 
         self.__dict__[i]=val[key.index(i)]
         if isinstance(val[key.index(i)],np.ndarray):
            self.mmd[i]=[self.__dict__[i].min(),self.__dict__[i].max()]

####
#### Method to compute derived quantities
####
   def _derivedv(self,varname):
      if varname == 'tix'   : self.tix    = self.pixx/self.ni
      if varname == 'tiy'   : self.tiy    = self.piyy/self.ni
      if varname == 'tiz'   : self.tiz    = self.pizz/self.ni
      if varname == 'ti'    : self.ti     = (self.tix+self.tiy+self.tiz)/3.
      if varname == 'tex'   : self.tex    = self.pexx/self.ne
      if varname == 'tey'   : self.tey    = self.peyy/self.ne
      if varname == 'tez'   : self.tez    = self.pezz/self.ne
      if varname == 'te'    : self.te     = (self.tex+self.tey+self.tez)/3.
      if varname == 'vix'   : self.vix    = self.jix/self.ni
      if varname == 'viy'   : self.viy    = self.jiy/self.ni
      if varname == 'viz'   : self.viz    = self.jiz/self.ni
      if varname == 'omix'  : self.omix   = af.pcurlx(self.viy,self.viz,dx=self.dx,dy=self.dy,dz=self.dz)
      if varname == 'omiy'  : self.omiy   = af.pcurly(self.viz,self.vix,dx=self.dx,dy=self.dy,dz=self.dz)
      if varname == 'omiz'  : self.omiz   = af.pcurlz(self.vix,self.viy,dx=self.dx,dy=self.dy,dz=self.dz)
      if varname == 'omi'   : pass
      if varname == 'ensti' : self.ensti  = self.omix**2+self.omiy**2+self.omiz**2
      if varname == 'dui'   : self.dui    = af.pdiv(self.vix,self.viy,self.viz,dx=self.dx,dy=self.dy,dz=self.dz)
      if varname == 'vex'   : self.vex    = -self.jex/self.ne
      if varname == 'vey'   : self.vey    = -self.jey/self.ne
      if varname == 'vez'   : self.vez    = -self.jez/self.ne
      if varname == 'omex'  : self.omex   = af.pcurlx(self.vey,self.vez,dx=self.dx,dy=self.dy,dz=self.dz)
      if varname == 'omey'  : self.omey   = af.pcurly(self.vez,self.vex,dx=self.dx,dy=self.dy,dz=self.dz)
      if varname == 'omez'  : self.omez   = af.pcurlz(self.vex,self.vey,dx=self.dx,dy=self.dy,dz=self.dz)
      if varname == 'ome'   : pass
      if varname == 'enste' : self.enste  = self.omex**2+self.omey**2+self.omez**2
      if varname == 'due'   : self.due    = af.pdiv(self.vex,self.vey,self.vez,dx=self.dx,dy=self.dy,dz=self.dz)
      if varname == 'cmx'   : self.cmx    = (self.vix+self.m_e*self.vex)/(1+self.m_e)
      if varname == 'cmy'   : self.cmy    = (self.viy+self.m_e*self.vey)/(1+self.m_e)
      if varname == 'cmz'   : self.cmz    = (self.viz+self.m_e*self.vez)/(1+self.m_e)
      if varname == 'omx'   : self.omx    = af.pcurlx(self.cmy,self.cmz,dx=self.dx,dy=self.dy,dz=self.dz)
      if varname == 'omy'   : self.omy    = af.pcurly(self.cmz,self.cmx,dx=self.dx,dy=self.dy,dz=self.dz)
      if varname == 'omz'   : self.omz    = af.pcurlz(self.cmx,self.cmy,dx=self.dx,dy=self.dy,dz=self.dz)
      if varname == 'om'    : pass
      if varname == 'enst'  : self.enst   = self.omx**2+self.omy**2+self.omz**2
      if varname == 'den'   : self.den    = self.ni+self.m_e*self.ne
      if varname == 'zpx'   : self.zpx    = self.bx/np.sqrt(self.den) + self.cmx
      if varname == 'zpy'   : self.zpy    = self.by/np.sqrt(self.den) + self.cmy
      if varname == 'zpz'   : self.zpz    = self.bz/np.sqrt(self.den) + self.cmz
      if varname == 'zmx'   : self.zmx    = self.bx/np.sqrt(self.den) - self.cmx
      if varname == 'zmy'   : self.zmy    = self.by/np.sqrt(self.den) - self.cmy
      if varname == 'zmz'   : self.zmz    = self.bz/np.sqrt(self.den) - self.cmz
      if varname == 'zpzm'  : pass

      if varname == 'pali'  : 
         tmp = af.pcurl(self.omix,self.omiy,self.omiz,dx=self.dx,dy=self.dy,dz=self.dz)
         self.pali   = 0.5*(tmp[0]**2+tmp[1]**2+tmp[2]**2); tmp=None
      if varname == 'pale'  : 
         tmp = af.pcurl(self.omex,self.omey,self.omez,dx=self.dx,dy=self.dy,dz=self.dz)
         self.pale   = 0.5*(tmp[0]**2+tmp[1]**2+tmp[2]**2); tmp=None
      if varname == 'pal'  : 
         tmp = af.pcurl(self.omx,self.omy,self.omz,dx=self.dx,dy=self.dy,dz=self.dz)
         self.pal    = 0.5*(tmp[0]**2+tmp[1]**2+tmp[2]**2); tmp=None
####
#### Method to close opened files.
####
   def fin(self):
      """
         close the run files.
      """
      for i in self.vars2l:
         self.__dict__[i+'f'].close()
####
#### Load energies for the run
####
   def loadenergies(self):
      """
         Loads the energies for the run object along with four different
         time series: 
         self.t -> Time in cyclotron units 
         self.tnl -> Time in units of Nominal nonlinear time 
                     based on initial energy
         self.ltnl -> Local Nonlinear time throughout the simulation
         self.ta -> Turbulence age based on the local nonlinear time.
      """
      self.evars=['t', 'eges', 'ebx' , 'eby' , 'ebz' , 'eex' , 'eey' , 'eez' ,\
      'eem' , 'ekix', 'ekiy', 'ekiz', 'ekex', 'ekey', 'ekez', 'ekin', 'eifx', \
      'eify', 'eifz', 'eefx', 'eefy', 'eefz', 'eipx', 'eipy', 'eipz', 'eepx', \
      'eepy', 'eepz']
      data=np.loadtxt(self.rundir+'/Energies.dat')
      for i in self.evars:
         self.__dict__[i]=data[:,self.evars.index(i)]
      self.eb0=0.5*(self.b0x**2+self.b0y**2+self.b0z**2)
      self.eb =self.ebx +self.eby +self.ebz
      self.eip=self.eipx+self.eipy+self.eipz
      self.eep=self.eepx+self.eepy+self.eepz
      self.eif=self.eifx+self.eify+self.eifz
      self.eef=self.eefx+self.eefy+self.eefz
      self.ee =self.eex +self.eey +self.eez
      self.edz=self.eb-self.eb0+self.eif+self.eef
     #self.tnl=self.t*np.sqrt(2*self.edz[0])*2*np.pi/self.lx
      self.tnl=self.t*np.sqrt(4*self.edz[0])*2*np.pi/self.lx
      self.ltnl=self.lx/(np.sqrt(self.edz)*4*np.pi)
      self.ta=np.zeros(len(self.eb))
      for i in range(1,len(self.eb)):
         self.ta[i]=self.ta[i-1]+self.dt*2./(self.ltnl[i-1]+self.ltnl[i])

###
### Method to load parameters
###
def loadparams(paramfile):
   params={}
   def _convert(val):
       constructors = [int, float, str]
       for c in constructors:
           try:
               return c(val)
           except ValueError:
               pass

   with open(paramfile) as f: 
      content = f.readlines()

   for item in content:
      if '#define' in item and item[0] != '!':
         if len(item.split()) > 2:
            key = item.split()[1]
            val = item.split()[2]
            val = _convert(item.split()[2])
         else:
            key = item.split()[1]
            val = True
         params[key] = val

   ## For hybrid code, set electron mass to extremely small 
   ## and speed of light to extremely large.
   if 'hybrid' in params: 
      params['c_2'] = 1e9
      if 'd_e2' in params:
         params['m_e'] = params['d_e2']
      else:
         params['m_e'] = 0.000545
   for i in ['b0x','b0y','b0z']:
      if i not in params:
         params[i]=0.
   
   # Set data type
   if 'eight_byte' in params: params['data_type']='d'
   elif 'four_byte' in params: params['data_type']='f'
   elif 'double_byte' in params: params['data_type']='bb'
   else: params['data_type']='b'
   #
   #
   if params['movie_header'] == '"movie2dC.h"':
      params['logvars']=['rho', 'jx','jy','jz', 'bx','by','bz', 'ex','ey','ez',
      'ne', 'jex','jey','jez', 'pexx','peyy','pezz','pexy','peyz','pexz',
      'ni', 'pixx','piyy','pizz','pixy','piyz','pixz']
   elif params['movie_header'] == '"movie4b.h"':
      params['logvars']=['rho', 'jx','jy','jz', 'bx','by','bz', 'ex','ey','ez',
      'ne','jex','jey','jez', 'pexx','peyy','pezz','pexz','peyz','pexy',
      'ni','jix','jiy','jiz', 'pixx','piyy','pizz' 'pixz','piyz','pixy']
   elif params['movie_header'] == '"movie2dD.h"':
      params['logvars']=['rho', 'jx','jy','jz', 'bx','by','bz', 'ex','ey','ez',
      'ne', 'jex','jey','jez', 'pexx','peyy','pezz','pexy','peyz','pexz',
      'ni', 'jix','jiy','jiz', 'pixx','piyy','pizz','pixy','piyz','pixz']
   elif params['movie_header'] == '"movie3dHeat.h"':
      params['logvars']=['rho', 'jx','jy','jz', 'bx','by','bz', 'ex','ey','ez',
      'ne', 'jex','jey','jez', 'pexx','peyy','pezz','pexy', 'peyz', 'pexz',
      'ni', 'pixx','piyy','pizz','pixy', 'piyz', 'pixz', 'epar1','epar2',
      'epar3','eperp1','eperp2','eperp3', 'vpar1','vpar2','vpar3']
   elif params['movie_header'] == '"movie_pic3.0.h"':
      params['logvars']=['rho', 'jx','jy','jz', 'bx','by','bz', 'ex','ey','ez',
      'ne', 'jex','jey','jez', 'pexx','peyy','pezz','pexy','peyz','pexz',
      'ni', 'jix','jiy','jiz', 'pixx','piyy','pizz','pixy','piyz','pixz']
   else:
      print('='*80 + \
                '\t This particular moive headder has not been coded!\n'\
                '\t Talk to Tulasi to fit it, or fix it yourself.\n'\
                '\t I dont care, Im a computer not a cop'\
                '='*80)
   #
   #
   # Derive some others
   if 'n_movieout' in params:
      params['dtmovie']=params['n_movieout']*params['dt']
   else:
      params['dtmovie']=params['movieout']
   params['nx']=int(params['pex']*params['nx'])
   params['ny']=int(params['pey']*params['ny'])
   params['nz']=int(params['pez']*params['nz'])
   params['dx']=params['lx']/params['nx']
   params['dy']=params['ly']/params['ny']
   params['dz']=params['lz']/params['nz']
   params['xx']=np.linspace(0.,params['lx'],params['nx'])
   params['yy']=np.linspace(0.,params['ly'],params['ny'])
   params['zz']=np.linspace(0.,params['lz'],params['nz'])
   params['xxt']=np.linspace(0.,2*np.pi,params['nx'])
   params['yyt']=np.linspace(0.,2*np.pi,params['ny'])
   params['zzt']=np.linspace(0.,2*np.pi,params['nz'])
   if all([i in params for i in ['b0x','b0y','b0z']]):
      params['b0'] =np.sqrt(params['b0x']**2+params['b0y']**2+params['b0z']**2)
   params['betai'] = 2*params['n_0']*params['T_i']/params['b0']**2
   params['betae'] = 2*params['n_0']*params['T_e']/params['b0']**2
   params['nprocs'] = int(params['pex']*params['pey']*params['pez'])
   params['lambdae']=np.sqrt(params['T_e']/(params['n_0']*params['c_2']))
   if params['m_e'] != 0:
      params['wce']=params['b0']/params['m_e']
      params['vthe']=np.sqrt(2*params['T_e']/params['m_e'])
   else:
      params['wce']=1e9
      params['vthe']=1e9
   params['de']=np.sqrt(params['m_e'])
   params['rhoe']=params['vthe']/params['wce']
   params['wpe']=np.sqrt(params['c_2'])/params['de']
   params['wpi']=np.sqrt(params['c_2'])
   params['vthi']=np.sqrt(2*params['T_i'])
   params['wci']=params['b0']
   params['rhoi']=params['vthi']/params['wci']
   params['kgrid'] = max(np.pi/params['dx'], np.pi/params['dy'], np.pi/params['dz'])
   params['ca']=params['b0']/np.sqrt(params['n_0']) 
   
   return params
