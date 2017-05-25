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

#from pylab import *
import numpy as np
from commands import getstatusoutput as syscomout
from os.path import expandvars as opexp
   
class p3d(object):
   """p3d object:
         Tulasi Parashar's version to read data from
         P3D github version.

         Created on 08/22/2014
         Last modified on 09/09/2015
   """
   def __init__(self,shelldirname="",data_type="",filenum=""):
# If no rundir specified
      if len(shelldirname) == 0: # Start init prompt
         shelldirname = raw_input('Please enter the rundir: ') 
      self.rundir = opexp(shelldirname)
# If Git version of the code AND filenum not given
      if len(filenum) == 0:
         self.filenum=raw_input("Please enter the file number to load (e.g. 000): ")
      else:
         self.filenum=filenum
# If data type not specified
      if len(data_type) == 0:
         self.data_type=raw_input("Data Type? [(b)yte/ double (bb)yte/ (f)our byte/ (d)ouble precision] ")
      else:
         self.data_type = data_type
# Parameters to load
      self.params=['pex','pey','pez','nx','ny','nz','lx','ly','lz',\
      'n_movieout','dt','T_i','T_e','n_0','b0x','b0y','b0z','m_e']
# Read parameters from file
      for i in self.params:
         comm="awk '/^#define/ && / "+i+" / {print $3}' "+self.rundir+"/paramfile"
         if syscomout("grep "+i+" "+self.rundir+"/paramfile")[0] == 0:
            exec('self.'+i+'=float(syscomout(comm)[1])')
         else:
            exec('self.'+i+'=float(0.)')
      self.dtmovie=self.n_movieout*self.dt
# Derive some others
      self.nx=int(self.pex*self.nx); 
      self.ny=int(self.pey*self.ny); 
      self.nz=int(self.pez*self.nz)
      self.dx=self.lx/self.nx; 
      self.dy=self.ly/self.ny; 
      self.dz=self.lz/self.nz
      self.xx=np.linspace(0.,self.lx,self.nx); 
      self.yy=np.linspace(0.,self.ly,self.ny); 
      self.zz=np.linspace(0.,self.lz,self.nz)
      self.xxt=np.linspace(0.,2*np.pi,self.nx); 
      self.yyt=np.linspace(0.,2*np.pi,self.ny); 
      self.zzt=np.linspace(0.,2*np.pi,self.nz)
      self.B0 =np.sqrt(self.b0x**2+self.b0y**2+self.b0z**2)
      self.betai = 2*self.n_0*self.T_i/self.B0
      self.betae = 2*self.n_0*self.T_e/self.B0
      self.nprocs = int(self.pex*self.pey*self.pez)
# If byte or double byte data, open the log file.
      if self.data_type in ("b", "bb"):
         print self.rundir+'/staging/movie.log.'+self.filenum
         self.logfile=open(self.rundir+"/staging/movie.log."+self.filenum,"r")
         self.logvars=['rho', 'jx', 'jy', 'jz', 'bx', 'by', 'bz', 'ex', 'ey'\
         , 'ez', 'ne', 'jex', 'jey', 'jez', 'pexx', 'peyy', 'pezz', 'pexy', \
         'peyz', 'pexz', 'ni', 'jix', 'jiy', 'jiz', 'pixx', 'piyy', 'pizz', \
         'pixy', 'piyz', 'pixz']
         self.szl=np.size(self.logvars)
####
#### Method to print the parameters associated with the run
####
   def print_params(self):
      """
         A quick method to print the parameters and variables attached with
         the p3d run object.
      """
      for i in self.params:
         exec('print i," = ",self.'+i)
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
### Byte data ###############################################################
      if self.data_type == 'b':
         f.seek(timeslice*self.nx*self.ny*self.nz)
         field = np.fromfile(f,dtype='int8',count=self.nx*self.ny*self.nz)
         field = np.reshape(field,(self.nx,self.ny,self.nz),order='F')
         exec('minmax=self.'+v+'minmax[timeslice]')
         field = minmax[0]+(minmax[1]-minmax[0])*field/255.
### Double byte data #######################################################
      elif self.data_type == 'bb':
         f.seek(2*timeslice*self.nx*self.ny*self.nz)
         field = np.fromfile(f,dtype='int16',count=self.nx*self.ny*self.nz)
         field = np.reshape(field,(self.nx,self.ny,self.nz),order='F')
         exec('minmax=self.'+v+'minmax[timeslice]')
         field = minmax[0]+(minmax[1]-minmax[0])*(field+32678.)/65535.
### Four byte (single precision) data ######################################
      elif self.data_type == 'f':
         f.seek(4*timeslice*self.nx*self.ny*self.nz)
         field = np.fromfile(f,dtype='float32',count=self.nx*self.ny*self.nz)
         field = np.reshape(field,(self.nx,self.ny,self.nz),order='F')
### Double precision data ##################################################
      elif self.data_type == 'd':
         f.seek(8*timeslice*self.nx*self.ny*self.nz)
         field = np.fromfile(f,dtype='float64',count=self.nx*self.ny*self.nz)
         field = np.reshape(field,(self.nx,self.ny,self.nz),order='F')
      return field
####
#### Method to define the variables to load, create variables and
#### open corresponding files.
#### 
   def vars2load(self,v2l):
      """
         Define the variables to load, define corresponding numpy arrays &
         open the files
      """
      if len(v2l) == 1:
         if v2l[0] == 'min':
               self.vars2l=['bx','by','bz','jix','jiy','jiz','ni']
         elif v2l[0] == 'all':
               self.vars2l=['bx','by','bz','ex','ey','ez','jix','jiy','jiz','jex','jey'\
               ,'jez','jx','jy','jz','ni','pixx','piyy','pizz','pixy','pixz','piyz'\
               ,'pexx','pexy','pexz','peyy','peyz','pezz','ne','rho']
         else:
               self.vars2l=v2l
      else:
         self.vars2l = v2l
# Create arrays and open files
      for i in self.vars2l:
         exec('self.'+i+'=np.array(('+str(self.nx)+','+str(self.ny)+','+str(self.nz)+'))')
         exec('self.'+i+'f=open("'+self.rundir+'/staging/movie.'+i+'.'+self.filenum+'","rb")')
# If 'b' or 'bb' data type, load minmax for each loaded variable
      if self.data_type in ('b', 'bb'):
         d=np.loadtxt(self.logfile)
         for i in self.vars2l:
            exec('self.'+i+'minmax=d[self.logvars.index("'+i+'")::len(self.logvars),:]')
# Find out the number of slices for the open file
      exec('self.'+self.vars2l[0]+'f.seek(0,2)')
      exec('filesize=self.'+self.vars2l[0]+'f.tell()')
#     if self.data_type == 'b':
#        exec('self.numslices=self.'+self.vars2l[0]+'f.tell()/('+str(self.nx)+'*'+str(self.ny)+'*'+str(self.nz)+')')
#     if self.data_type == 'bb':
#        exec('self.numslices=self.'+self.vars2l[0]+'f.tell()/('+str(self.nx)+'*'+str(self.ny)+'*'+str(self.nz)+'*2)')
#     if self.data_type == 'f':
#        exec('self.numslices=self.'+self.vars2l[0]+'f.tell()/('+str(self.nx)+'*'+str(self.ny)+'*'+str(self.nz)+'*4)')
#     if self.data_type == 'd':
#        exec('self.numslices=self.'+self.vars2l[0]+'f.tell()/('+str(self.nx)+'*'+str(self.ny)+'*'+str(self.nz)+'*8)')
      if self.data_type == 'b':
         self.numslices=filesize/(self.nx*self.ny*self.nz)
      if self.data_type == 'bb':
         self.numslices=filesize/(self.nx*self.ny*self.nz*2)
      if self.data_type == 'f':
         self.numslices=filesize/(self.nx*self.ny*self.nz*4)
      if self.data_type == 'd':
         self.numslices=filesize/(self.nx*self.ny*self.nz*8)
####
#### Method to load time slices for the loaded variables.
####
   def loadslice(self,it):
      """
         Load the variables initialized by self.vars2load()
      """
      if self.data_type in ('b', 'bb'):
         for i in self.vars2l:
            exec('self.'+i+'=self.readslice(self.'+i+'f,'+str(it)+',"'+i+'")')
      else:
         for i in self.vars2l:
            exec('self.'+i+'=self.readslice(self.'+i+'f,'+str(it)+')')
      self.ard={}; self.mmd={}
      for i in self.vars2l:
         exec('self.ard["'+i+'"]=self.'+i)
         exec('self.mmd["'+i+'"]=[self.'+i+'.min(),self.'+i+'.max()]')
      self.time = it*self.dtmovie
####
#### Method to close opened files.
####
   def fin(self):
      """
         close the run files.
      """
      for i in self.vars2l:
         exec('self.'+i+'f.close()')
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
      #data=np.loadtxt('fort.9')
      data=np.loadtxt(self.rundir+'/Energies.dat')
      for i in range(np.size(self.evars)):
         exec('self.'+self.evars[i]+'=data[:,'+str(i)+']')
      self.eb0=0.5*(self.b0x**2+self.b0y**2+self.b0z**2)
      self.eb =self.ebx +self.eby +self.ebz
      self.eip=self.eipx+self.eipy+self.eipz
      self.eep=self.eepx+self.eepy+self.eepz
      self.eif=self.eifx+self.eify+self.eifz
      self.eef=self.eefx+self.eefy+self.eefz
      self.ee =self.eex +self.eey +self.eez
      self.edz=self.eb-self.eb0+self.eif+self.eef
      self.tnl=self.t*np.sqrt(2*self.edz[0])*2*np.pi/self.lx
      self.ltnl=self.lx/(np.sqrt(self.edz)*4*np.pi)
      self.ta=np.zeros(len(self.eb))
      for i in range(1,len(self.eb)):
         self.ta[i]=self.ta[i-1]+self.dt*2./(self.ltnl[i-1]+self.ltnl[i])
#     ## Compute Qi, Qe time series
#     dt=self.t[11]-self.t[10]
#     import AnalysisFunctions as af
#     zi,tt,eips=af.fitpoly(rc.t,rc.eip,5,0,rc.t.max(),numpts=500)
#     ze,tt,eeps=af.fitpoly(rc.t,rc.eep,5,0,rc.t.max(),numpts=500)
#     self.qi=np.zeros(len(self.t))
#     self.qe=np.zeros(len(self.t))
#     for i in range(1,len(self.t)-1):
#        self.qi[i]=(eips[i+1]-eips[i-1])/dt
#        self.qe[i]=(eeps[i+1]-eeps[i-1])/dt
