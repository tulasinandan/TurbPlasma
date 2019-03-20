#!/home/tulasi/bin/python
###########################################################################
##
## The class to handle f3d simulation data in python. An object is created 
## that has the run parameters attached with it. The basic usage of the 
##  class can be done in the following ways.
##
## Example 1:
## > from f3d import f3d
## > r1=f3d('OTV')
## > r1.vars2load(['bx','by','bz'])
## > for i in range(20):
## >     r1.loadslice(i)
## >     DO WHATEVER!!!
## > r1.fin()
##
##
## Example 2:
## > from f3d import f3d
## > r1=f3d('OTV')
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

from pylab import *
from subprocess import getstatusoutput as syscomout
from os.path import expandvars as opexp
   
class f3d(object):
   """f3d object:
         Tulasi Parashar's version.
         Created on 08/22/2014
         Last modified on 08/23/2014
   """
   def __init__(self,shelldirname="",data_type=""):
# If no rundir specified
      if len(shelldirname) == 0: # Start init prompt
         shelldirname = input('Please enter the rundir: ') 
      self.rundir = opexp(shelldirname)
# If data type not specified
      if len(data_type) == 0:
         self.data_type=input("Data Type? [(b)yte/ (d)ouble precision] ")
      else:
         self.data_type = data_type
# Parameters to load
      self.readparams=['Nx0_tot', 'Ny0_tot', 'Nz0_tot', 'ALX_TOT', 'ALY_TOT',\
      'ALZ_TOT', 'DTP', 'NOUTMOV', 'NPTS', 'SD', 'TEMPF','BX0','BY0','BZ0',\
      'DENLOBE']
      self.params=['nx', 'ny', 'nz', 'lx', 'ly', 'lz', 'dt', 'NOUTMOV', \
      'npts', 'SD', 'Temp','b0x','b0y','b0z','n_0']
# Read parameters from file
      for i in range(len(self.readparams)):
         comm="awk '/^ "+self.readparams[i]+" / {print $3}' "+self.rundir+"/paramfile"
         if syscomout(comm)[1].endswith(','):
            print(self.params[i], '=', float(syscomout(comm)[1][:-1]))
         else:
            print(self.params[i], '=', float(syscomout(comm)[1]))
         if syscomout("grep "+self.readparams[i]+" "+self.rundir+"/paramfile")[0] == 0:
            if syscomout(comm)[1].endswith(','):
               exec('self.'+self.params[i]+'=float(syscomout(comm)[1][:-1])')
            else:
               exec('self.'+self.params[i]+'=float(syscomout(comm)[1])')
         else:
            exec('self.'+self.params[i]+'=float(0.)')
# Derive some others
      self.nx=int(self.nx)
      self.ny=int(self.ny)
      self.nz=int(self.nz)
      self.dx=self.lx/self.nx; 
      self.dy=self.ly/self.ny; 
      self.dz=self.lz/self.nz
      self.xx=np.linspace(0.,self.lx,self.nx); 
      self.yy=np.linspace(0.,self.ly,self.ny); 
      self.zz=np.linspace(0.,self.lz,self.nz)
      self.xxt=np.linspace(0.,2*pi,self.nx); 
      self.yyt=np.linspace(0.,2*pi,self.ny); 
      self.zzt=np.linspace(0.,2*pi,self.nz)
      self.beta = 2*self.n_0*self.Temp/(self.b0x**2+self.b0y**2+self.b0z**2)

####
#### Method to print the parameters associated with the run
####
   def print_params(self):
      """
         A quick method to print the parameters and variables attached with
         the f3d run object.
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
            self.vars2l=['bx','by','bz','curpx','curpy','curpz','den']
            self.arrs2l=['bx','by','bz','jix','jiy','jiz','ni']
         elif v2l[0] == 'all':
            self.vars2l=['bx','by','bz','curpx','curpy','curpz','curx',\
            'cury','curz','pfi','psi','den']
            self.arrs2l=['bx','by','bz','jix','jiy','jiz','jx',\
            'jy','jz','pfi','psi','ni']
#        else:
#              self.vars2l=v2l
#     else:
#        self.vars2l = v2l
# Create arrays and open files
      for i in range(len(self.vars2l)):
         exec('self.'+self.arrs2l[i]+'=np.array(('+str(self.nx)+','+str(self.ny)+','+str(self.nz)+'))')
         exec('self.'+self.vars2l[i]+'f=open("'+self.rundir+'/'+self.vars2l[i]+'f","rb")')
# If 'b' or 'bb' data type, load minmax for each loaded variable
      if self.data_type in ('b', 'bb'):
         d=loadtxt(self.logfile)
         for i in self.vars2l:
            exec('self.'+i+'minmax=d[self.logvars.index("'+i+'")::len(self.logvars),:]')
# Find out the number of slices for the open file
      exec('self.'+self.vars2l[0]+'f.seek(0,2)')
      if self.data_type == 'b':
         exec('self.numslices=self.'+self.vars2l[0]+'f.tell()/('+str(self.nx)+'*'+str(self.ny)+'*'+str(self.nz)+')')
      if self.data_type == 'd':
         exec('self.numslices=self.'+self.vars2l[0]+'f.tell()/('+str(self.nx)+'*'+str(self.ny)+'*'+str(self.nz)+'*8)')
####
#### Method to load time slices for the loaded variables.
####
   def loadslice(self,it):
      """
         Load the variables initialized by self.vars2load()
      """
      if self.data_type in ('b'):
         for i in range(len(self.vars2l)):
            exec('self.'+self.arrs2l[i]+'=self.readslice(self.'+self.vars2l[i]+'f,'+str(it)+',"'+i+'")')
      else:
         for i in range(len(self.vars2l)):
            exec('self.'+self.arrs2l[i]+'=self.readslice(self.'+self.vars2l[i]+'f,'+str(it)+')')
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
      print('NEED TO IMPLEMENT IT')
