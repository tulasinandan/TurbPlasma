import os
import sys
import numpy as np
from os.path import basename, realpath
from scipy.ndimage import gaussian_filter as gf
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.ticker import AutoMinorLocator
from mpl_toolkits.axes_grid1 import make_axes_locatable
sys.path.insert(0,os.environ['HOME']+'/AJGAR/TurbPlasma')
#
# THE FOLLOWING ROUTINE READS A GENERIC 3D ARRAY FROM A
# DIRECT ACCESS UNFORMATTED FILE AT A GIVEN SNAPSHOT
#
def readsu(filename,nx,ny,nz,timeslice):
   ff=open(filename,mode='rb')
   ff.seek(8*timeslice*nx*ny*nz)
   field = np.fromfile(ff,dtype='float64',count=nx*ny*nz).reshape(field,(nx,ny,nz),order='F')
   ff.close()
   return field

def minmax(a):
   return a.min(),a.max()

def sminmax(a,prec=3):
   exec('min="%.'+str(prec)+'f"%a.min()')
   exec('max="%.'+str(prec)+'f"%a.max()')
   return '('+min+','+max+')' 

def create_object():
   try:
      obj=sys.argv[1]
   except:
      obj=input("At least tell me what kind of code: ")

   exec('from '+obj+' import '+obj)
   exec('myobj = '+obj)
   rc = myobj(*sys.argv[2:])
   return rc

def ask_for_steps(nslices):
   print('There are '+str(nslices)+' time snapshots')
   print('Please provide the first, last and step of snapshots to read:')
   slc=input().split()
   bs=int(slc[0]); fs=int(slc[1]); step=int(slc[2])
   return bs,fs,step


def getpgcmap(cmp='BuRd'):
   import pyqtgraph as pg
   import numpy as np
   if cmp =='bryw':
      STEPS = np.array([0.0, 0.2, 0.6, 1.0])
      CLRS =           ['k', 'r', 'y', 'w']
      ## Create a ColorMap
      clrmp = pg.ColorMap(STEPS, np.array([pg.colorTuple(pg.Color(c)) for c in CLRS]))
      ## Get the LookupTable
      lut = clrmp.getLookupTable()
   elif cmp == 'TrmBlk':
      pos = np.array([0.0, 0.5, 1.0])
      color = np.array([[0,0,0,255], [255,128,0,255], [255,255,0,255]], dtype=np.ubyte)
      map = pg.ColorMap(pos, color)
      lut = map.getLookupTable(0.0, 1.0, 256)
   elif cmp == 'RdBu':
      pos = np.array([0.0,0.5,1.0])
      color = np.array([[255,0,0,0],[255,255,255,255],[0,0,255,0]],dtype=np.ubyte)
      map = pg.ColorMap(pos,color)
      lut = map.getLookupTable(0.0,1.0,256)
   elif cmp == 'BuRd':
      pos = np.array([0.0,0.5,1.0])
      color = np.array([[0,0,255,0],[255,255,255,255],[255,0,0,0]],dtype=np.ubyte)
      map = pg.ColorMap(pos,color)
      lut = map.getLookupTable(0.0,1.0,256)
   return lut

def compute2didx(extar,slc):
   x1=np.argmin(np.abs(extar[0]-slc[0]))
   x2=np.argmin(np.abs(extar[0]-slc[1]))
   y1=np.argmin(np.abs(extar[1]-slc[2]))
   y2=np.argmin(np.abs(extar[1]-slc[3]))
   if len(extar) == 3:
   	z1=np.argmin(np.abs(extar[2]-slc[4]))
   	z2=np.argmin(np.abs(extar[2]-slc[5]))
   	if x1==x2: IDX=np.s_[x1,y1:y2,z1:z2]
   	elif y1==y2: IDX=np.s_[x1:x2,y1,z1:z2]
   	elif z1==z2: IDX=np.s_[x1:x2,y1:y2,z1]
   else:
        IDX=np.s_[x1:x2,y1:y2]
   return IDX

def compute1didx(extar,slc):
   x1=np.argmin(np.abs(extar[0]-slc[0]))
   x2=np.argmin(np.abs(extar[0]-slc[1]))
   if len(extar) == 2:
      y1=np.argmin(np.abs(extar[1]-slc[2]))
      y2=np.argmin(np.abs(extar[1]-slc[3]))
      if x1==x2: 
         IDX=np.s_[x1,y1:y2]
      elif y1==y2: 
         IDX=np.s_[x1:x2,y1]
   if len(extar) == 3:
      z1=np.argmin(np.abs(extar[2]-slc[4]))
      z2=np.argmin(np.abs(extar[2]-slc[5]))
      if (x1==x2 and y1==y2): IDX=np.s_[x1,y1,z1:z2]
      if (y1==y2 and z1==z2): IDX=np.s_[x1:x2,y1,z1]
      if (x1==x2 and z1==z2): IDX=np.s_[x1,y1:y2,z1]
   else:
      IDX=np.s_[x1:x2]
   return IDX

def imss(fdic,
        key,
        cut=None,
        ax=None,
        extent=None,
        cbar=None,
        smth=None,
        nolabels=None,
        lblsz='small',
        **kwargs):
    """
    A wrapper function for imshow to do most 
    tedious stuff for my simulations
    """
    old_ax = plt.gca() # Get Current Axis
    if ax is None: 
        ax = old_ax
    else:
        plt.sca(ax)    # Set Current Axis

    if type(key) is str: plt_val = fdic[key]
    else               : plt_val = key

    if smth is not None:
      plt_val = gf(plt_val,sigma=smth) 

    if cut is None:
      if len(plt_val.shape) == 2:
         IDX=np.s_[:,:]
      else:
         IDX=np.s_[:,:,0]
    else:
      if len(plt_val.shape) == 2:
         IDX=compute2didx([fdic['xx'],fdic['yy']],cut)
      else:
         IDX=compute2didx([fdic['xx'],fdic['yy'],fdic['zz']],cut)

# Use the dict values of xx and yy to set extent
    ext = [fdic['xx'][IDX[0]][0],
           fdic['xx'][IDX[0]][-1],
           fdic['yy'][IDX[1]][0],
           fdic['yy'][IDX[1]][-1]]

    if 'cmap' in kwargs: cmap=kwargs.pop('cmap')
    else:                      cmap='PuOr'

    im = ax.imshow(plt_val[IDX].T,
                           origin='low',
                           extent=ext,
                           cmap=cmap,            # I just love this color map
                           aspect='equal',
                           **kwargs)

    if extent is not None:
        ax.set_xlim(extent[:2])
        ax.set_ylim(extent[2:])

    ax.autoscale(False)

    if nolabels is None:
      ax.set_xlabel(r'$X (d_i)$',size=lblsz)
      ax.set_ylabel(r'$Y (d_i)$',size=lblsz)

    ax.xaxis.set_tick_params(which='both',labelsize=lblsz)
    #minorLocator = AutoMinorLocator()           # Note the second call is so that the minor x ticks are not
    #ax.xaxis.set_minor_locator(minorLocator)    # the same as the y ticks

    ax.yaxis.set_tick_params(which='both',labelsize=lblsz)
    #minorLocator = AutoMinorLocator()
    #ax.yaxis.set_minor_locator(minorLocator)

    plt.minorticks_on()
    plt.sca(old_ax)

    # Code to implement for a cbar
    if cbar:
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", "3%", pad="1.5%")
        plt.colorbar(im, cax=cax)

        cax.xaxis.set_tick_params(which='both',labelsize=lblsz)
        cax.yaxis.set_tick_params(which='both',labelsize=lblsz)

        plt.draw()
        return im,cax

    else:
        return im

def iplt(fdic,
        key,
        cut=None,
        ax=None,
        cbar=None,
        smth=None,
        nolabels=None,
        lblsz='small',
        **kwargs):
    """
    A wrapper function for imshow to do most 
    tedious stuff for my simulations
    """
    old_ax = plt.gca() # Get Current Axis
    if ax is None: 
        ax = old_ax
    else:
        plt.sca(ax)    # Set Current Axis

    if type(key) is str: plt_val = fdic[key]
    else               : plt_val = key

    if cut is None:
      if len(plt_val.shape) == 1:
         IDX=np.s_[:]
      elif len(plt_val.shape) == 2:
         IDX=np.s_[:,0]
      else:
         IDX=np.s_[:,0,0]
    else:
      if len(plt_val.shape) == 1:
         IDX=compute1didx([fdic['xx']],cut)
      elif len(plt_val.shape) == 2:
         IDX=compute1didx([fdic['xx'],fdic['yy']],cut)
      else:
         IDX=compute1didx([fdic['xx'],fdic['yy'],fdic['zz']],cut)

    im = ax.plot(fdic['xx'],plt_val[IDX], **kwargs)

    ax.autoscale(False)

    if nolabels is None:
      ax.set_xlabel(r'$X (d_i)$',size=lblsz)
      ax.set_ylabel(key,size=lblsz)

    ax.xaxis.set_tick_params(which='both',labelsize=lblsz)
    ax.yaxis.set_tick_params(which='both',labelsize=lblsz)

    plt.minorticks_on()
    plt.sca(old_ax)

    return im

def calc_dep(v2lu,primitives,derived):
   # Find out the dependencies and append them before the variable
   lgcl=[False]*len(v2lu)
   # If there is even a single False, go through the loop
   while not all(lgcl):
   # For all values in the list
      for i in v2lu:
   # Find the index of the item in list
         idx=v2lu.index(i)
   # If it is a primitive item, set the dependency satisfied flag to true
         if i in primitives:
            lgcl[idx]=True
   # If it is a derived quantity, check all its dependencied for existence in the list
         elif i in derived:
   # create empty logical and insertion lists to insert dependencies before item i
            tmplgcl=[]; insrt=[]
            for j in derived[i]:
   # If the dependency of the derived quantity is in list, don't do anything
               if j in v2lu:
                  pass
               else:
   # Else append the insertion list
                  insrt.append(j)
   # If the inserted quantity is a primitive, set the dependency satisfied flag to true, else false
                  if j in primitives:
                     tmplgcl.append(True)
                  else:
                     tmplgcl.append(False)
   # Insert the newly found quantities into the list _before_ i
            v2lu=v2lu[:idx]+insrt+v2lu[idx:]
   # Set the dependency satisfied flag for i to be True
            lgcl[idx]=True
   # Insert the dependency satisfied flags for the newly inserted quantities
            lgcl=lgcl[:idx]+tmplgcl+lgcl[idx:]
   return v2lu
