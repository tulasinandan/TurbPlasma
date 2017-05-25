#!/usr/bin/env python
import numpy as np
import pyqtgraph as pg
from subs import *
from pyqtgraph.Qt import QtGui, QtCore

rc=create_object()
variables=raw_input("Variables to load, e.g. all, min, bx by bz: ").split()
rc.vars2load(variables)
rcd=rc.__dict__
bs,fs,step=ask_for_steps(rc.numslices)
smooth=raw_input("Smooth data? y/[n] ")
if smooth == 'y':
   from scipy.ndimage import gaussian_filter as gf
   numsmooth=int(raw_input("How much? "))
cmp=raw_input("Cmap? [BuRd], RdBu, TrmBlk, bryw: ")
if cmp == '':
   cmp='BuRd'
pgcmp = getpgcmap(cmp=cmp)

# Create the window
app=QtGui.QApplication([])
win=pg.GraphicsWindow(title="Multiplot-Test")
win.resize(1000,600)
pg.setConfigOptions(antialias=True)

#Create the canvas
pltdict={}; imgdict={}
for j in rc.vars2l:
   idx=rc.vars2l.index(j)
   if idx < 4:
      haha=[]
   elif np.mod(idx,4) == 0:
      win.nextRow()
   exec('p'+j+'=win.addPlot()')
   exec('pltdict["'+j+'"]=p'+j)
   if (rc.ny > 1):
      exec('img'+j+'=pg.ImageItem()')
      exec('imgdict["'+j+'"]=img'+j)
      pltdict[j].addItem(imgdict[j])
#     pltdict[j].setAspectLocked()
win.show()

# Loop for plottting multiple time slices
for it in range(bs,fs,step):
   print 'Reading time slice ', it
   rc.loadslice(it)
# Make plots
   for j in rc.vars2l:
      if (rc.ny==1 and rc.nz==1):
         pltdict[j].plot(rc.xx,rcd[j][:,0,0],clear=True)
      else:
         if smooth == 'y':
            imgdict[j].setImage(gf(rcd[j][:,::-1,0].T,sigma=numsmooth),clear=True,lut=pgcmp)
         else:
            imgdict[j].setImage(rcd[j][:,::-1,0].T,clear=True,lut=pgcmp)
      pltdict[j].setTitle(j+' '+sminmax(rcd[j]))
   pg.QtGui.QApplication.processEvents()
   haha=raw_input('Hello?')
print 'All done!'
rc.fin()
