#!/usr/bin/env python
#
# Code to show multiple variable time slices on one plot
# Equivalent of multigray.pro from the historic P3D IDL 
# codes. It uses the p3d object to load the run data.
#
# Right now, it works best from the ipython shell. It has
# some issues if run directly from the command line. I'll
# maybe fix that one day but until then, use it from the
# ipython shell.
#
#                 Tulasi Nandan Parashar
#                 2015/09/09
########################################################
from pylab import *
from subs import *

rc = create_object()

variables=raw_input("Variables to load, e.g. all, min, bx by bz: ").split()
rc.vars2load(variables)
rcd=rc.__dict__

movieout = raw_input('Output (m)ovies or (s)how plots? [default: s]: ')
if movieout == '':
   movieout = 's'
cmp = raw_input("cmap? e.g. gray,seismic,RdBu,jet,hot,spectral,hot_r etc. [default RdBu]: ")
if cmp == '':
   cmp = 'RdBu'

bs,fs,step = ask_for_steps(rc.numslices)

# find out the screen layout
# plots will be shown on a nplotx,nplotx grid
if np.modf(sqrt(len(rc.vars2l)))[0] == 0:
   nplotx=int(np.modf(sqrt(len(rc.vars2l)))[1])
else: 
   nplotx=int(np.modf(sqrt(len(rc.vars2l)))[1]+1)
f,axarr=subplots(nplotx,nplotx,sharex='col',sharey='row')
figtitle=f.suptitle('Multiplots')

# Loop for plottting multiple time slices
for it in range(bs,fs,step):
   print 'Reading time slice ', it
   rc.loadslice(it)

   for j in rc.vars2l:
      l,m=rc.vars2l.index(j)/nplotx,mod(rc.vars2l.index(j),nplotx)
      axarr[l,m].clear()
      if (rc.ny==1 and rc.nz==1):
        #axarr[l,m].plot(rc.xx,rcd[j][:,0,0])
         iplt(rcd,j,ax=axarr[l,m])
      else:
         imss(rcd,j,interpolation="none",ax=axarr[l,m],nolabels=1)
         if l < nplotx-1:
            setp(axarr[l,m].get_xticklabels(),visible=False)
         if m > 0:
            setp(axarr[l,m].get_yticklabels(),visible=False)
      axarr[l,m].set_title(j+' '+sminmax(rcd[j]),size="x-small")
      axarr[l,m].set_adjustable("box-forced")
   figtitle.set_text('t='+"%.3f"%rc.time)

   if movieout == 's':
      f.show()
      draw()
      haha=raw_input('Hello?')
   elif movieout == 'm':
      imagefilename='movie.%04d.'% it
      savefig(imagefilename+rc.dirname+'.png',bbox_inches='tight')
   else:
      print 'Please make a valid choice'
rc.fin()
