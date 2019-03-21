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

def multigray(rc,variables,bs,fs,step,movieout,cmp):
    import matplotlib.pylot as plt
    import numpy as np
    from TurbPlasma.Utilities.subs import iplt,imss
    rc.vars2load(variables)
    rcd=rc.__dict__
    
    # find out the screen layout
    # plots will be shown on a nplotx,nplotx grid
    if np.modf(np.sqrt(len(rc.vars2l)))[0] == 0:
       nplotx=int(np.modf(np.sqrt(len(rc.vars2l)))[1])
    else: 
       nplotx=int(np.modf(sqrt(len(rc.vars2l)))[1]+1)
    f,axarr=plt.subplots(nplotx,nplotx,sharex='col',sharey='row')
    figtitle=f.suptitle('Multiplots')
    
    # Loop for plottting multiple time slices
    for it in range(bs,fs,step):
       print('Reading time slice ', it)
       rc.loadslice(it)
    
       for j in rc.vars2l:
          l,m=rc.vars2l.index(j)/nplotx,np.mod(rc.vars2l.index(j),nplotx)
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
          plt.draw()
          haha=input('Hello?')
       elif movieout == 'm':
          imagefilename='movie.%04d.'% it
          plt.savefig(imagefilename+rc.dirname+'.png',bbox_inches='tight')
       else:
          print('Please make a valid choice')

if __name__=="__main__":
    from TurbPlasma.Utilities.subs import create_object
    rc = create_object()
    variables=input("Variables to load, e.g. all, min, bx by bz: ").split()
    bs,fs,step = ask_for_steps(rc.numslices)
    movieout = input('Output (m)ovies or (s)how plots? [default: s]: ')
    if movieout == '':
       movieout = 's'
    cmp = input("cmap? e.g. gray,seismic,RdBu,jet,hot,spectral,hot_r etc. [default RdBu]: ")
    if cmp == '':
       cmp = 'RdBu'
    multigray(rc,variables,bs,fs,step,movieout,cmp)
    rc.fin()
