#!/usr/bin/env python
import matplotlib.pyplot as plt

def eplt(rc,timeunit='c',sv=None,**kwargs):
   rc.loadenergies(); rcd=rc.__dict__

   pltvars=['eges','eb','eip','eep','eif','eef','ee','ev']
   pltlbls=['$\delta E_{tot}$','$\delta E_{B}$','$\delta E_{th_i}$',\
            '$\delta E_{th_e}$','$\delta E_{if}$','$\delta E_{ef}$',\
            '$\delta E_{E}$','$\delta E_{v}$']
   pltclrs=['b','k','g','m','c','y','#800000','r']
   
   if timeunit == 'c':
      time = rc.t; timelabel='$t \omega_{ci}$'
   elif timeunit == 'nl':
      time=rc.ta; timelabel='$\hat{t}_{nl}$'
   
   for i in pltvars:
      if i in rcd:
         idx=pltvars.index(i)
         plt.plot(time,rcd[i]-rcd[i][0],label=pltlbls[idx],color=pltclrs[idx],**kwargs)
   plt.legend(loc='best',ncol=3,fancybox=True,framealpha=0.2)
   plt.axis('tight')
   plt.xlabel(timelabel,size='x-large')
   plt.ylabel('$\delta E$',size='x-large')
   plt.title('Change in energy for '+rc.dirname+'($\\beta_p=$'+str(rc.betai)+',$\\beta_e$='+str(rc.betae)+',$L_x$='+str(rc.lx)+')')
   if sv is not None:
      plt.savefig('eplots.'+rc.dirname+'.'+sv,dpi=300, bbox_inches='tight')
   plt.show()

if __name__=="__main__":
   from subs import create_object
   rc=create_object()
   eplt(rc,timeunit='nl',sv='png',linewidth=2)
