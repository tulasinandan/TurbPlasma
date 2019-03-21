#!/usr/bin/env python
import matplotlib.pyplot as plt

def eplt(rc,timeunit='wci',sv=None,**kwargs):
   rc.loadenergies(); rcd=rc.__dict__

   if timeunit == 'wci':
      time = rc.t; timelabel='$t \omega_{ci}$'
   elif timeunit == 'wpi':
      time = rc.t*rc.wpi; timelabel='$t \omega_{pi}$'
   elif timeunit == 'wpe':
      time = rc.t*rc.wpe; timelabel='$t \omega_{pe}$'
   elif timeunit == 'nl':
      time=rc.ta; timelabel='$\hat{t}_{nl}$'
   
   pltvars=['eges','eb','eip','eep','eif','eef','ee','ev']
   pltlbls=['$\delta E_{tot}$','$\delta E_{B}$','$\delta E_{th_i}$',\
            '$\delta E_{th_e}$','$\delta E_{if}$','$\delta E_{ef}$',\
            '$\delta E_{E}$','$\delta E_{v}$']
   pltclrs=['b','k','g','m','c','y','#800000','r']
   
   for i in pltvars:
      if i in rcd:
         idx=pltvars.index(i)
         plt.plot(time,rcd[i]-rcd[i][0],label=pltlbls[idx],color=pltclrs[idx],**kwargs)
   plt.legend(loc='best',ncol=3,fancybox=True,framealpha=0.2)
   plt.axis('tight')
   plt.xlabel(timelabel,size='x-large')
   plt.ylabel('$\delta E$',size='x-large')
   plt.title('Change in energy for '+rc.dirname+'($\\beta_p=$'+str(rc.betai)+',$\\beta_e$='+str(rc.betae)+',$L_x$='+str(rc.lx)+')')
   if sv:
      plt.savefig('eplots.'+rc.dirname+'.'+sv,dpi=300, bbox_inches='tight')
   plt.show()

def turbeplt(rc,sv=None,**kwargs):
   rc.loadenergies()
   plt.clf()
   plt.plot(rc.tnl,rc.eb+rc.ee+rc.eif+rc.eef-rc.eb0,label='Fluctuation')
   plt.plot(rc.tnl,rc.eep+rc.eip-rc.eep[0]-rc.eip[0],label='Thermal')
   plt.plot(rc.tnl,rc.eb+rc.ee+rc.eif+rc.eef+rc.eip+rc.eep-rc.eb0-rc.eep[0]-rc.eip[0],label='Total')
   plt.plot(rc.tnl,rc.eb-rc.eb0,alpha=0.3,label='Eb')
   plt.plot(rc.tnl,rc.eif+rc.eef,alpha=0.3,label='Flow')
   plt.legend(loc='best',fancybox=True,framealpha=0.2)
   plt.xlabel(r'$\tau_{nl}$')
   plt.ylabel(r'$E(\tau)$')
   if sv is not None:
      plt.savefig('Turb-Energetics'+rc.dirname+'.'+sv,dpi=300,bbox_inches='tight')
   plt.show()

if __name__=="__main__":
   from TurbPlasma.Utilities.subs import create_object
   rc=create_object()
   eplt(rc,timeunit='nl',sv='png',linewidth=2)
