#!/usr/bin/env python
from pylab import *
from p3d import p3d
from os.path import basename
if len(sys.argv) == 1:
   rdir=raw_input("Physics directory, thermal directory, code_type, data_type: ")
else:
   rdir=sys.argv[1]
foo = []
for j in range(4):
   if  j < len(sys.argv)-1:
      foo.append(sys.argv[j+1])
   else:
      foo.append('')

rc = p3d(foo[0],foo[2],foo[3])
rc.loadenergies()
rt = p3d(foo[1],foo[2],foo[3])
rt.loadenergies()

evars=rc.evars+['eip','eep']
for i in evars[1:]:
#For example  ee  =rc.ee  -rt.ee + rt.ee[0]
   exec('rc.'+i+'=rc.'+i+'-rt.'+i+'+rt.'+i+'[0]')

plot(rc.t,rc.eges-rc.eges[0], label='$\delta E_{tot}$')
plot(rc.t,rc.eb  -rc.eb[0]  , label='$\delta E_{b}$')
plot(rc.t,rc.eip -rc.eip[0] , label='$\delta E_{ip}$')
plot(rc.t,rc.eep -rc.eep[0] , label='$\delta E_{ep}$')
plot(rc.t,rc.eif -rc.eif[0] , label='$\delta E_{if}$')
plot(rc.t,rc.eef -rc.eef[0] , label='$\delta E_{ef}$')
plot(rc.t,rc.ee  -rc.ee[0]  , label='$\delta E_{e}$')
legend(loc='best',ncol=3,fancybox=True,framealpha=0.2)
#legend(loc='best', fancybox=True, framealpha=0.5)
axis('tight')
xlabel('$\hat{t}_{nl}$',size='x-large')
ylabel('$\delta E$',size='x-large')
title('Thermally Corrected change in energy for '+basename(rdir)+'($\\beta_p=$'+str(rc.betai)+',$\\beta_e$='+str(rc.betae)+',$L_x$='+str(rc.lx)+')',size='small')
savefig('TC-eplots.'+basename(rdir)+'.png',dpi=300, bbox_inches='tight')
show()
