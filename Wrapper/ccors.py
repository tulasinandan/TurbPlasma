import numpy as np
import matplotlib.pyplot as plt
import AnalysisFunctions as af
lddata=input('Load Data again? ')
if lddata == 'y':
   from subs import create_object
   rc=create_object()
   rc.vars2load(['all'])
   rcd=rc.__dict__
   rc.loadenergies()
   nltime=np.float(input("What nonlinear time? "))
   slc=np.int(rc.t[np.argmin(np.abs(rc.ta-nltime))]/rc.dtmovie)
   smt=float(input("How smooth? "))
   rc.loadslice(slc,smth=smt)
   rc.computevars(['tempi','omi','tempe','ome'])
   rc.addattr(['Ap','dti','tite','dte'],[0.5*(rc.tix+rc.tiy)/rc.tiz,rc.ti-np.mean(rc.ti),rc.ti*rc.T_e/(rc.te*rc.T_i),rc.te-np.mean(rc.te)])
   rc.addattr(['dAp','dtite'],[rc.Ap-1.,rc.tite-1.0])
   gj=af.pgrad(rc.jz,rc.dx,rc.dy,rc.dz,smth=1)
   bdgj = rc.bx*gj[0]+rc.by*gj[1]+rc.bz*gj[2]

def correlation(ar1,ar2,ax=0,dx=1.):
   nlen=np.shape(ar1)[ax]/2
   r=np.zeros(nlen);corr=np.zeros(nlen)
   ar1=ar1-np.mean(ar1)
   ar2=ar2-np.mean(ar2)
   mamb=np.mean(ar1)*np.mean(ar2)
   meanab=np.mean(ar1*ar2)
   rmsab =np.sqrt(np.mean(ar1**2)*np.mean(ar2**2))
   for i in range(nlen):
      ars=np.roll(ar2,i,axis=ax)
      corr[i]=np.mean(ar1*ars)#-mamb
      r[i]=i*dx
   corr = corr/rmsab
   return r,corr

def plot_cor(a,b,nama,namb):
   import matplotlib.pyplot as plt
   rs,cors=correlation(a,b)
   ru,coru=correlation(np.abs(a),np.abs(b))
   plt.clf()
   plt.plot(rs,cors/rmsabs,label='signed')
   plt.plot(ru,coru/rmsabu,label='unsigned')
   plt.xlabel('$\Delta_x$')
   plt.title('Correlation between '+nama+' and '+namb)
   plt.legend(loc='best')
   plt.show()

#lst = ['omzi','jz','dti','dAp','dtite','omze','dte']
#lst = ['omzi','jz','dti','dAp']
lst=input("What variables to compute correlations for? ").split()
from itertools import combinations 
l=list(combinations(lst,2)); sz=len(l)

d = {}
for i in range(sz):
   a = l[i][0]; b = l[i][1]
   print('Computing correlations for',a,b, '...')
   rs,d['c-'+a+'s-'+b+'s'] = correlation(rcd[a],rcd[b])
   rs,d['c-'+a+'s-'+b+'u'] = correlation(rcd[a],np.abs(rcd[b]))
   rs,d['c-'+a+'u-'+b+'s'] = correlation(np.abs(rcd[a]),rcd[b])
   rs,d['c-'+a+'u-'+b+'u'] = correlation(np.abs(rcd[a]),np.abs(rcd[b]))

d['r']=rs*rc.dx

np.save('cors.npy',d)
