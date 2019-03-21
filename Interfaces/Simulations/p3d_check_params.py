#!/usr/bin/env python
import numpy as np
from TurbPlasma.Interfaces.Simulations.p3d import loadparams
from TurbPlasma.Analysis.Simulations import tfps
from sys import argv as argv

def printwarn(something="SOMETHING"):
   print() 
   print()
   print()
   print('!!!!!!!!!!!!!!!! WARNING !!!!!!!!!!!!!!!!!')
   print('!!!!!!!!!!!!!!!! WARNING !!!!!!!!!!!!!!!!!')
   print(something,' UNDER-RESOLVED')
   print()
   print()
   print()

def check_params(paramfile):
    dd=loadparams(paramfile)
    
    if 'nu' not in dd: dd['nu'] = 0.
    if 'eta' not in dd: dd['eta'] = 0.
    if 'chi' not in dd: dd['chi'] = 0.
    tn=3.e-07
    
    #does time step sample all the frequencies appropriately?
    print()
    print()
    print('############ NUMERICAL PARAMETERS ###############')
    print('Lx=','%.3f'%dd['lx'],'\t Ly=','%.3f'%dd['ly'],'\t Lz=','%.3f'%dd['lz'],'\t ppg=',dd['ppg'])
    print('Nx=',dd['nx'],'\t Ny=', dd['ny'],'\t Nz=',dd['nz'],'\t On',dd['nprocs'],' procs')
    print('b0=','%.3f'%dd['b0'],'\t dx=','%.3f'%dd['dx'],'\t dy=','%.3f'%dd['dy'],'\t dz=','%.3f'%dd['dz'])
    print('me/mi =','%.5f'%dd['m_e'],'\t nu=','%.5f'%dd['nu'],'\t eta=','%.5f'%dd['eta'],'\t kgrid=','%.3f'%dd['kgrid'])
    print('dt=','%.6f'%dd['dt'],'\t substeps=',dd['substeps'])
    print('Expected Runtime is: \n(xNsteps for time to request):',tn*dd['nx']*dd['ny']*dd['nz']*2*dd['ppg']/(dd['nprocs']*60.),'minutes/(simulation step)')
    print()
    print('############# PLASMA PARAMETERS #################')
    print('beta_i=','%.3f'%dd['betai'],'\t betae=','%.3f'%dd['betae'],'\t me/mi=','%.5f'%dd['m_e'])
    print('rhoi=','%.3f'%dd['rhoi'],'\t rhoe=','%.3f'%dd['rhoe'],'\t Debey Length=','%.3f'%dd['lambdae'])
    print('Va=','%.3f'%dd['ca'],'\t Vthi=','%.3f'%dd['vthi'],'\t Vthe=','%.3f'%dd['vthe'],'\t c=','%.3f'%np.sqrt(dd['c_2']))
    print('wci=','%.3f'%dd['wci'],'\t wce=','%.3f'%dd['wce'], 'wpi=','%.3f'%dd['wpi'],'\t wpe=','%.3f'%dd['wpe'])
    
    print()
    print('######### RESOLUTION ESTIMATES #################')
    print() 
    print('##### TWO FLUID WAVES ==>')
    f,s = tfps(beta=dd['betai']+dd['betae'], ca=dd['ca'], de2=dd['m_e'], theta=0., kk=dd['kgrid'])
    tmp = dd['dt']*max(s)/(dd['substeps']*min([dd['dx'],dd['dy'],dd['dz']]))
    tmpnu=dd['nu']*dd['kgrid']**3/max(s)
    print('At',0.*90,': F/W=','%.5f'%s[0],' A/KAW=','%.5f'%s[1],' S/A=','%.5f'%s[2])
    print('Courant factor:',' %.3f'%tmp,', nu*k^3/Cf=','%.5f'%tmpnu)
    if tmp >= 0.9: printwarn("PARALLEL TWO FLUID WAVES")
    print('--')
    
    f,s = tfps(beta=dd['betai']+dd['betae'], ca=dd['ca'], de2=dd['m_e'], theta=45, kk=dd['kgrid'])
    tmp = dd['dt']*max(s)/(dd['substeps']*min([dd['dx'],dd['dy'],dd['dz']]))
    tmpnu=dd['nu']*dd['kgrid']**3/max(s)
    print('At',0.5*90.,': F/W=','%.5f'%s[0],' A/KAW=','%.5f'%s[1],' S/A=','%.5f'%s[2])
    print('Courant factor:',' %.3f'%tmp,', nu*k^3/Cf=','%.5f'%tmpnu)
    if tmp >= 0.9: printwarn("OBLIQUE TWO FLUID WAVES")
    print('--')
    
    f,s = tfps(beta=dd['betai']+dd['betae'], ca=dd['ca'], de2=dd['m_e'], theta=80, kk=dd['kgrid'])
    tmp = dd['dt']*max(s)/(dd['substeps']*min([dd['dx'],dd['dy'],dd['dz']]))
    tmpnu=dd['nu']*dd['kgrid']**3/max(s)
    print('At',0.88*90.,': F/W=','%.5f'%s[0],' A/KAW=','%.5f'%s[1],' S/A=','%.5f'%s[2])
    print('Courant factor:',' %.3f'%tmp,', nu*k^3/Cf=','%.5f'%tmpnu)
    if tmp >= 0.9: printwarn("HIGHLY OBLIQUE TWO FLUID WAVES")
    print('--')
    
    print()
    print('##### MISC SPEEDS ==>')
    #
    tmp = dd['dt']*dd['ca']/(2*dd['substeps']*min([dd['dx'],dd['dy'],dd['dz']]))
    if tmp >= 0.9: printwarn("ALFVEN WAVE")
    print('Courant factor for Alfven waves is: ','%.3f'%tmp)
    #
    tmp=dd['vthi']*dd['dt']/(min([dd['dx'],dd['dy'],dd['dz']]))
    if tmp >= 0.9: printwarn("PROTON THERMAL SPEED")
    print('Courant factor for proton thermal speed is: ','%.3f'%tmp)
    #
    if 'hybrid' in dd:
       print()
    else:
    #
       tmp=dd['vthe']*dd['dt']/(min([dd['dx'],dd['dy'],dd['dz']]))
       if tmp >= 0.9: printwarn("ELECTRON THERMAL SPEED")
       print('Courant factor for electron thermal speed is: ','%.3f'%tmp)
    #
       tmp=np.sqrt(dd['c_2'])*dd['dt']/(2*dd['substeps']*dd['dx'])
       if tmp >= 0.9: printwarn("SPEED OF LIGHT")
       print('Courant factor for speed of light is: ','%.3f'%tmp)
    print()
    print('##### TIMES ==>')
    #
    tmp=1./(dd['dt']*dd['wci']) 
    if tmp < 30: printwarn("PROTON CYCLOTRON TIME")
    print('%.3f'%tmp,' time steps /wci')
    #
    if 'hybrid' in dd:
       print()
    else:
    #
       tmp=1./(dd['dt']*dd['wce'])
       if tmp < 10: printwarn("ELECTRON CYCLOTRON TIME")
       print('%.3f'%tmp,' time steps /wce')
    #
       tmp=1./(dd['dt']*dd['wpi'])
       if tmp < 20: printwarn("PROTON PLASMA TIME")
       print('%.3f'%tmp,' time steps /wpi')
    #
       tmp=1./(dd['dt']*dd['wpe']) 
       if tmp < 10: printwarn("ELECTRON PLASMA TIME")
       print('%.3f'%tmp,' time steps /wpe')
    print()
    print('##### LENGTHS ==>')
    #
    tmp=1/min([dd['dx'],dd['dy'],dd['dz']])
    if tmp < 6 : printwarn("PROTON INERTIAL LENGTH")
    print('At max ','%.3f'%tmp,' grid points /di')
    #
    tmp=dd['rhoi']/min([dd['dx'],dd['dy'],dd['dz']])
    if tmp < 10: printwarn("PROTON GYRO LENGTH")
    print('At max ','%.3f'%tmp,' grid points /rhoi')
    #
    if 'hybrid' in dd:
       print()
    else:
    #
       tmp=dd['lambdae']/min([dd['dx'],dd['dy'],dd['dz']])
       if tmp < 1./np.pi: printwarn("DEBYE LENGTH")
       print('At max ','%.3f'%tmp,' grid points /lde')
    #
       tmp=dd['rhoe']/min([dd['dx'],dd['dy'],dd['dz']])
       if tmp < 6 : printwarn("ELECTRON GYRO LENGTH")
       print('At max ','%.3f'%tmp,' grid points /rhoe')
    #
       tmp=dd['de']/min([dd['dx'],dd['dy'],dd['dz']])
       if tmp < 6 : printwarn("ELECTRON INERTIAL LENGTH")
       print('At max ','%.3f'%tmp,' grid points /de')
    print()

if __name__=="__main__":
    try:
       argv[1]
    except NameError:
       print("well, what paramfile do you want to use?")
       paramfile=raw_input("Param file: ")
    else:
       paramfile=argv[1]

    check_params(paramfile)
