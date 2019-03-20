import sys
import os
sys.path.insert(0,os.environ['HOME']+'/AJGAR/TurbPlasma')
import pandas as pd
import numpy as np
import AnalysisFunctions as af
from . import Time_Series_Analysis as tsa

def create_df_ascii(dirname,chop_ends=100,hampel_filter=None,lowpass_freq=None,spacecraft=1):
   try:
      print('Creating Data Frame from '+dirname)
   except:
      print("No directory given, exiting")
      return
   spacecraft = str(spacecraft)
   cadence=os.path.basename(dirname)

   if cadence == 'B':
      suffix=spacecraft+'_resB1.dat'
   elif cadence == 'electrons':
      suffix=spacecraft+'_resNe_1.dat'
   elif cadence == 'ions':
      suffix=spacecraft+'_resNi_1.dat'
   elif cadence == 'B_hmp':
      suffix=spacecraft+'_resB1_hmp.dat'
   elif cadence == 'B_hmp_lp':
      suffix=spacecraft+'_resB1_hmp_lp.dat'
   elif cadence == 'electrons_hmp':
      suffix=spacecraft+'_resNe_1_hmp.dat'
   elif cadence == 'electrons_hmp_lp':
      suffix=spacecraft+'_resNe_1_hmp_lp.dat'
   elif cadence == 'ions_hmp':
      suffix=spacecraft+'_resNi_1_hmp.dat'
   elif cadence == 'ions_hmp_lp':
      suffix=spacecraft+'_resNi_1_hmp_lp.dat'
   

   b =np.loadtxt(dirname+'/B'  +suffix) 
   dt=b[2,0]-b[1,0]
   nt=len(b[:,0]); print('Number of entries nt=',nt)
   nt=nt if nt%2==0 else nt-1
#     
   if chop_ends is not None:
      b=b[chop_ends:nt,:][:-chop_ends,:] 
      vi=np.loadtxt(dirname+'/Vi_'+suffix)[chop_ends:nt,:][:-chop_ends,:] 
      ve=np.loadtxt(dirname+'/Ve_'+suffix)[chop_ends:nt,:][:-chop_ends,:] 
      ni=np.loadtxt(dirname+'/Ni_'+suffix)[chop_ends:nt,:][:-chop_ends,:] 
      ne=np.loadtxt(dirname+'/Ne_'+suffix)[chop_ends:nt,:][:-chop_ends,:] 
      R =np.loadtxt(dirname+'/R'  +suffix)[chop_ends:nt,:][:-chop_ends,:] 
   else:
      b=b[:nt,:]
      vi=np.loadtxt(dirname+'/Vi_'+suffix)[:nt,:]
      ve=np.loadtxt(dirname+'/Ve_'+suffix)[:nt,:]
      ni=np.loadtxt(dirname+'/Ni_'+suffix)[:nt,:]
      ne=np.loadtxt(dirname+'/Ne_'+suffix)[:nt,:]
      R =np.loadtxt(dirname+'/R'  +suffix)[:nt,:]

#
   if hampel_filter is not None:
      from .Hampel import spectral_hampel_scalar, spectral_hampel_vector
      print('Running Spectral Hampel filter on ...')
      print('... p+ velocity')
      vi[:,1],vi[:,2],vi[:,3]=spectral_hampel_vector(vi[:,1],vi[:,2],vi[:,3],dt,\
               hampel_filter['p+'][0],hampel_filter['p+'][1]) 
      print('... e- velocity')
      ve[:,1],ve[:,2],ve[:,3]=spectral_hampel_vector(ve[:,1],ve[:,2],ve[:,3],dt,\
               hampel_filter['e-'][0],hampel_filter['e-'][1]) 
      print('... p+ density')
      ni[:,1] = spectral_hampel_scalar(ni[:,1],dt,hampel_filter['p+'][0],hampel_filter['p+'][1])
      print('... e- density')
      ne[:,1] = spectral_hampel_scalar(ne[:,1],dt,hampel_filter['e-'][0],hampel_filter['e-'][1])

      if 'chop_edges' in hampel_filter:
         chop_edges=hampel_filter['chop_edges'] 
         b = b[chop_edges:,:][:-chop_edges,:] 
         vi=vi[chop_edges:,:][:-chop_edges,:] 
         ve=ve[chop_edges:,:][:-chop_edges,:] 
         ni=ni[chop_edges:,:][:-chop_edges,:] 
         ne=ne[chop_edges:,:][:-chop_edges,:] 
         R = R[chop_edges:,:][:-chop_edges,:] 


   if lowpass_freq is not None and type(lowpass_freq) is not dict:
      temp=lowpass_freq
      lowpass_freq={}
      for i in 'vi','ve','ni','ne':
         lowpass_freq[i]=temp
      temp=None
   if lowpass_freq is not None:
      print('Applying low-pass filter on ...')
      print('vi at ',str(lowpass_freq['vi']),' ... and ...')
      print('ve at ',str(lowpass_freq['ve']),' ... ')
      for i in range(1,4): 
         vi[:,i]=tsa.lowpass_sharp(vi[:,i],dt,lowpass_freq['vi'])
         ve[:,i]=tsa.lowpass_sharp(ve[:,i],dt,lowpass_freq['ve'])
      print('ni at ',str(lowpass_freq['ni']),' ... ')
      ni[:,1]=tsa.lowpass_sharp(ni[:,1],dt,lowpass_freq['ni'])
      print('ne at ',str(lowpass_freq['ne']),' ... ')
      ne[:,1]=tsa.lowpass_sharp(ne[:,1],dt,lowpass_freq['ne'])

   df=pd.DataFrame({'t':R[:,0],'x':R[:,1],'y':R[:,2],'z':R[:,3],'bx':b[:,1],\
                    'by':b[:,2],'bz':b[:,3],'vix':vi[:,1],'viy':vi[:,2],\
                    'viz':vi[:,3],'ni':ni[:,1],'vex':ve[:,1],'vey':ve[:,2],\
                    'vez':ve[:,3],'ne':ne[:,1]})
   return df

