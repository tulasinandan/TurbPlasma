#!/usr/bin/env python
import random as rnd
import numpy as np
import numpy.fft as nf
import pandas as pd

pi = np.pi

def ms_inc(df1,df2):
   df=df1-df2
   df['dr']=np.sqrt(df['x']**2+df['y']**2+df['z']**2)
   return df

