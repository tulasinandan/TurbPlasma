#!/usr/bin/env python
import sys, os
sys.path.insert(0,os.environ['HOME']+'/AJGAR/TurbPlasma')
import random as rnd
import numpy as np
import numpy.fft as nf
import AnalysisFunctions as af
import pandas as pd
#import seaborn as sns
import matplotlib.pyplot as plt
from . import ftsa
import numba

pi = np.pi

def ms_inc(df1,df2):
   df=df1-df2
   df['dr']=np.sqrt(df['x']**2+df['y']**2+df['z']**2)
   return df

