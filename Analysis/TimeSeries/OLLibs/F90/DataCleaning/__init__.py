from .dc import *

def autobad_df(df,columns,*args):
   if type(columns) not in [list, tuple]: 
      columns=columns.split()
   for i in columns:
      df[i]=autobad(df[i],*args)
   return df
def hampel_df(df,columns,*args):
   if type(columns) not in [list, tuple]: 
      columns=columns.split()
   for i in columns:
      df[i]=hampel(df[i],*args)
   return df
def fill_gaps_df(df,*args,columns):
   if type(columns) not in [list, tuple]: 
      columns=columns.split()
   for i in columns:
      df[i]=fill_gaps(df['time'],*args,columns)
   return df
   
