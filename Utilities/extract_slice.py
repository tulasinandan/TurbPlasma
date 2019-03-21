#!/usr/bin/env python
def create_slice():
    import numpy as np
    from TurbPlasma.Utilities.subs import create_object
    from scipy.ndimage import gaussian_filter as gf
    # Create the P3D-Old Object
    rc=create_object()
    # Ask for variables to extract
    vars2ext=input("Which variables to extract? e.g. all or bx by bz etc. ").split()
    # Set those variables to be loaded in P3D object
    rc.vars2load(vars2ext)
    # Ask for time slice to read
    slice2ext=int(input("Which slice out of "+str(rc.numslices)+"? "))
    # Ask if want to smooth data
    smooth=int(input("How much smoothing (integer, 0 for none)? "))
    if smooth == '': smooth=0
    # Load the time slice
    rc.loadslice(slice2ext)
    
    # Write the variables to a file
    for i in rc.vars2l:
       filename=rc.dirname+"."+i+"."+str(slice2ext)+".dat"
       print(filename)
       gf(rc.__dict__[i],sigma=smooth).tofile(filename)

if __name__=="__main__":
    create_slice()
