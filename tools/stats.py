#!/usr/bin/env python
import sys
import numpy as np

filename = sys.argv[1]

Cycle, Rmax, R2 = np.loadtxt(filename, usecols = (3, 4, 5), unpack = True)

print ("========================================================")
print ("Statistics for file: ", filename)
print ("========================================================")
print ("Max element residual")
print ("--------------------")
print ("MINIMUM, Cycle#: ", np.nanmin(Rmax), ", ", Cycle[np.nanargmin(Rmax)])
print ("MAXIMUM, Cycle#: ", np.nanmax(Rmax), ", " , Cycle[np.nanargmax(Rmax)])
print ("MEAN: ", np.nanmean(Rmax))
print ("STD: ", np.nanstd(Rmax))
print ("")
print ("Frobenius norm squared residual")
print ("--------------------------------")
print ("MINIMUM, Cycle#: ", np.nanmin(R2), ", ", Cycle[np.nanargmin(R2)])
print ("MAXIMUM, Cycle#: ", np.nanmax(R2), ", " , Cycle[np.nanargmax(R2)])
print ("MEAN: ", np.nanmean(R2))
print ("STD: ", np.nanstd(R2))
print ("")
