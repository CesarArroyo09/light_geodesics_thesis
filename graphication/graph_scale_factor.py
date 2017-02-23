#!/usr/bin/env python

#Import packages
import numpy as np
import matplotlib.pyplot as plt

#Open data file
f1 = open('/home/cesar/light_geodesics_thesis/scale_factor.dat','r')

#Read file
lines1 = f1.readlines()

#Close file
f1.close()

#Variables to store information
cosmictime = []
scalefactor = []

#Scan rows
for line in lines1:
    p = line.split()
    cosmictime.append(float(p[0]))
    scalefactor.append(float(p[2]))


cosmictimeg = np.array(cosmictime)
scalefactorg = np.array(scalefactor)

#Plot data
plt.figure(1)
plt.xlabel(r'Cosmic time (1 time unit $\sim 1.4 \, Gyr$)', fontsize = 14)
plt.ylabel(r'Scale factor $a(t)$', fontsize = 14)
plt.title(r'Scale factor vs Cosmic time ($\Omega_{\Lambda} = 0.7$ and $\Omega_m = 0.3$)', fontsize = 16)
plt.plot(cosmictimeg,scalefactorg,'r-',linewidth = 2.5)
#plt.legend(loc='best')
plt.show()
