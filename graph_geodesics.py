#!/usr/bin/env python

#Import packages
import numpy as np
import matplotlib.pyplot as plt

#Open data file
f1 = open('/home/cesar/light_geodesics_thesis/frw/geodesic_solution_gamma_zero.dat','r')
f2 = open('/home/cesar/light_geodesics_thesis/perturbed_minkowski/geodesic_solution_minkowski_hernquist_poster.dat','r')
f3 = open('/home/cesar/light_geodesics_thesis/frw/geodesic_solution_poster.dat','r')

#Read file
lines1 = f1.readlines()
lines2 = f2.readlines()
lines3 = f3.readlines()

#Close file
f1.close()
f2.close()
f3.close()

#Variables to store information
radiusfrw = []
difftfrw = []
difftfrwfrw = []
radiusmink = []
difftmink = []

radiusisw = []
isw =[]

#Scan rows
for line in lines1:
    p = line.split()
    radiusfrw.append(float(p[2]))
    difftfrw.append(float(p[7]))
    difftfrwfrw.append(float(p[8]))

for line in lines2:
    p = line.split()
    radiusmink.append(float(p[2]))
    difftmink.append(float(p[7]))

for line in lines3:
    p = line.split()
    radiusisw.append(float(p[2]))
    isw.append(float(p[10]))

radiusisw2 = []
isw2 = []
length = len(radiusisw)

for i in range(length):
    if(radiusisw[i] > 4000.0):
        radiusisw2.append(radiusisw[i])
        isw2.append(isw[i])

iswabs = np.abs(isw2)

radiusfrwg = np.array(radiusfrw)
difftfrwg = np.array(difftfrw)
difftfrwfrwg = np.array(difftfrwfrw)
radiusminkg = np.array(radiusmink)
difftminkg = np.array(difftmink)

radiusisw2g = np.array(radiusisw2)
isw2g = np.array(iswabs)


#Plot data
plt.figure(1)
plt.xlabel('Radius (kpc/h)')
plt.ylabel(r'$\Delta T/T$')
plt.title('Fractional changes in energy for FRW with Hernquist perturbation')
plt.plot(radiusfrwg,difftfrwg,'r-',label='FRW with static perturbation')
plt.plot(radiusfrwg,difftfrwfrwg,'g--',label='FRW')
plt.legend(loc='best')
plt.show()

plt.figure(2)
plt.xlabel('Radius (kpc/h)')
plt.title(r'Fractional change in energy for Minkowski with Hernquist perturbation')
plt.plot(radiusminkg, difftminkg, 'b-', label='Minkowski with static perturbation')
plt.show()

plt.figure(3)
plt.xlabel('Radius (kpc/h)')
plt.title('Temperature fluctuaction for FRW spacetime with Hernquist sphere')
plt.yscale('log')
plt.plot(radiusisw2g, isw2g, 'g-')
plt.show()
