#!/usr/bin/env python

#Import packages
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

#Open data file
#f1 = open('/home/cesar/light_geodesics_thesis/frw/geodesic_solution_gamma_zero.dat','r')
#f2 = open('/home/cesar/light_geodesics_thesis/perturbed_minkowski/geodesic_solution_minkowski_hernquist_poster.dat','r')
#f3 = open('/home/cesar/light_geodesics_thesis/frw/geodesic_solution_poster.dat','r')
f = open('/home/cesar/light_geodesics_thesis/frw/geodesic_solution_poster.dat','r')


#Read file
#lines1 = f1.readlines()
#lines2 = f2.readlines()
#lines3 = f3.readlines()
lines = f.readlines()

#Close file
#f1.close()
#f2.close()
#f3.close()
f.close()

#Variables to store information
radius = []
difftfrw = []
difftpert = []

#radiusfrw = []
#difftfrw = []
#difftfrwfrw = []
#radiusmink = []
#difftmink = []

#radiusisw = []
#isw =[]

#Scan rows
for line in lines:
    p = line.split()
    radius.append(float(p[2]))
    difftfrw.append(float(p[7]))
    difftpert.append(float(p[8]))

#for line in lines1:
#    p = line.split()
#    radiusfrw.append(float(p[2]))
#    difftfrw.append(float(p[7]))
#    difftfrwfrw.append(float(p[8]))

#for line in lines2:
#    p = line.split()
#    radiusmink.append(float(p[2]))
#    difftmink.append(float(p[7]))

#for line in lines3:
#    p = line.split()
#    radiusisw.append(float(p[2]))
#    isw.append(float(p[10]))

#radiusisw2 = []
#isw2 = []
#length = len(radiusisw)

#for i in range(length):
#    if(radiusisw[i] > 4000.0):
#        radiusisw2.append(radiusisw[i])
#        isw2.append(isw[i])

#iswabs = np.abs(isw2)

#radiusfrwg = np.array(radiusfrw)
#difftfrwg = np.array(difftfrw)
#difftfrwfrwg = np.array(difftfrwfrw)
#radiusminkg = np.array(radiusmink)
#difftminkg = np.array(difftmink)

#radiusisw2g = np.array(radiusisw2)
#isw2g = np.array(iswabs)

radiusg = np.array(radius)
difftfrwg = np.array(difftfrw)
difftpertg = np.array(difftpert)


#Plot data
plt.figure(1)
plt.xlabel('Radius (kpc/h)', fontsize = 18)
plt.ylabel(r'$\Delta E/E$', fontsize = 24)
#plt.title('Fractional changes in energy for FRW with static Hernquist perturbation')
plt.plot(radiusg,difftpertg,'k--',label='FRW', linewidth = 2.5)
plt.plot(radiusg,difftfrwg,'r-',label='FRW with evolving perturbation', linewidth = 2.5)
plt.xlim([-2000,8000])
plt.legend(loc='best')
current_axis = plt.gca()
current_axis.add_patch(patches.Rectangle((7800,-0.003),200,0.0001,fill=None))
plt.show()

plt.figure(2)
plt.xlabel('Radius (kpc/h)', fontsize = 18)
plt.ylabel(r'$\Delta E/E$', fontsize = 24)
#plt.title('Fractional changes in energy for FRW with static Hernquist perturbation')
plt.plot(radiusg,difftpertg,'k--',label='FRW', linewidth = 2.5)
plt.plot(radiusg,difftfrwg,'r-',label='FRW with evolving perturbation', linewidth = 2.5)
plt.xlim([7800,8000])
plt.ylim([-0.003,-0.0029])
plt.legend(loc='best')
plt.show()

#plt.figure(1)
#plt.xlabel('Radius (kpc/h)', fontsize = 16)
#plt.ylabel(r'$\Delta E/E$', fontsize = 20)
#plt.title('Fractional changes in energy for FRW with static Hernquist perturbation')
#plt.plot(radiusfrwg,difftfrwg,'r-',label='FRW with static perturbation', linewidth = 2.5)
#plt.plot(radiusfrwg,difftfrwfrwg,'k--',label='FRW', linewidth = 2.5)
#plt.legend(loc='best')
#plt.show()

#plt.figure(2)
#plt.xlabel('Radius (kpc/h)', fontsize = 16)
#plt.ylabel(r'$\Delta E/E$', fontsize = 20)
#plt.title(r'Fractional change in energy for Minkowski with Hernquist perturbation')
#plt.plot(radiusminkg, difftminkg, 'b-', label='Minkowski with static perturbation', linewidth = 2.5)
#plt.legend(loc='best')
#plt.show()

#plt.figure(3)
#plt.xlabel('Radius (kpc/h)')
#plt.title('Temperature fluctuaction for FRW spacetime with Hernquist sphere')
#plt.yscale('log')
#plt.plot(radiusisw2g, isw2g, 'g-')
#plt.show()