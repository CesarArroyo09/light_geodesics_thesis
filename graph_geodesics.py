#!/usr/bin/env python

#Import packages
import numpy as np
import matplotlib.pyplot as plt

#Open data file
f1 = open('geodesic_solution.dat','r')

#Read file
lines1 = f1.readlines()

#Close file
f1.close()

#Variables to store information
radius = []
difft = []
difftfrw = []
difference = []

#Scan rows
for line in lines1:
    p = line.split()
    radius.append(float(p[2]))
    difft.append(float(p[11]))
    difftfrw.append(float(p[12]))
    difference.append(float(p[13]))

radiusg = np.array(radius)
difftg = np.array(difft)
difftfrwg = np.array(difftfrw)
differenceg = np.array(difference)



#Plot data
plt.figure(1)
plt.xlabel('Radius')
plt.ylabel(r'$\Delta T/T$')
plt.title('Fractional changes in energy (step size = 0.01)')
plt.plot(radiusg,difftg,'r-',label='Program')
plt.plot(radiusg,difftfrwg,'g--',label='Theoretical')
plt.legend(loc='best')
plt.show()

plt.figure(2)
plt.xlabel('Radius')
plt.title(r'${(\Delta T/T)}_{program} - {(\Delta T/T)}_{theoretical}$')
plt.plot(radiusg,differenceg)
plt.show()
