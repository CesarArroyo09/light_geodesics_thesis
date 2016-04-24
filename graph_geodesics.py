#!/usr/bin/env python

#Import packages
import numpy as np
import matplotlib.pyplot as plt

#Open data file
f1 = open('geodesic_solution.dat','r')
f2 = open('geodesic_solution_rk.dat','r')

#Read file
lines1 = f1.readlines()
lines2 = f2.readlines()

#Close file
f1.close()
f2.close()

#Variables to store information
lamb1 = []
lamb2 = []
x01 = []
x02 = []
x2 = []
p2 = []
energy1 = []
v1 = []
energy2 = []
v2 = []

#Scan rows
for line in lines1:
    p = line.split()
    lamb1.append(float(p[0]))
    x01.append(float(p[1]))
    x2.append(float(p[3]))
    p2.append(float(p[7]))
    energy1.append(float(p[9]))
    v1.append(float(p[10]))

for line in lines2:
    p = line.split()
    lamb2.append(float(p[0]))
    x02.append(float(p[1]))
    energy2.append(float(p[9]))
    v2.append(float(p[10]))

lamb1g = np.array(lamb1)
x01g = np.array(x01)
x2g = np.array(x2)
p2g = np.array(p2)
energy1g = np.array(energy1)
v1g = np.array(v1)

lamb2g = np.array(lamb2)
x02g = np.array(x02)
energy2g = np.array(energy2)
v2g = np.array(v2)



#Plot data
plt.figure(1)
plt.subplot(211)
#plt.xlabel('Affine parameter')
plt.ylabel(r'$x^2$')
plt.title('y coordinate as a function of lambda')
plt.plot(lamb1g,x2g,'ro',markersize=6)
plt.subplot(212)
plt.xlabel('Affine parameter')
plt.ylabel(r'$p^2$')
plt.title('Momentum in y direction as function of lambda')
plt.plot(lamb1g,p2g,'ro',markersize=6)
plt.show()

plt.figure(2)
plt.subplot(211)
plt.ylabel('Energy')
plt.title(r'Energy as a function of $x^0=t$ for Euler method')
plt.plot(x01g,energy1g,'ro',markersize=6)
plt.subplot(212)
plt.xlabel(r'$x^0$')
plt.ylabel('Energy')
plt.title(r'Energy as a function of $x^0=ct$ for Runge-Kutta method')
plt.plot(x02g,energy2g,'ro',markersize=6)
plt.show()

plt.figure(3)
plt.subplot(211)
plt.ylabel(r'g_{\mu\nu} p^{\mu} p^{\nu}')
plt.title('Violation of the null geodesic condition for Euler method')
plt.plot(lamb1g,v1g,'ro',markersize=6)
plt.subplot(212)
plt.xlabel('Affine parameter')
plt.ylabel(r'g_{\mu\nu} p^{\mu} p^{\nu}')
plt.title('Violation of the null geodesic condition for Runge-Kutta method')
plt.plot(lamb2g,v2g,'ro',markersize=6)
plt.show()
