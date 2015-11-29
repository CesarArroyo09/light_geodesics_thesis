/*This program evaluates the differential equations for a photon's geodesic in a perturbed spacetime*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define A 1.0
#define G 1.0
#define M 1.0

/*Function for the gravitational potential to be used*/
double potential(double x1, double x2, double x3)
{
  return G*M/(sqrt(A*A + x1*x1 + x2*x2 +x3*x3));
}

/*First partial derivative respect to the ith coordinate.
First argument xi denotes the ith coordinate.
The xj's denote the other coordinates. The order doesn't matter since they appear in the same way in the equation.*/
double ith_der_potential(double xi, double xj1, double xj2)
{
  return -G*M*xi/(pow(A*A+xi*xi+xj1*xj1+xj2*xj2, 1.5));
}

/*This is the function of the first differential equation for the geodesics.
l^{dot} = f1(variables of the problem)*/
double geodesic_equation_1(double aprime, double a, double p0, double p1, double p2, double p3, double eta, double x1, double x2, double x3)
{
  double f = (aprime/a)*pow(p0,2) + 2*p0*(ith_der_potential(x1,x2,x3)*p1 + ith_der_potential(x2,x1,x3)*p2 + ith_der_potential(x3,x1,x2)*p3) + (aprime/a)*(1-4*potential(x1,x2,x3))*(pow(p1,2) + pow(p2,2) + pow(p3,2));
  return f;
}

/*Function for the second differential equation for the geodesics.*/
double geodesic_equation_i(double aprime, double a, double pi, double pj1, double pj2, double pj3, double eta, double xi, double xj1, double xj2)
{
  double f = ith_der_potential(xi,xj1,xj2)*pow(p0,2) + 2*(aprime/a)*po*pi - ith_der_potential(xi,xj1,xj2)*(pow(pi,2)-pow(pj1,2)-pow(pj2,2)) - 2*ith_der_potential(xj1,xi,xj2)*pi*pj1 - 2*ith_der_potential(xj2,xi,xj1)*pi*pj2;
  return f;
}
