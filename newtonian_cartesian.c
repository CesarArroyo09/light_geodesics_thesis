/*This program evaluates the differential equations for a photon's geodesic in a perturbed spacetime in CARTESIAN coordinates.
This program solves the particular case for the Minkowski perturbed spacetime with metric: $g_{ab} = {\eta}_{ab} + h_{ab}$. Where $h_{ab}$ corresponds to the perturbation in the Newtonian weak field limit. A Plummer potential with adequate parameters have been used to simulate the perturbation.
The equations are written in the form $\frac{d(x or p)^{\alpha}}{d\lambda}=f(x^{\alpha},p^{\alpha})$.
Where $p^{\alpha}={\dot{x}}^{\alpha}$ and the indice $\alpha$ runs from 0 to 3.
The coordinates for the photon's geodesics are: (ct,x,y,z) = (x0,x1,x2,x3).*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define A 1.0     //Distance parameter of the perturbations
#define G 43007.01     //Gravitational constant
#define M 15000.0     //Mass of the perturbation
#define C 299792.458  //Speed of light
#define NLINES 4000000 //Number of lines in frw.dat file
#define DLAMBDA 0.001   //Geodesics parameter step

typedef long double mydbl;

/*Function for the gravitational potential used. Potential for Plummer model.*/
mydbl potential(mydbl x1, mydbl x2, mydbl x3)
{
  return -G*M/(sqrtl(A*A + x1*x1 + x2*x2 +x3*x3));
}

/*First partial derivative of potential respect to the ith coordinate.
First argument xi denotes the ith coordinate.
The xj's denote the other coordinates. The order doesn't matter since they appear in the same way in the equation.*/
mydbl ith_der_potential(mydbl xi, mydbl xj1, mydbl xj2)
{
  return G*M*xi/(powl(A*A+xi*xi+xj1*xj1+xj2*xj2, 1.5));
}

/*Function of the 0th momentum component differential equation for the geodesics.
${p0}^{dot} = f0(x^{\alpha},p^{\alpha})$.*/
mydbl geodesic_equation_0(mydbl p0, mydbl p1, mydbl p2, mydbl p3, mydbl x0, mydbl x1, mydbl x2, mydbl x3)
{
  mydbl f = -2.0*p0*(1/(C*C))*(p1*ith_der_potential(x1, x2, x3) + p2*ith_der_potential(x2, x1, x3) + p3*ith_der_potential(x3, x2, x1));
  return f;
}

/*Function for the ith momentum component differential equation for the geodesics, where i=1,2,3.
xi is the ith coordinate, pi is the ith momentum.
The other position and momentum quantities are denoted xj1, xj2, pj1 and pj2. Quantities xj1 and xj2 appear in symmetric way in this function so order doesn't matter, the same way for pj1 and pj2. Nevertheless, when calling this function the position of xj1 respecto to the 'x' quantities in arguments should be the same that position of pj1 respect to the 'p' quantities in the argument.*/
mydbl geodesic_equation_i(mydbl p0, mydbl pi, mydbl pj1, mydbl pj2, mydbl x0, mydbl xi, mydbl xj1, mydbl xj2)
{
  mydbl f =  (1/powl(C,2))*( 2*pi * (pi*ith_der_potential(xi,xj1,xj2) + pj1*ith_der_potential(xj1,xi,xj2) + pj2*ith_der_potential(xj2,xj1,xi)) - ith_der_potential(xi,xj1,xj2) * (p0*p0 + pi*pi + pj1*pj1 + pj2*pj2) ) ;
  return f;
}

void runge_kutta_4(mydbl *x0, mydbl *x1, mydbl *x2, mydbl *x3, mydbl *p0, mydbl *p1, mydbl *p2, mydbl *p3, mydbl *lambda)
{
  /*Increment in the variables of the differential equation we want to solve*/
  mydbl dx0, dx1, dx2, dx3, dp0, dp1, dp2, dp3;

  /*dxi = (k1,i + 2*k2,j + 2*k3,j + k4,j)/6. In this sections the ki,j are declared with i=1,2,3,4.*/
  mydbl k1x0, k1x1, k1x2, k1x3, k1p0, k1p1, k1p2, k1p3;
  mydbl k2x0, k2x1, k2x2, k2x3, k2p0, k2p1, k2p2, k2p3;
  mydbl k3x0, k3x1, k3x2, k3x3, k3p0, k3p1, k3p2, k3p3;
  mydbl k4x0, k4x1, k4x2, k4x3, k4p0, k4p1, k4p2, k4p3;

  /*This section calculates the k1 quantities*/
  k1x0 = *p0*DLAMBDA; k1x1 = *p1*DLAMBDA; k1x2 = *p2*DLAMBDA; k3x0 = *p3*DLAMBDA;
  k1p0 = geodesic_equation_0(*p0, *p1, *p2, *p3, *x0, *x1, *x2, *x3)*DLAMBDA;
  k1p1 = geodesic_equation_i(*p0, *p1, *p2, *p3, *x0, *x1, *x2, *x3)*DLAMBDA;
  k2p2 = geodesic_equation_i(*p0, *p2, *p1, *p3, *x0, *x2, *x1, *x3)*DLAMBDA;
  k3p3 = geodesic_equation_i(*p0, *p3, *p2, *p1, *x0, *x3, *x2, *x1)*DLAMBDA;

  /*This section calculates the k2 quantities*/
  k2x0 = DLAMBDA*(*p0 + 0.5*k1x0); k2x1 = DLAMBDA*(*p1 + 0.5*k1x1); k2x2 = DLAMBDA*(*p2 + 0.5*k1x2); k2x3 = DLAMBDA*(*p3 + 0.5*k1x3);
  k2p0 = DLAMBDA*geodesic_equation_0(*p0 + 0.5*k1p0, *p1 + 0.5*k1p1, *p2 + 0.5*k1p2, *p3 + 0.5*k1p3, *x0 + 0.5*k1x0, *x1 + 0.5*k1x1, *x2 + 0.5*k1x2, *x3 + 0.5*k1x3);
  k2p1 = DLAMBDA*geodesic_equation_i(*p0 + 0.5*k1p0, *p1 + 0.5*k1p1, *p2 + 0.5*k1p2, *p3 + 0.5*k1p3, *x0 + 0.5*k1x0, *x1 + 0.5*k1x1, *x2 + 0.5*k1x2, *x3 + 0.5*k1x3);
  k2p2 = DLAMBDA*geodesic_equation_i(*p0 + 0.5*k1p0, *p2 + 0.5*k1p2, *p1 + 0.5*k1p1, *p3 + 0.5*k1p3, *x0 + 0.5*k1x0, *x2 + 0.5*k1x2, *x1 + 0.5*k1x1, *x3 + 0.5*k1x3);
  k2p3 = DLAMBDA*geodesic_equation_i(*p0 + 0.5*k1p0, *p3 + 0.5*k1p3, *p2 + 0.5*k1p2, *p1 + 0.5*k1p1, *x0 + 0.5*k1x0, *x3 + 0.5*k1x3, *x2 + 0.5*k1x2, *x1 + 0.5*k1x1);
  
  /*This section calculates the k3 quantities*/
  k3x0 = DLAMBDA*(*p0 + 0.5*k2x0); k3x1 = DLAMBDA*(*p1 + 0.5*k2x1); k3x2 = DLAMBDA*(*p2 + 0.5*k2x2); k3x3 = DLAMBDA*(*p3 + 0.5*k2x3);
  k3p0 = DLAMBDA*geodesic_equation_0(*p0 + 0.5*k2p0, *p1 + 0.5*k2p1, *p2 + 0.5*k2p2, *p3 + 0.5*k2p3, *x0 + 0.5*k2x0, *x1 + 0.5*k2x1, *x2 + 0.5*k2x2, *x3 + 0.5*k2x3);
  k3p1 = DLAMBDA*geodesic_equation_i(*p0 + 0.5*k2p0, *p1 + 0.5*k2p1, *p2 + 0.5*k2p2, *p3 + 0.5*k2p3, *x0 + 0.5*k2x0, *x1 + 0.5*k2x1, *x2 + 0.5*k2x2, *x3 + 0.5*k2x3);
  k3p2 = DLAMBDA*geodesic_equation_i(*p0 + 0.5*k2p0, *p2 + 0.5*k2p2, *p1 + 0.5*k2p1, *p3 + 0.5*k2p3, *x0 + 0.5*k2x0, *x2 + 0.5*k2x2, *x1 + 0.5*k2x1, *x3 + 0.5*k2x3);
  k3p3 = DLAMBDA*geodesic_equation_i(*p0 + 0.5*k2p0, *p3 + 0.5*k2p3, *p2 + 0.5*k2p2, *p1 + 0.5*k2p1, *x0 + 0.5*k2x0, *x3 + 0.5*k2x3, *x2 + 0.5*k2x2, *x1 + 0.5*k2x1);

  /*This section calculates the k4 quantities*/
  k4x0 = DLAMBDA*(*p0 + k3x0); k4x1 = DLAMBDA*(*p1 + k3x1); k4x2 = DLAMBDA*(*p2 + k3x2); k4x3 = DLAMBDA*(*p3 + k3x3);
  k4p0 = DLAMBDA*geodesic_equation_0(*p0 + k3p0, *p1 + k3p1, *p2 + k3p2, *p3 + k3p3, *x0 + k3x0, *x1 + k3x1, *x2 + k3x2, *x3 + k3x3);
  k4p1 = DLAMBDA*geodesic_equation_i(*p0 + k3p0, *p1 + k3p1, *p2 + k3p2, *p3 + k3p3, *x0 + k3x0, *x1 + k3x1, *x2 + k3x2, *x3 + k3x3);
  k4p2 = DLAMBDA*geodesic_equation_i(*p0 + k3p0, *p2 + k3p2, *p1 + k3p1, *p3 + k3p3, *x0 + k3x0, *x2 + k3x2, *x1 + k3x1, *x3 + k3x3);
  k4p3 = DLAMBDA*geodesic_equation_i(*p0 + k3p0, *p3 + k3p3, *p2 + k3p2, *p1 + k3p1, *x0 + k3x0, *x3 + k3x3, *x2 + k3x2, *x1 + k3x1);

  /*Calculation of the increments*/
  dx0 = (k1x0 + 2.0*k2x0 + 2.0*k3x0 + k4x0)/6.0;
  dx1 = (k1x1 + 2.0*k2x1 + 2.0*k3x1 + k4x1)/6.0;
  dx2 = (k1x2 + 2.0*k2x2 + 2.0*k3x2 + k4x2)/6.0;
  dx3 = (k1x3 + 2.0*k2x3 + 2.0*k3x3 + k4x3)/6.0;
  dp0 = (k1p0 + 2.0*k2p0 + 2.0*k3p0 + k4p0)/6.0;
  dp1 = (k1p1 + 2.0*k2p1 + 2.0*k3p1 + k4p1)/6.0;
  dp2 = (k1p2 + 2.0*k2p2 + 2.0*k3p2 + k4p2)/6.0;
  dp3 = (k1p3 + 2.0*k2p3 + 2.0*k3p3 + k4p3)/6.0;

  /*New values of the variables of the differential equation. Since we are using pointers, when called the routine the value of variable change.*/
  *x0 = *x0 + dx0; *x1 = *x1 + dx1; *x2 = *x2 + dx2; *x3 = *x3 + dx3;
  *p0 = *p0 + dp0; *p1 = *p1 + dp1; *p2 = *p2 + dp2; *p3 = *p3 + dp3;
  
  /*Increment of parameter of geodesics*/
  *lambda = *lambda + DLAMBDA;
}

/*To set the initial value of p1, it must hold $g_{\mu\nu}p^{\mu}p^{\nu} = 0$.
This factor multiplies p0 to guarantee that p1 fulfill the null geodesic condition.*/
mydbl condition_factor(mydbl x1, mydbl x2, mydbl x3)
{
  mydbl g = 1.0 + 2.0*potential(x1,x2,x3)/(C*C);
  return g;
}

/*$cp^{0}$ multiplied by this factor allows to obtain the energy for a local inertial observer in this spacetime.*/
mydbl energy_factor(mydbl x1, mydbl x2, mydbl x3)
{
  mydbl g = 1.0 + potential(x1,x2,x3)/(C*C);
  return g;
}

/*Violation of null geodesics condition $g_{\mu\nu}p^{\mu}p^{\nu} = 0$.*/
mydbl violation(mydbl x1, mydbl x2, mydbl x3, mydbl p0, mydbl p1, mydbl p2, mydbl p3)
{
  mydbl v = -(1.0+2.0*potential(x1,x2,x3)/(C*C))*p0*p0 + (1.0-2.0*potential(x1,x2,x3)/(C*C))*(p1*p1 + p2*p2 + p3*p3);
  return v;
}

int main(void)
{
  int i;            //For array manipulation

  /*Initial conditions*/
  mydbl t = 0.0, x1 = -40.0, x2 = 0.0, x3 = 0.0, p0 = 1.0e-1, p1, p2 = 0.0, p3 = 0.0, lambda = 0.0, energy, v, difft, energy1;
  p1 = condition_factor(x1,x2,x3)*p0;
  energy1 = C*energy_factor(x1,x2,x3)*p0;
  v = violation(x1,x2,x3,p0,p1,p2,p3);
  difft = (energy1 - energy1)/energy1;

  /*Pointer to file where solution of differential equation will be saved.*/
  FILE *geodesic;
  geodesic = fopen("geodesic_solution.dat","w");

  /*Write line of initial values in file*/
  fprintf(geodesic,"%.18Lf %.18Lf %.18Lf %.18Lf %.18Lf %.18Lf %.18Lf %.18Lf %.18Lf %.18Lf %.18Lf %.18Lf\n", lambda, t, x1, x2, x3, p0, p1, p2, p3, energy1, v, difft);

  /*Solution of the differential equation*/
  for(i=0; i<(NLINES); i++)
    {
      runge_kutta_4(&t, &x1, &x2, &x3, &p0, &p1, &p2, &p3, &lambda);
      energy = C*energy_factor(x1,x2,x3)*p0;
      v = violation(x1,x2,x3,p0,p1,p2,p3);
      difft = (energy - energy1)/energy1;
      fprintf(geodesic,"%.18Lf %.18Lf %.18Lf %.18Lf %.18Lf %.18Lf %.18Lf %.18Lf %.18Lf %.18Lf %.18Lf %.18Lf\n", lambda, t, x1, x2, x3, p0, p1, p2, p3, energy, v, difft);
    }

  /** Releasing all used space in memory **/
  fclose(geodesic); //Close file storing the results
}
