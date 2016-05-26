/*This program evaluates the differential equations for a photon's geodesic in a perturbed Minkowski spacetime in SPHERICAL coordinates restricted to RADIAL MOTION only.*/

/*This program solves the particular case for the Minkowski perturbed spacetime with metric: $g_{ab} = {\eta}_{ab} + h_{ab}$. Where $h_{ab}$ corresponds to the perturbation in the Newtonian weak field limit. A Plummer potential with adequate parameters have been used to simulate the perturbation.
The equations are written in the form $\frac{d(x or p)^{\alpha}}{d\lambda}=f(x^{\alpha},p^{\alpha})$.
Where $p^{\alpha}={\dot{x}}^{\alpha}$ and the indice $\alpha$ runs from 0 to 3.
This system can be restricted by initial conditions in $p^{\theta} = p^{\phi} = 0$ to radial motion only so $\theta$ and $\phi$ are fixed in evolution.
The coordinates for the photon's geodesics are then: (ct,r) = (x0,x1).*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define A 10.0     //Distance parameter of the perturbations
#define G 43007.01     //Gravitational constant
#define M 15000.0     //Mass of the perturbation
#define C 299792.458  //Speed of light
#define NSTEPS 40000000 //Number of steps for solving geodesics
#define NLINES 100000 //Number of lines in geodesic_solution.dat file
#define DLAMBDA 0.01   //Geodesics parameter step

typedef long double mydbl;

/*Function for the gravitational potential to be used. Potential for Plummer model.*/
mydbl potential(mydbl r)
{
  return -G*M/(sqrtl(A*A + r*r));
}

/*Derivative of potential respect to radial coordinate.*/
mydbl der_potential(mydbl r)
{
  return G*M*r/(powl(A*A+r*r, 1.5));
}

/*Function of the 0th momentum component differential equation for the geodesics.
${p0}^{dot} = f0(x^{\alpha},p^{\alpha})$.*/
mydbl geodesic_equation_0(mydbl p0, mydbl pr, mydbl r)
{
  mydbl f = -(2.0/(C*C))*der_potential(r)*p0*pr/(1.0 + 2.0*potential(r)/(C*C));
  return f;
}

/*Function of the 1th (radial) momentum component differential equation for the geodesics.
${p1}^{dot} = f1(x^{\alpha},p^{\alpha})$.*/
mydbl geodesic_equation_r(mydbl p0, mydbl pr, mydbl r)
{
  mydbl f =  (der_potential(r)/(C*C))*(pr*pr - p0*p0)/(1.0 - 2.0*potential(r)/(C*C));
  return f;
}

/*Function for solving the geodesics differential equations using 4th order Runge-Kutta method.
Arguments are pointers so variables in that memory addresses are changed every time this function is called.*/
void runge_kutta_4(mydbl *x0, mydbl *x1, mydbl *p0, mydbl *p1, mydbl *lambda)
{
  /*Increments in the variables of the differential equation we want to solve*/
  mydbl dx0, dx1, dp0, dp1;

  /*dxi = (k1,i + 2*k2,j + 2*k3,j + k4,j)/6. In this sections the ki,j are declared with i=1,2,3,4.*/
  mydbl k1x0, k1x1, k1p0, k1p1;
  mydbl k2x0, k2x1, k2p0, k2p1;
  mydbl k3x0, k3x1, k3p0, k3p1;
  mydbl k4x0, k4x1, k4p0, k4p1;

  /*This section calculates the k1 quantities*/
  k1x0 = *p0*DLAMBDA; k1x1 = *p1*DLAMBDA;
  k1p0 = DLAMBDA*geodesic_equation_0(*p0, *p1, *x1);
  k1p1 = DLAMBDA*geodesic_equation_r(*p0, *p1, *x1);

  /*This section calculates the k2 quantities*/
  k2x0 = DLAMBDA*(*p0 + 0.5*k1p0); k2x1 = DLAMBDA*(*p1 + 0.5*k1p1);
  k2p0 = DLAMBDA*geodesic_equation_0(*p0 + 0.5*k1p0, *p1 + 0.5*k1p1, *x1 + 0.5*k1x1);
  k2p1 = DLAMBDA*geodesic_equation_r(*p0 + 0.5*k1p0, *p1 + 0.5*k1p1, *x1 + 0.5*k1x1);

  /*This section calculates the k3 quantities*/
  k3x0 = DLAMBDA*(*p0 + 0.5*k2p0); k3x1 = DLAMBDA*(*p1 + 0.5*k2p1);
  k3p0 = DLAMBDA*geodesic_equation_0(*p0 + 0.5*k2p0, *p1 + 0.5*k2p1, *x1 + 0.5*k2x1);
  k3p1 = DLAMBDA*geodesic_equation_r(*p0 + 0.5*k2p0, *p1 + 0.5*k2p1, *x1 + 0.5*k2x1);

  /*This section calculates the k4 quantities*/
  k4x0 = DLAMBDA*(*p0 + k3p0); k4x1 = DLAMBDA*(*p1 + k3p1);
  k4p0 = DLAMBDA*geodesic_equation_0(*p0 + k3p0, *p1 + k3p1, *x1 + k3x1);
  k4p1 = DLAMBDA*geodesic_equation_r(*p0 + k3p0, *p1 + k3p1, *x1 + k3x1);

  /*Calculation of the increments*/
  dx0 = (k1x0 + 2.0*k2x0 + 2.0*k3x0 + k4x0)/6.0;
  dx1 = (k1x1 + 2.0*k2x1 + 2.0*k3x1 + k4x1)/6.0;
  dp0 = (k1p0 + 2.0*k2p0 + 2.0*k3p0 + k4p0)/6.0;
  dp1 = (k1p1 + 2.0*k2p1 + 2.0*k3p1 + k4p1)/6.0;

  /*New values of the variables of the differential equation. Since we are using pointers, when called the routine the value of variable change.*/
  *x0 = *x0 + dx0; *x1 = *x1 + dx1;
  *p0 = *p0 + dp0; *p1 = *p1 + dp1;
  
  /*Increment of parameter of geodesics*/
  *lambda = *lambda + DLAMBDA;

}

/*To set the initial value of pr, it must hold $g_{\mu\nu}p^{\mu}p^{\nu} = 0$.
This factor multiplies p0 to guarantee that p1 fulfill the null geodesic condition.*/
mydbl condition_factor(mydbl r)
{
  return sqrtl((1.0+2.0*potential(r)/(C*C))/(1.0 - 2.0*potential(r)/(C*C)));
}

/*$cp^{0}$ multiplied by this factor allows to obtain the energy for a local inertial observer in this spacetime.*/
mydbl energy_factor(mydbl r)
{
  mydbl g = sqrtl(1.0 + 2.0*potential(r)/(C*C));
  return g;
}

/*Violation of null geodesics condition $g_{\mu\nu}p^{\mu}p^{\nu} = 0$.*/
mydbl violation(mydbl r, mydbl p0, mydbl pr)
{
  mydbl f = -(1.0+2.0*potential(r)/(C*C))*p0*p0 + (1.0-2.0*potential(r)/(C*C))*(pr*pr);
  return f;
}

int main(void)
{
  /*Initial conditions*/
  mydbl x0 = 0.0, r = -200.0, p0 = 1.0e-3, pr, lambda = 0.0, energy1, energy, v, difft;
  pr = condition_factor(r)*p0;
  energy1 = C*energy_factor(r)*p0;
  v = violation(r, p0, pr);
  difft = (energy1 - energy1)/energy1;
  
  /*Pointer to file where solution of differential equation will be saved.*/
  FILE *geodesic;
  geodesic = fopen("geodesic_solution.dat","w");

  /*Write line of initial values in file*/
  fprintf(geodesic, "%16.8Le %16.8Le %16.8Le %16.8Le %16.8Le %16.8Le %16.8Le %16.8Le\n", lambda, x0, r, p0, pr, energy1, v, difft);

  long long int ii;   //For counting steps

  /*Solution of the differential equation*/
  for(ii=0; ii< NSTEPS; ii++)
    {
      runge_kutta_4(&x0, &r, &p0, &pr, &lambda);
      if((ii%(NSTEPS/NLINES)) == 0)
	{
	  energy = C*energy_factor(r)*p0;
	  v = violation(r, p0, pr);
	  difft = (energy - energy1)/energy1;
	  fprintf(geodesic, "%16.8Le %16.8Le %16.8Le %16.8Le %16.8Le %16.8Le %16.8Le %16.8Le\n", lambda, x0, r, p0, pr, energy, v, difft);
	}
    }

  /** Releasing all used space in memory **/
  fclose(geodesic); //Close file storing the results
}
