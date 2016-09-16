/*This program evaluates the differential equations for a photon's geodesic in a perturbed Minkowski spacetime in SPHERICAL coordinates restricted to RADIAL MOTION with a Runge-Kutta-Fehlberg 45 method.*/

/*This program solves the particular case for the Minkowski perturbed spacetime with metric: $g_{ab} = {\eta}_{ab} + h_{ab}$. Where $h_{ab}$ corresponds to the perturbation in the Newtonian weak field limit. A Plummer potential with adequate parameters have been used to simulate the perturbation.
The equations are written in the form $\frac{d(x or p)^{\alpha}}{d\lambda}=f(x^{\alpha},p^{\alpha})$.
Where $p^{\alpha}={\dot{x}}^{\alpha}$ and the indice $\alpha$ runs from 0 to 3.
This system can be restricted by initial conditions in $p^{\theta} = p^{\phi} = 0$ to radial motion only so $\theta$ and $\phi$ are fixed in evolution.
The coordinates for the photon's geodesics are then: (ct,r) = (x0,x1).*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/* PHYSICAL CONSTANTS AND PARAMETERS */
#define A 10.0     //Distance parameter of the perturbations
#define G 43007.01     //Gravitational constant
#define M 15000.0     //Mass of the perturbation
#define C 299792.458  //Speed of light

/* PROGRAM PARAMETERS */
#define NSTEPS 100000000 //Number of steps for solving geodesics
#define NLINES 100000 //Number of lines in geodesic_solution.dat file
#define DLAMBDA 0.01   //Geodesics parameter step

//We define mydbl as a generic variable type which can be changed when needed to: float, double, long double
typedef long double mydbl;

/* RUNGE-KUTTA-FEHLBERG METHOD CONSTANTS */
const mydbl c2 = 0.25, c3 = 3.0/8.0, c4 = 12.0/13.0, c5 = 1.0, c6 = 0.5; //Multiplying increment in t
const mydbl a21 = 0.25, a31 = 3.0/32.0, a32 = 9.0/32.0, a41 = 1932.0/2197.0, a42 = -7200.0/2197.0, a43 = 7296.0/2197.0, a51 = 439.0/216.0, a52 = -8.0, a53 = 3680.0/513.0, a54 = -845.0/4104.0, a61 = -8.0/27.0, a62 = 2.0, a63 = -3544.0/2565.0, a64 = 1859.0/4104.0, a65 = -11.0/40.0; //This adds to the approximated solution in step i. First script means K value and second which K is multiplying
const mydbl r1 = 1.0/360.0, r3 = -128.0/4275.0, r4 = -2197.0/75240.0, r5 = 1.0/50.0, r6 = 2.0/55.0; //Coefficients to calculate R
const mydbl b1 = 25.0/216.0, b3 = 1408.0/2565.0, b4 = 2197.0/4104.0, b5 = -1.0/5.0;

/* FUNCTIONS OF THE PROGRAM */

/*Find and return the maximum value for objects in a vector*/
mydbl maxValue(mydbl arr[], size_t size)
{
  size_t i;
  mydbl maxValue = arr[0];
  
  for(i = 1; i < size; i++)
    {
      if(arr[i] > maxValue)
	{
	  maxValue = arr[i];
	}
    }
  return maxValue;
}

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

/*Function for solving the geodesic's differential equations using a Runge-Kutta-Fehlber method with adaptative stepsize.
Arguments are pointer so variables in that memory addresses are changed every time this function is called. There is no need to define auxiliary variables to store the results of computations.*/
void runge_kutta_fehlberg(mydbl *x0, mydbl *x1, mydbl *p0, mydbl *p1, mydbl *lambda, mydbl *dlambda, mydbl lambda_endpoint, mydbl tol, mydbl dlambda_max, mydbl dlambda_min, FILE *fp, mydbl energy1)
{
  /*Define the k quantities for all the variables*/
  mydbl k1x0, k1x1, k1p0, k1p1;
  mydbl k2x0, k2x1, k2p0, k2p1;
  mydbl k3x0, k3x1, k3p0, k3p1;
  mydbl k4x0, k4x1, k4p0, k4p1;
  mydbl k5x0, k5x1, k5p0, k5p1;
  mydbl k6x0, k6x1, k6p0, k6p1;
  
  /*Define the errors for every differential equation, a vector that carries them and a variable for storing the biggest error*/
  mydbl r[4], rp0, rp1, rx0, rx1, rmax;

  /*Define the increments for the values to be solved in the differential equations*/
  mydbl dp0, dp1, dx0, dx1;

  /*Define some relevant quantities to be printed in a file*/
  mydbl energy, difft, v;

  /*Relevant quantity to determine the new stepsize*/
  mydbl delta;

  while(*lambda < lambda_endpoint)
    {
      /*Check if the independent variable plus the actual stepsize is still in the desired interval*/
      if((*lambda + *dlambda) > lambda_endpoint)
	{
	  *dlambda = lambda_endpoint - *lambda;
	}

      /*Calculate the k quantities for all the variables*/
      /*This section calculates the k1 quantities*/
      k1x0 = (*dlambda)*(*p0); k1x1 = (*dlambda)*(*p1);
      k1p0 = *dlambda*geodesic_equation_0(*p0, *p1, *x1);
      k1p1 = *dlambda*geodesic_equation_r(*p0, *p1, *x1);
      
      /*This section calculates the k2 quantities*/
      k2x0 = *dlambda*(*p0 + a21*k1p0); k2x1 = *dlambda*(*p1 + a21*k1p1);
      k2p0 = *dlambda*geodesic_equation_0(*p0 + a21*k1p0, *p1 + a21*k1p1, *x1 + a21*k1x1);
      k2p1 = *dlambda*geodesic_equation_r(*p0 + a21*k1p0, *p1 + a21*k1p1, *x1 + a21*k1x1);
      
      /*This section calculates the k3 quantities*/
      k3x0 = *dlambda*(*p0 + a31*k1p0 + a32*k2p0); k3x1 = *dlambda*(*p1 + a31*k1p1 + a32*k2p1);
      k3p0 = *dlambda*geodesic_equation_0(*p0 + a31*k1p0 + a32*k2p0, *p1 + a31*k1p1 + a32*k2p1, *x1 + a31*k1x1 + a32*k2x1);
      k3p1 = *dlambda*geodesic_equation_r(*p0 + a31*k1p0 + a32*k2p0, *p1 + a31*k1p1 + a32*k2p1, *x1 + a31*k1x1 + a32*k2x1);
      
      /*This section calculates the k4 quantities*/
      k4x0 = *dlambda*(*p0 + a41*k1p0 + a42*k2p0 + a43*k3p0); k4x1 = *dlambda*(*p1 + a41*k1p1 + a42*k2p1 + a43*k3p1);
      k4p0 = *dlambda*geodesic_equation_0(*p0 + a41*k1p0 + a42*k2p0 + a43*k3p0, *p1 + a41*k1p1 + a42*k2p1 + a43*k3p1, *x1 + a41*k1x1 + a42*k2x1 + a43*k3x1);
      k4p1 = *dlambda*geodesic_equation_r(*p0 + a41*k1p0 + a42*k2p0 + a43*k3p0, *p1 + a41*k1p1 + a42*k2p1 + a43*k3p1, *x1 + a41*k1x1 + a42*k2x1 + a43*k3x1);

      /*This section calculates the k5 quantities*/
      k5x0 = *dlambda*(*p0 + a51*k1p0 + a52*k2p0 + a53*k3p0 + a54*k4p0); k5x1 = *dlambda*(*p1 + a51*k1p1 + a52*k2p1 + a53*k3p1 + a54*k4p1);
      k5p0 = *dlambda*geodesic_equation_0(*p0 + a51*k1p0 + a52*k2p0 + a53*k3p0 + a54*k4p0, *p1 + a51*k1p1 + a52*k2p1 + a53*k3p1 + a54*k4p1, *x1 + a51*k1x1 + a52*k2x1 + a53*k3x1 + a54*k4x1);
      k5p1 = *dlambda*geodesic_equation_r(*p0 + a51*k1p0 + a52*k2p0 + a53*k3p0 + a54*k4p0, *p1 + a51*k1p1 + a52*k2p1 + a53*k3p1 + a54*k4p1, *x1 + a51*k1x1 + a52*k2x1 + a53*k3x1 + a54*k4x1);

      /*This section calculates the k6 quantities*/
      k6x0 = *dlambda*(*p0 + a61*k1p0 + a62*k2p0 + a63*k3p0 + a64*k4p0 + a65*k5p0); k6x1 = *dlambda*(*p1 + a61*k1p1 + a62*k2p1 + a63*k3p1 + a64*k4p1 + a65*k5p1);
      k6p0 = *dlambda*geodesic_equation_0(*p0 + a61*k1p0 + a62*k2p0 + a63*k3p0 + a64*k4p0 + a65*k5p0, *p1 + a61*k1p1 + a62*k2p1 + a63*k3p1 + a64*k4p1 + a65*k5p1, *x1 + a61*k1x1 + a62*k2x1 + a63*k3x1 + a64*k4x1 + a65*k5x1);
      k6p1 = *dlambda*geodesic_equation_r(*p0 + a61*k1p0 + a62*k2p0 + a63*k3p0 + a64*k4p0 + a65*k5p0, *p1 + a61*k1p1 + a62*k2p1 + a63*k3p1 + a64*k4p1 + a65*k5p1, *x1 + a61*k1x1 + a62*k2x1 + a63*k3x1 + a64*k4x1 + a65*k5x1);

      /*Determination of the biggest error between the errors for each differential equation*/
      rp0 = fabsl(r1*k1p0 + r3*k3p0 + r4*k4p0 + r5*k5p0 + r6*k6p0)/(*dlambda);
      rp1 = fabsl(r1*k1p1 + r3*k3p1 + r4*k4p1 + r5*k5p1 + r6*k6p1)/(*dlambda);
      rx0 = fabsl(r1*k1x0 + r3*k3x0 + r4*k4x0 + r5*k5x0 + r6*k6x0)/(*dlambda);
      rx1 = fabsl(r1*k1x1 + r3*k3x1 + r4*k4x1 + r5*k5x1 + r6*k6x1)/(*dlambda);
      r[0] = rp0; r[1] = rp1; r[2] = rx0; r[3] = rx1;
      rmax = maxValue(r, 4);

      /*Check whether the worst approximation is less than the introduced tolerance*/
      if(rmax <= tol)
	{
	  *lambda = *lambda + *dlambda;
	  
	  dp0 = b1*k1p0 + b3*k3p0 + b4*k4p0 + b5*k5p0;
	  dp1 = b1*k1p1 + b3*k3p1 + b4*k4p1 + b5*k5p1;
	  dx0 = b1*k1x0 + b3*k3x0 + b4*k4x0 + b5*k5x0;
	  dx1 = b1*k1x1 + b3*k3x1 + b4*k4x1 + b5*k5x1;
	  *p0 = *p0 + dp0; *p1 = *p1 + dp1; *x0 = *x0 + dx0; *x1 = *x1 + dx1;
	  
	  /*Relevant quantities to save in the .dat file. These only are calculated if the approximation is accepted*/
	  energy = C*energy_factor(*x1)*(*p0);
  	  difft = (energy - energy1)/energy1;
  	  v = violation(*x1, *p0, *p1);
  	  fprintf(fp, "%16.8Le %16.8Le %16.8Le %16.8Le %16.8Le %16.8Le %16.8Le %16.8Le %16.8Le %16.8Le\n", *lambda, *dlambda, *x0, *x1, *p0, *p1, energy, v, difft, rmax);
	}

      /*Calculate a new stepsize*/
      delta = 0.84*powl(tol/rmax, 0.25);

      if(delta <= 0.1)
	{
	  *dlambda = *dlambda*0.1;
	}
      else if(delta >= 4.0)
	{
	  *dlambda = *dlambda*4.0;
	}
      else
	{
	  *dlambda = *dlambda*delta;
	}

      /*Check whether the stepsize is too big or too small*/
      if(*dlambda > dlambda_max)
	{
	  *dlambda = dlambda_max;
	}
      else if(*dlambda < dlambda_min)
	{
	  printf("Minimum value for stepsize dlambda exceeded.\n");
	  break;
	}
    }
}

int main(void)
{
  /***SOLVES GEODESIC EQUATIONS FOR PERTURBED FRW UNIVERSE WITH STATIC PLUMMER POTENTIAL ***/
  
  /*Parameters for RKF method*/
  mydbl dlambda_max = 1.0e-03, dlambda_min = 1.0e-05, tol = 1.0e-07, lambdaEndPoint = 1.0e+04;

  /*Initial conditions*/
  mydbl x0 = 0.0, r = -500.0, p0 = 1.0e-2, pr, lambda = 0.0, energy1, energy, v, difft, dlambda = dlambda_max;
  pr = condition_factor(r)*p0;
  energy1 = C*energy_factor(r)*p0;
  v = violation(r, p0, pr);
  difft = (energy1 - energy1)/energy1;
  
  /*Pointer to file where solution of differential equation will be saved.*/
  FILE *geodesic;
  geodesic = fopen("geodesic_solution.dat","w");

  /*Write line of initial values in file*/
  fprintf(geodesic, "%16.8Le %16.8Le %16.8Le %16.8Le %16.8Le %16.8Le %16.8Le %16.8Le %16.8Le %16.8Le\n", lambda, dlambda, x0, r, p0, pr, energy1, v, difft, (mydbl)0.0);

  /*Solution of the differential equation*/
  runge_kutta_fehlberg(&x0, &r, &p0, &pr, &lambda, &dlambda, lambdaEndPoint, tol, dlambda_max, dlambda_min, geodesic, energy1);

  /** Releasing all used space in memory **/
  fclose(geodesic); //Close file storing the results
}
