/*This program evaluates the differential equations for a photon's geodesic in a perturbed spacetime in SPHERICAL coordinates and restricted to RADIAL MOTION with a Runge-Kutta-Fehlberg 45 method.*/

/*This program solves the particular case for flat FRW perturbed spacetime with metric: $g_{ab} = {[g]}_{ab} + h_{ab}$. Where ${[g]}_{ab}$ corresponds to the flat FRW metric and $h_{ab}$ corresponds to the perturbation in the conformal Newtonian gauge. A Plummer potential with adequate parameters have been used to simulate the perturbation.
The equations are written in the form $\frac{d(x or p)^{\alpha}}{d\lambda}=f(x^{\alpha},p^{\beta})$ and the indices $\alpha$ and $\beta$ run from 0 to 1 since the motion is only radial.
Where $p^{\alpha}={\dot{x}}^{\alpha}$.
The coordinates for the photon's geodesics are then: (ct,r) = (x0,x1).*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_errno.h>          //GSL error management module
#include <gsl/gsl_spline.h>         //GSL interpolation module

/* PHYSICAL CONSTANTS AND PARAMETERS */
#define A 10.0     //Distance parameter of the perturbations
#define G 43007.01     //Gravitational constant
#define M 0.0     //Mass of the perturbation
#define C 299792.458  //Speed of light

/* PROGRAM PARAMETERS (ONLY USING NLINES FRW) */
#define NLINES 100000 //Number of lines in geodesic_solution.dat file
#define NSTEPS 75000000 //Number of steps for solving geodesics
#define NLINESFRW 10000 //Number of lines in frw.dat file
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

/*Interpolation of function inside *spline object evaluated at an abscisa 'x'.
Argument *spline is a pointer to a spline object which stores the type of interpolation to be made. x is the independent variable where the function is evaluated. *acc is a pointer to a lookup object for interpolations.*/
double interpolator(gsl_spline *spline, double x, gsl_interp_accel *acc)
{
  double a = gsl_spline_eval(spline, x, acc);  //Interpolates data to abcisa x using method in spline and acceleration object acc

  return a;  //Return value of interpolated function at x
}

/*Function for the gravitational potential to be used. Potential for Plummer model.*/
mydbl potential(mydbl r)
{
  return -G*M/(sqrtl(A*A + r*r));
}

/*Derivative of potential respecto to radial coordinate.*/
mydbl der_potential(mydbl r)
{
  return G*M*r/(powl(A*A+r*r, 1.5));
}

/*Function of the 0th momentum component differential equation for the geodesics.
${p0}^{dot} = f0(x^{\alpha},p^{\alpha})$.*/
mydbl geodesic_equation_0(gsl_spline *spline1, gsl_interp_accel *acc1, gsl_spline *spline2, gsl_interp_accel *acc2, mydbl p0, mydbl pr, mydbl x0, mydbl r)
{
  double t = (double)(1.0*x0/C);
  mydbl a = (mydbl) 1.0*interpolator(spline1, t, acc1);
  mydbl adot = (mydbl) 1.0*interpolator(spline2, t, acc2);
  mydbl f = -2.0*der_potential(r)*p0*pr/(C*C + 2.0*potential(r)) - (1.0 - 2.0*potential(r)/(C*C))*(a*adot)*(pr*pr)/(C + 2.0*potential(r)/C);
  return f;
}

/*Function of the 1th (radial) momentum component differential equation for the geodesics.
${p1}^{dot} = f1(x^{\alpha},p^{\alpha})$.*/
mydbl geodesic_equation_r(gsl_spline *spline1, gsl_interp_accel *acc1, gsl_spline *spline2, gsl_interp_accel *acc2, mydbl p0, mydbl pr, mydbl x0, mydbl r)
{
  double t = (double)(1.0*x0/C);
  mydbl a = (mydbl) 1.0*interpolator(spline1, t, acc1);
  mydbl adot = (mydbl) 1.0*interpolator(spline2, t, acc2);
  mydbl f = - (der_potential(r)*p0*p0)/(a*a*(C*C - 2.0*potential(r))) - (2.0*adot*p0*pr)/(C*a) + der_potential(r)*(pr*pr)/(C*C - 2.0*potential(r));
  return f;
}

/*To set the initial value of pr, it must hold $g_{\mu\nu}p^{\mu}p^{\nu} = 0$.
This factor multiplies p0 to guarantee that p1 fulfill the null geodesic condition.*/
mydbl condition_factor(mydbl r, double a)
{
  return (mydbl)(1.0/a)*sqrtl((1.0+2.0*potential(r)/(C*C))/(1.0 - 2.0*potential(r)/(C*C)));
}

/*$cp^{0}$ multiplied by this factor allows to obtain the energy for a local inertial observer in this spacetime.*/
mydbl energy_factor(mydbl r)
{
  mydbl g = sqrtl(1.0 + 2.0*potential(r)/(C*C));
  return g;
}

/*Violation of null geodesics condition $g_{\mu\nu}p^{\mu}p^{\nu} = 0$.*/
mydbl violation(mydbl r, mydbl p0, mydbl pr, double a)
{
  mydbl f = -(1.0+2.0*potential(r)/(C*C))*p0*p0 + (mydbl)(1.0*a*a)*(1.0-2.0*potential(r)/(C*C))*pr*pr;
  return f;
}

/*Function for solving the geodesic's differential equations using a Runge-Kutta-Fehlber method with adaptative stepsize.
Arguments are pointer so variables in that memory addresses are changed every time this function is called. There is no need to define auxiliary variables to store the results of computations.*/
void runge_kutta_fehlberg(gsl_spline *spline1, gsl_interp_accel *acc1, gsl_spline *spline2, gsl_interp_accel *acc2, mydbl *x0, mydbl *x1, mydbl *p0, mydbl *p1, mydbl *lambda, mydbl *dlambda, mydbl lambda_endpoint, mydbl tol, mydbl dlambda_max, mydbl dlambda_min, FILE *fp, double aem, mydbl energy1)
{
  /*Define the k quantities for all the variables*/
  mydbl k1x0, k1x1, k1p0, k1p1;
  mydbl k2x0, k2x1, k2p0, k2p1;
  mydbl k3x0, k3x1, k3p0, k3p1;
  mydbl k4x0, k4x1, k4p0, k4p1;
  mydbl k5x0, k5x1, k5p0, k5p1;
  mydbl k6x0, k6x1, k6p0, k6p1;
  
  /*Define a vector that carries the errors for each differential equation and a variable for storing the biggest error*/
  mydbl r[4], rmax;

  /*Define the increments for the values to be solved in the differential equations*/
  mydbl dp0, dp1, dx0, dx1;

  /*Define some relevant quantities to be printed in a file*/
  mydbl energy, difft, t, v, difference;
  double aobs, difftfrw;

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
      k1p0 = *dlambda*geodesic_equation_0(spline1, acc1, spline2, acc2, *p0, *p1, *x0, *x1);
      k1p1 = *dlambda*geodesic_equation_r(spline1, acc1, spline2, acc2, *p0, *p1, *x0, *x1);
      
      /*This section calculates the k2 quantities*/
      k2x0 = *dlambda*(*p0 + a21*k1p0); k2x1 = *dlambda*(*p1 + a21*k1p1);
      k2p0 = *dlambda*geodesic_equation_0(spline1, acc1, spline2, acc2, *p0 + a21*k1p0, *p1 + a21*k1p1, *x0 + a21*k1x0, *x1 + a21*k1x1);
      k2p1 = *dlambda*geodesic_equation_r(spline1, acc1, spline2, acc2, *p0 + a21*k1p0, *p1 + a21*k1p1, *x0 + a21*k1x0, *x1 + a21*k1x1);
      
      /*This section calculates the k3 quantities*/
      k3x0 = *dlambda*(*p0 + a31*k1p0 + a32*k2p0); k3x1 = *dlambda*(*p1 + a31*k1p1 + a32*k2p1);
      k3p0 = *dlambda*geodesic_equation_0(spline1, acc1, spline2, acc2, *p0 + a31*k1p0 + a32*k2p0, *p1 + a31*k1p1 + a32*k2p1, *x0 + a31*k1x0 + a32*k2x0, *x1 + a31*k1x1 + a32*k2x1);
      k3p1 = *dlambda*geodesic_equation_r(spline1, acc1, spline2, acc2, *p0 + a31*k1p0 + a32*k2p0, *p1 + a31*k1p1 + a32*k2p1, *x0 + a31*k1x0 + a32*k2x0, *x1 + a31*k1x1 + a32*k2x1);
      
      /*This section calculates the k4 quantities*/
      k4x0 = *dlambda*(*p0 + a41*k1p0 + a42*k2p0 + a43*k3p0); k4x1 = *dlambda*(*p1 + a41*k1p1 + a42*k2p1 + a43*k3p1);
      k4p0 = *dlambda*geodesic_equation_0(spline1, acc1, spline2, acc2, *p0 + a41*k1p0 + a42*k2p0 + a43*k3p0, *p1 + a41*k1p1 + a42*k2p1 + a43*k3p1, *x0 + a41*k1x0 + a42*k2x0 + a43*k3x0, *x1 + a41*k1x1 + a42*k2x1 + a43*k3x1);
      k4p1 = *dlambda*geodesic_equation_r(spline1, acc1, spline2, acc2, *p0 + a41*k1p0 + a42*k2p0 + a43*k3p0, *p1 + a41*k1p1 + a42*k2p1 + a43*k3p1, *x0 + a41*k1x0 + a42*k2x0 + a43*k3x0, *x1 + a41*k1x1 + a42*k2x1 + a43*k3x1);

      /*This section calculates the k5 quantities*/
      k5x0 = *dlambda*(*p0 + a51*k1p0 + a52*k2p0 + a53*k3p0 + a54*k4p0); k5x1 = *dlambda*(*p1 + a51*k1p1 + a52*k2p1 + a53*k3p1 + a54*k4p1);
      k5p0 = *dlambda*geodesic_equation_0(spline1, acc1, spline2, acc2, *p0 + a51*k1p0 + a52*k2p0 + a53*k3p0 + a54*k4p0, *p1 + a51*k1p1 + a52*k2p1 + a53*k3p1 + a54*k4p1, *x0 + a51*k1x0 + a52*k2x0 + a53*k3x0 + a54*k4x0, *x1 + a51*k1x1 + a52*k2x1 + a53*k3x1 + a54*k4x1);
      k5p1 = *dlambda*geodesic_equation_r(spline1, acc1, spline2, acc2, *p0 + a51*k1p0 + a52*k2p0 + a53*k3p0 + a54*k4p0, *p1 + a51*k1p1 + a52*k2p1 + a53*k3p1 + a54*k4p1, *x0 + a51*k1x0 + a52*k2x0 + a53*k3x0 + a54*k4x0, *x1 + a51*k1x1 + a52*k2x1 + a53*k3x1 + a54*k4x1);

      /*This section calculates the k6 quantities*/
      k6x0 = *dlambda*(*p0 + a61*k1p0 + a62*k2p0 + a63*k3p0 + a64*k4p0 + a65*k5p0); k6x1 = *dlambda*(*p1 + a61*k1p1 + a62*k2p1 + a63*k3p1 + a64*k4p1 + a65*k5p1);
      k6p0 = *dlambda*geodesic_equation_0(spline1, acc1, spline2, acc2, *p0 + a61*k1p0 + a62*k2p0 + a63*k3p0 + a64*k4p0 + a65*k5p0, *p1 + a61*k1p1 + a62*k2p1 + a63*k3p1 + a64*k4p1 + a65*k5p1, *x0 + a61*k1x0 + a62*k2x0 + a63*k3x0 + a64*k4x0 + a65*k5x0, *x1 + a61*k1x1 + a62*k2x1 + a63*k3x1 + a64*k4x1 + a65*k5x1);
      k6p1 = *dlambda*geodesic_equation_r(spline1, acc1, spline2, acc2, *p0 + a61*k1p0 + a62*k2p0 + a63*k3p0 + a64*k4p0 + a65*k5p0, *p1 + a61*k1p1 + a62*k2p1 + a63*k3p1 + a64*k4p1 + a65*k5p1, *x0 + a61*k1x0 + a62*k2x0 + a63*k3x0 + a64*k4x0 + a65*k5x0, *x1 + a61*k1x1 + a62*k2x1 + a63*k3x1 + a64*k4x1 + a65*k5x1);

      /*Determination of the biggest error between the errors for each differential equation*/
      r[0] = fabsl(r1*k1p0 + r3*k3p0 + r4*k4p0 + r5*k5p0 + r6*k6p0)/(*dlambda);
      r[1] = fabsl(r1*k1p1 + r3*k3p1 + r4*k4p1 + r5*k5p1 + r6*k6p1)/(*dlambda);
      r[2] = fabsl(r1*k1x0 + r3*k3x0 + r4*k4x0 + r5*k5x0 + r6*k6x0)/(*dlambda);
      r[3] = fabsl(r1*k1x1 + r3*k3x1 + r4*k4x1 + r5*k5x1 + r6*k6x1)/(*dlambda);
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
  	  t = *x0/C;
  	  aobs = interpolator(spline1, (double)(1.0*t), acc1);
  	  v = violation(*x1, *p0, *p1, aobs);
  	  difftfrw = (aem/aobs) - 1.0;
  	  difference = difft - (mydbl)(1.0*difftfrw);
  	  fprintf(fp, "%16.8Le %16.8Le %16.8Le %16.8Le %16.8Le %16.8Le %16.8Le %16.8Le %16.8Le %16.8e %16.8Le %16.8Le\n", *lambda, *dlambda, *x0, *x1, *p0, *p1, energy, v, difft, difftfrw, difference, rmax);
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
  /***READ SCALE FACTOR DATA AND PREPARE OBJECTS FOR INTERPOLATION ***/
  int i;            //For array manipulation
  
  /*Pointer to scale_factor.data file*/
  FILE *frw;        
  frw = fopen("../scale_factor.dat","r");

  /*Variables and arrays to read the data*/
  double cosmictime[NLINESFRW], conftime, scale_factor[NLINESFRW], der_scale_factor[NLINESFRW];

  /*Reading the data*/
  for(i=0; i<NLINESFRW; i++)
    {
      fscanf(frw,"%lf %lf %lf %lf", &cosmictime[i], &conftime, &scale_factor[i], &der_scale_factor[i]);
    }

  /*Free space in memory*/
  fclose(frw);
  
  /*** Initializes objects for interpolation. 1 is for interpolation of scale factor, 2 is for interpolation of derivative of scale factor ***/

  /*Allocate space in memory*/
  gsl_interp_accel *acc1 = gsl_interp_accel_alloc();  //Acceleration type object (for index lookup)
  gsl_interp_accel *acc2 = gsl_interp_accel_alloc();
  gsl_spline *spline1 = gsl_spline_alloc(gsl_interp_cspline, NLINESFRW);  //Spline type object (define interpolation type and space in memory, works for both)
  gsl_spline *spline2 = gsl_spline_alloc(gsl_interp_cspline, NLINESFRW);

  /*Initializes objects for interpolation*/
  gsl_spline_init(spline1, cosmictime, scale_factor, NLINESFRW);  //Initializes spline object for data cosmictime, scale_factor of size NLINES
  gsl_spline_init(spline2, cosmictime, der_scale_factor, NLINESFRW);  //Initializes spline object for data cosmictime, der_scale_factor of size NLINES

  /************************************************************************************/
  
  /***SOLVES GEODESIC EQUATIONS FOR PERTURBED FRW UNIVERSE WITH STATIC PLUMMER POTENTIAL ***/

  /*Parameters for RKF method*/
  mydbl dlambda_max = 1.0e-01, dlambda_min = 1.0e-05, tol = 1.0e-10, lambdaEndPoint = 1.0e+05;

  /*Initial conditions*/
  mydbl ti = 7.0 ,x0, r = -500.0, p0 = 1.0e-3, pr, lambda = 0.0, energy1, energy, v, difft, difference, dlambda = dlambda_max;
  double difftfrw, aem, aobs;
  x0 = C*ti;
  aem = interpolator(spline1, (double)(1.0*ti), acc1);
  pr = condition_factor(r, aem)*p0;
  energy1 = C*energy_factor(r)*p0;
  v = violation(r, p0, pr, aem);
  difft = (energy1 - energy1)/energy1;
  difftfrw = (aem/aem) - 1.0;
  difference = difft - (mydbl)(1.0*difftfrw);
  
  /*Pointer to file where solution of differential equation will be saved.*/
  FILE *geodesic;
  geodesic = fopen("geodesic_solution.dat","w");

  /*Write line of initial values in file*/
  fprintf(geodesic, "%16.8Le %16.8Le %16.8Le %16.8Le %16.8Le %16.8Le %16.8Le %16.8Le %16.8Le %16.8e %16.8Le %16.8Le\n", lambda, dlambda, x0, r, p0, pr, energy1, v, difft, difftfrw, difference, (mydbl)0.0);

  /*Solution of the differential equation*/
  runge_kutta_fehlberg(spline1, acc1, spline2, acc2, &x0, &r, &p0, &pr, &lambda, &dlambda, lambdaEndPoint, tol, dlambda_max, dlambda_min, geodesic, aem, energy1);

  /************************************************************************************/

  /***RELEASING ALL SPACE USED IN MEMORY***/
  fclose(geodesic); //Close file storing the results
  gsl_spline_free(spline1);  //Free memory of spline object
  gsl_spline_free(spline2);
  gsl_interp_accel_free(acc1);  //Free memory of accel object
  gsl_interp_accel_free(acc2);
}
