/*This program evaluates the differential equations for a photon's geodesic in a perturbed spacetime in SPHERICAL coordinates and restricted to photon RADIAL MOTION only.*/

/*This program solves the particular case for flat FRW perturbed spacetime with metric: $g_{ab} = {[g]}_{ab} + h_{ab}$. Where ${[g]}_{ab}$ corresponds to the flat FRW metric and $h_{ab}$ corresponds to the perturbation in the conformal Newtonian gauge. A Plummer potential with adequate parameters have been used to simulate the perturbation.
The equations are written in the form $\frac{d(x or p)^{\alpha}}{d\lambda}=f(x^{\alpha},p^{\alpha})$ and the indice $\alpha$ runs from 0 to 1 since the motion is only radial.
Where $p^{\alpha}={\dot{x}}^{\alpha}$.
The coordinates for the photon's geodesics are then: (ct,r) = (x0,x1).*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_errno.h>          //GSL error management module
#include <gsl/gsl_spline.h>         //GSL interpolation module

#define A 10.0     //Distance parameter of the perturbations
#define G 43007.01     //Gravitational constant
#define M 0.0     //Mass of the perturbation
#define C 299792.458  //Speed of light
#define NLINES 100000 //Number of lines in geodesic_solution.dat file
#define NSTEPS 60000000 //Number of steps for solving geodesics
#define NLINESFRW 10000 //Number of lines in frw.dat file
#define DLAMBDA 0.01   //Geodesics parameter step

typedef long double mydbl;

/*Interpolation of scale factor at time eta.
Argument *spline is a pointer to a spline object which stores the type of interpolation to be made. eta is value of cosmic time to be evaluated. *acc is a pointer to a lookup object for interpolations.*/
double interp_scale_factor(gsl_spline *spline, double eta, gsl_interp_accel *acc)
{
  double a = gsl_spline_eval(spline, eta, acc);  //Interpolates data to abcisa eta using method in spline and acceleration object acc

  return a;  //Return scale factor at a time eta
}

/*Interpolation of derivative of scale factor at time eta.
Argument *spline is a pointer to a spline object which stores the type of interpolation to be made. eta is value of cosmic time to be evaluated. *acc is a pointer to a lookup object for interpolations.*/
double interp_der_scale_factor(gsl_spline *spline, double eta, gsl_interp_accel *acc)
{
  double adot = gsl_spline_eval(spline, eta, acc);  //Interpolates data to abcisa eta using method in spline and acceleration object acc

  return adot;  //Return derivative of scale factor at a time eta
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
  mydbl a = (mydbl) 1.0*interp_scale_factor(spline1, t, acc1);
  mydbl adot = (mydbl) 1.0*interp_der_scale_factor(spline2, t, acc2);
  mydbl f = -(2.0/(C*C))*der_potential(r)*p0*pr/(1.0 + 2.0*potential(r)/(C*C)) - (1.0 - 2.0*potential(r)/(C*C))*(a*adot/C)*(pr*pr)/(1.0 + 2.0*potential(r)/(C*C));
  return f;
}

/*Function of the 1th (radial) momentum component differential equation for the geodesics.
${p1}^{dot} = f1(x^{\alpha},p^{\alpha})$.*/
mydbl geodesic_equation_r(gsl_spline *spline1, gsl_interp_accel *acc1, gsl_spline *spline2, gsl_interp_accel *acc2, mydbl p0, mydbl pr, mydbl x0, mydbl r)
{
  double t = (double)(1.0*x0/C);
  mydbl a = (mydbl) 1.0*interp_scale_factor(spline1, t, acc1);
  mydbl adot = (mydbl) 1.0*interp_der_scale_factor(spline2, t, acc2);
  mydbl f = - (der_potential(r)*p0*p0)/(a*a*C*C*(1.0 - 2.0*potential(r)/(C*C))) - (2.0*adot*p0*pr)/(C*a) + (der_potential(r)/(C*C))*(pr*pr)/(1.0 - 2.0*potential(r)/(C*C));
  return f;
}

/*Function for solving the geodesics differential equations using 4th order Runge-Kutta method.
Arguments are pointer so variables in that memory addresses are changed every time this function is called.*/
void runge_kutta_4(gsl_spline *spline1, gsl_interp_accel *acc1, gsl_spline *spline2, gsl_interp_accel *acc2, mydbl *x0, mydbl *x1, mydbl *p0, mydbl *p1, mydbl *lambda)
{
  /*Increment in the variables of the differential equation we want to solve*/
  mydbl dx0, dx1, dp0, dp1;

  /*dxi = (k1,j + 2*k2,j + 2*k3,j + k4,j)/6. In this sections the ki,j are declared with i=1,2,3,4.*/
  mydbl k1x0, k1x1, k1p0, k1p1;
  mydbl k2x0, k2x1, k2p0, k2p1;
  mydbl k3x0, k3x1, k3p0, k3p1;
  mydbl k4x0, k4x1, k4p0, k4p1;

  /*This section calculates the k1 quantities*/
  k1x0 = *p0*DLAMBDA; k1x1 = *p1*DLAMBDA;
  k1p0 = DLAMBDA*geodesic_equation_0(spline1, acc1, spline2, acc2, *p0, *p1, *x0, *x1);
  k1p1 = DLAMBDA*geodesic_equation_r(spline1, acc1, spline2, acc2, *p0, *p1, *x0, *x1);

  /*This section calculates the k2 quantities*/
  k2x0 = DLAMBDA*(*p0 + 0.5*k1p0); k2x1 = DLAMBDA*(*p1 + 0.5*k1p1);
  k2p0 = DLAMBDA*geodesic_equation_0(spline1, acc1, spline2, acc2, *p0 + 0.5*k1p0, *p1 + 0.5*k1p1, *x0 + 0.5*k1x0, *x1 + 0.5*k1x1);
  k2p1 = DLAMBDA*geodesic_equation_r(spline1, acc1, spline2, acc2, *p0 + 0.5*k1p0, *p1 + 0.5*k1p1, *x0 + 0.5*k1x0, *x1 + 0.5*k1x1);

  /*This section calculates the k3 quantities*/
  k3x0 = DLAMBDA*(*p0 + 0.5*k2p0); k3x1 = DLAMBDA*(*p1 + 0.5*k2p1);
  k3p0 = DLAMBDA*geodesic_equation_0(spline1, acc1, spline2, acc2, *p0 + 0.5*k2p0, *p1 + 0.5*k2p1, *x0 + 0.5*k2x0, *x1 + 0.5*k2x1);
  k3p1 = DLAMBDA*geodesic_equation_r(spline1, acc1, spline2, acc2, *p0 + 0.5*k2p0, *p1 + 0.5*k2p1, *x0 + 0.5*k2x0, *x1 + 0.5*k2x1);

  /*This section calculates the k4 quantities*/
  k4x0 = DLAMBDA*(*p0 + k3p0); k4x1 = DLAMBDA*(*p1 + k3p1);
  k4p0 = DLAMBDA*geodesic_equation_0(spline1, acc1, spline2, acc2, *p0 + 0.5*k3p0, *p1 + 0.5*k3p1, *x0 + 0.5*k3x0, *x1 + 0.5*k3x1);
  k4p1 = DLAMBDA*geodesic_equation_r(spline1, acc1, spline2, acc2, *p0 + 0.5*k3p0, *p1 + 0.5*k3p1, *x0 + 0.5*k3x0, *x1 + 0.5*k3x1);

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

int main(void)
{
  /***READ SCALE FACTOR DATA AND PREPARE OBJECTS FOR INTERPOLATION ***/
  int i;            //For array manipulation
  
  /*Pointer to frw.data file*/
  FILE *frw;        
  frw = fopen("frw.dat","r");

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

  /*Initial conditions*/
  mydbl ti = 7.0 ,x0, r = -200.0, p0 = 1.0e-3, pr, lambda = 0.0, energy1, energy, v, difft, difference;
  double difftfrw, aem, aobs;
  x0 = C*ti;
  aem = interp_scale_factor(spline1, (double)(1.0*ti), acc1);
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
  fprintf(geodesic, "%16.8Le %16.8Le %16.8Le %16.8Le %16.8Le %16.8Le %16.8Le %16.8Le %16.8e %16.8Le\n", lambda, x0, r, p0, pr, energy1, v, difft, difftfrw, difference);

  long long int ii;

  /*Solution of the differential equation*/
  for(ii=0; ii<(1+ NSTEPS); ii++)
    {
      runge_kutta_4(spline1, acc1, spline2, acc2, &x0, &r, &p0, &pr, &lambda);
      if((ii%(NSTEPS/NLINES)) == 0)
	{
	  energy = C*energy_factor(r)*p0;
	  difft = (energy - energy1)/energy1;
	  ti = x0/C;
	  aobs = interp_scale_factor(spline1, (double)(1.0*ti), acc1);
	  v = violation(r, p0, pr, aobs);
	  difftfrw = (aem/aobs) - 1.0;
	  difference = difft - (mydbl)(1.0*difftfrw);
	  fprintf(geodesic, "%16.8Le %16.8Le %16.8Le %16.8Le %16.8Le %16.8Le %16.8Le %16.8Le %16.8e %16.8Le\n", lambda, x0, r, p0, pr, energy, v, difft, difftfrw, difference);
	} 
    }

  /************************************************************************************/

  /*** Releasing all used space in memory ***/
  fclose(geodesic); //Close file storing the results
  gsl_spline_free(spline1);  //Free memory of spline object
  gsl_spline_free(spline2);
  gsl_interp_accel_free(acc1);  //Free memory of accel object
  gsl_interp_accel_free(acc2);
}
