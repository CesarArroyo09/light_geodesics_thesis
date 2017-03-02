/*This program evaluates the differential equations for a photon's geodesic in a FRW perturbed spacetime in SPHERICAL coordinates and restricted to photon RADIAL MOTION only. This implementation solves a system with three differential equations where the independent variable is the radial comoving coordinate.*/

/*This program solves the particular case for flat FRW perturbed spacetime with metric: $g_{ab} = {[g]}_{ab} + h_{ab}$. Where ${[g]}_{ab}$ corresponds to the flat FRW metric and $h_{ab}$ corresponds to the perturbation in the conformal Newtonian gauge. A Plummer potential or Hernquist potential with adequate parameters can be used to simulate the perturbation.
The indendent variable is the radial comoving coordinate r.
The dependent variables are $x^0=ct$, $k^0$ and $k^r$.
The equations are written in the form $\frac{dy^{\alpha}}{dr}=f(r,y^{\alpha})$ where $y^{\alpha}$ can be equal to $x^0=ct$, $k^0$ and $k^r$ depending on where $\alpha$ equals to 0,1 or 2.
The basic definition for the momentum of the photon is $k^{\alpha}=frac{x^{\alpha}}{dr}$.*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_errno.h>          //GSL error management module
#include <gsl/gsl_spline.h>         //GSL interpolation module

/* PHYSICAL CONSTANTS AND PARAMETERS */
#define A 10.0     //Distance parameter of the perturbations
#define G 43007.01     //Gravitational constant
#define M 15000.0     //Mass of the perturbation
#define C 299792.458  //Speed of light
#define GAMMA 0.8    //Parameter for halo mass evolution

/* PROGRAM PARAMETERS */
#define NLINES 100000 //Number of lines in geodesic_solution.dat file
#define NSTEPS 80000000 //Number of steps for solving geodesics
#define NLINESFRW 10000 //Number of lines in frw.dat file
#define DR 0.0001   //Geodesics parameter step

typedef long double mydbl;

/*Interpolation of function inside *spline object evaluated at an abscisa 'x'.
Argument *spline is a pointer to a spline object which stores the type of interpolation to be made. x is the independent variable where the function is evaluated. *acc is a pointer to a lookup object for interpolations.*/
double interpolator(gsl_spline *spline, double x, gsl_interp_accel *acc)
{
  double a = gsl_spline_eval(spline, x, acc);  //Interpolates data to abcisa x using method in spline and acceleration object acc
  return a;  //Return value of interpolated function at x
}

/*Mass evolution. Function providing the mass at a given z. Mass at present time is obtained at z=0.*/
mydbl mass(mydbl a)
{
  mydbl z = 1.0/a - 1.0;
  return M*expl(-GAMMA*z);
}

/*Plummer model. Provided two variables pointers, this function stores the potential and the derivative of the potential for the Plummer model in these variables at a given radius.*/
void plummer_model(mydbl a, mydbl adot, mydbl r, mydbl *potential, mydbl *der_potential_radial, mydbl *der_potential_temporal, mydbl *dx0dr)
{
  mydbl rsign = copysignl(1.0,r);
  *potential = -G*mass(a)/(sqrtl(A*A + r*r));
  *der_potential_radial = G*mass(a)*r/(powl(A*A+r*r, 1.5));
  *der_potential_temporal = *potential*GAMMA*adot/(a*a);
  *dx0dr = rsign*a*sqrtl((C*C - 2.0*(*potential))/(C*C + 2.0*(*potential)));
}

/*Hernquist model. Provided two variables pointers, this function stores the potential and the derivative of the potential for the Hernquist model in these variables at a given radius.*/
void hernquist_model(mydbl a, mydbl adot, mydbl r, mydbl *potential, mydbl *der_potential_radial, mydbl *der_potential_temporal, mydbl *dx0dr)
{
  double rsign = copysign(1.0,(double)(1.0*r));
  mydbl absr = fabsl(r);
  if(absr > 2000.0)
    {
      *potential = 0.0;
      *der_potential_radial = 0.0;
      *der_potential_temporal = 0.0;
    }
  else if(rsign == -1.0)
    {
      *potential = -G*mass(a)/(A + absr);
      *der_potential_radial = -G*mass(a)/powl(absr+A,2.0);
      *der_potential_temporal = *potential*GAMMA*adot/(a*a);
    }
  else if(rsign == 1.0)
    {
      *potential = -G*mass(a)/(A + absr);
      *der_potential_radial = G*mass(a)/powl(absr+A,2.0);
      *der_potential_temporal = *potential*GAMMA*adot/(a*a);
    }
  *dx0dr = rsign*a*sqrtl((C*C - 2.0*(*potential))/(C*C + 2.0*(*potential)));
}

/*Function of the $x^0$ differential equation for the geodesics.*/
mydbl geodesic_equation_x0(gsl_spline *spline1, gsl_interp_accel *acc1, mydbl x0, mydbl r, void (*model)(mydbl, mydbl, mydbl, mydbl *, mydbl *, mydbl *, mydbl *))
{
  double t = (double)(1.0*x0/C);
  mydbl a = (mydbl) 1.0*interpolator(spline1, t, acc1);
  mydbl adot = 0.0;
  mydbl dx0dr;
  mydbl potential, der_potential_radial, der_potential_temporal;
  (*model)(a, adot, r, &potential, &der_potential_radial, &der_potential_temporal, &dx0dr);
  return dx0dr;
}

/*Function of the $k^0$ differential equation for the geodesics.*/
mydbl geodesic_equation_k0(gsl_spline *spline1, gsl_interp_accel *acc1, gsl_spline *spline2, gsl_interp_accel *acc2, mydbl k0, mydbl kr, mydbl x0, mydbl r, void (*model)(mydbl, mydbl, mydbl, mydbl *, mydbl *, mydbl *, mydbl *))
{
  double t = (double)(1.0*x0/C);
  mydbl a = (mydbl) 1.0*interpolator(spline1, t, acc1);
  mydbl adot = (mydbl) 1.0*interpolator(spline2, t, acc2);
  mydbl potential, der_potential_radial, der_potential_temporal, dx0dr;
  (*model)(a, adot, r, &potential, &der_potential_radial, &der_potential_temporal, &dx0dr);
  mydbl f = -2.0*k0*(der_potential_temporal*dx0dr/C + der_potential_radial)/(C*C + 2.0*potential) + der_potential_temporal*(dx0dr*k0 + a*a*kr)/(C*C*C) - (a*adot*kr/C)*((C*C - 2.0*potential)/(C*C + 2.0*potential));
  return f;
}

/*Function of the $k^r$ differential equation for the geodesics.*/
mydbl geodesic_equation_kr(gsl_spline *spline1, gsl_interp_accel *acc1, gsl_spline *spline2, gsl_interp_accel *acc2, mydbl k0, mydbl kr, mydbl x0, mydbl r, void (*model)(mydbl, mydbl, mydbl, mydbl *, mydbl *, mydbl *, mydbl *))
{
  double t = (double)(1.0*x0/C);
  mydbl a = (mydbl) 1.0*interpolator(spline1, t, acc1);
  mydbl adot = (mydbl) 1.0*interpolator(spline2, t, acc2);
  mydbl potential, der_potential_radial, der_potential_temporal, dx0dr;
  (*model)(a, adot, r, &potential, &der_potential_radial, &der_potential_temporal, &dx0dr);
  mydbl f = (2.0*k0/C)*(der_potential_temporal/(C*C - 2.0*potential) - adot/a) + (der_potential_radial/(C*C - 2.0*potential))*(kr - k0*dx0dr/(a*a));
  return f;
}

/*Function for solving the geodesics differential equations using 4th order Runge-Kutta method.
Arguments are pointer so variables in that memory addresses are changed every time this function is called.*/
void runge_kutta_4(gsl_spline *spline1, gsl_interp_accel *acc1, gsl_spline *spline2, gsl_interp_accel *acc2, mydbl *k0, mydbl *kr, mydbl *x0, mydbl *r, void (*model)(mydbl, mydbl, mydbl, mydbl *, mydbl *, mydbl *, mydbl *))
{
  /*Increment in the variables of the differential equation we want to solve*/
  mydbl dx0, dk0, dkr;

  /*dxi = (k1,j + 2*k2,j + 2*k3,j + k4,j)/6. In this sections the ki,j are declared with i=1,2,3,4.*/
  mydbl k1x0, k1k0, k1kr;
  mydbl k2x0, k2k0, k2kr;
  mydbl k3x0, k3k0, k3kr;
  mydbl k4x0, k4k0, k4kr;

  /*This section calculates the k1 quantities*/
  k1x0 = DR*geodesic_equation_x0(spline1, acc1, *x0, *r, (*model));
  k1k0 = DR*geodesic_equation_k0(spline1, acc1, spline2, acc2, *k0, *kr, *x0, *r, (*model));
  k1kr = DR*geodesic_equation_kr(spline1, acc1, spline2, acc2, *k0, *kr, *x0, *r, (*model));

  /*This section calculates the k2 quantities*/
  k2x0 = DR*geodesic_equation_x0(spline1, acc1, *x0 + 0.5*k1x0, *r, (*model));
  k2k0 = DR*geodesic_equation_k0(spline1, acc1, spline2, acc2, *k0 + 0.5*k1k0, *kr + 0.5*k1kr, *x0 + 0.5*k1x0, *r, (*model));
  k2kr = DR*geodesic_equation_kr(spline1, acc1, spline2, acc2, *k0 + 0.5*k1k0, *kr + 0.5*k1kr, *x0 + 0.5*k1x0, *r, (*model));

  /*This section calculates the k3 quantities*/
  k3x0 = DR*geodesic_equation_x0(spline1, acc1, *x0 + 0.5*k2x0, *r, (*model));
  k3k0 = DR*geodesic_equation_k0(spline1, acc1, spline2, acc2, *k0 + 0.5*k2k0, *kr + 0.5*k2kr, *x0 + 0.5*k2x0, *r, (*model));
  k3kr = DR*geodesic_equation_kr(spline1, acc1, spline2, acc2, *k0 + 0.5*k2k0, *kr + 0.5*k2kr, *x0 + 0.5*k2x0, *r, (*model));

  /*This section calculates the k4 quantities*/
  k4x0 = DR*geodesic_equation_x0(spline1, acc1, *x0 + 0.5*k3x0, *r, (*model));
  k4k0 = DR*geodesic_equation_k0(spline1, acc1, spline2, acc2, *k0 + k3k0, *kr + k3kr, *x0 + k3x0, *r, (*model));
  k4kr = DR*geodesic_equation_kr(spline1, acc1, spline2, acc2, *k0 + k3k0, *kr + k3kr, *x0 + k3x0, *r, (*model));

  /*Calculation of the increments*/
  dx0 = (k1x0 + 2.0*k2x0 + 2.0*k3x0 + k4x0)/6.0;
  dk0 = (k1k0 + 2.0*k2k0 + 2.0*k3k0 + k4k0)/6.0;
  dkr = (k1kr + 2.0*k2kr + 2.0*k3kr + k4kr)/6.0;

  /*New values of the variables of the differential equation. Since we are using pointers, when called the routine the value of variable change.*/
  *x0 = *x0 + dx0; *k0 = *k0 + dk0; *kr = *kr + dkr;
  
  /*Increment of the radial coordinate*/
  *r = *r + DR;
}

/*To set the initial value of kr, it must hold $g_{\mu\nu}k^{\mu}k^{\nu} = 0$.
This factor multiplies k0 to guarantee that kr fulfill the null geodesic condition.*/
mydbl condition_factor(mydbl r, double a, void (*model)(mydbl, mydbl, mydbl, mydbl *, mydbl *, mydbl *, mydbl *))
{
  mydbl potential, der_potential_radial, der_potential_temporal, dx0dr;
  (*model)((mydbl)(1.0*a), 0.0, r, &potential, &der_potential_radial, &der_potential_temporal, &dx0dr);
  return (mydbl)(1.0/a)*sqrtl((C*C + 2.0*potential)/(C*C - 2.0*potential));
}

/*$ck^0$ multiplied by this factor allows to obtain the energy for a local inertial observer in this spacetime.*/
mydbl energy_factor(mydbl r, double a, void (*model)(mydbl, mydbl, mydbl, mydbl *, mydbl *, mydbl *, mydbl *))
{
  mydbl potential, der_potential_radial, der_potential_temporal, dx0dr;
  (*model)((mydbl)(1.0*a), 0.0, r, &potential, &der_potential_radial, &der_potential_temporal, &dx0dr);
  mydbl g = sqrtl(1.0 + 2.0*potential/(C*C));
  return g;
}

/*Violation of null geodesics condition $g_{\mu\nu}k^{\mu}k^{\nu} = 0$.*/
mydbl violation(mydbl r, mydbl k0, mydbl kr, double a, void (*model)(mydbl, mydbl, mydbl, mydbl *, mydbl *, mydbl *, mydbl *))
{
  mydbl potential, der_potential_radial, der_potential_temporal, dx0dr;
  (*model)((mydbl)(1.0*a), 0.0, r, &potential, &der_potential_radial, &der_potential_temporal, &dx0dr);
  mydbl f = -(1.0+2.0*potential/(C*C))*k0*k0 + (mydbl)(1.0*a*a)*(1.0-2.0*potential/(C*C))*kr*kr;
  return f;
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
  
  /***SOLVES GEODESIC EQUATIONS FOR PERTURBED FRW UNIVERSE WITH A GIVEN POTENTIAL MODEL ***/

  /*Choose model to use for solving geodesics*/
  void (*model_in_use)(mydbl, mydbl, mydbl, mydbl *, mydbl *, mydbl *, mydbl *);
  model_in_use = plummer_model;

  /*Initial conditions*/
  mydbl ti = 9.0 ,x0, r = -4000.0, k0 = 1.0e-03, kr, energy1, energy, v, difft, difference, pfrw, isw;
  double difftfrw, aem, aobs;
  x0 = C*ti;
  aem = interpolator(spline1, (double)(1.0*ti), acc1);
  kr = condition_factor(r, aem, model_in_use)*k0;
  energy1 = C*energy_factor(r, aem, model_in_use)*k0;
  v = violation(r, k0, kr, aem, model_in_use);
  difft = (energy1 - energy1)/energy1;
  difftfrw = (aem/aem) - 1.0;
  difference = difft - (mydbl)(1.0*difftfrw);
  pfrw = ((mydbl)(1.0*aem)*1.0e-03/((mydbl)(1.0*aem)));
  isw = k0/pfrw - 1.0;

  /*Pointer to file where solution of differential equation will be saved.*/
  FILE *geodesic;
  geodesic = fopen("geodesic_solution.dat","w");

  /*Write line of initial values in file*/
  fprintf(geodesic, "%16.8Le %16.8Le %16.8Le %16.8Le %16.8Le %16.8Le %16.8Le %16.8e %16.8Le %16.8Le\n", x0, r, k0, kr, energy1, v, difft, difftfrw, difference, isw);

  long long int ii;

  /*Solution of the differential equation*/
  for(ii=0; ii< NSTEPS; ii++)
    {
      runge_kutta_4(spline1, acc1, spline2, acc2, &k0, &kr, &x0, &r, model_in_use);
      if((ii%(NSTEPS/NLINES)) == 0)
	{
	  ti = x0/C;
	  aobs = interpolator(spline1, (double)(1.0*ti), acc1);
	  energy = C*energy_factor(r, aobs, model_in_use)*k0;
	  v = violation(r, k0, kr, aobs, model_in_use);
	  difft = (energy - energy1)/energy1;
	  difftfrw = (aem/aobs) - 1.0;
	  difference = difft - (mydbl)(1.0*difftfrw);
	  pfrw = (mydbl)(1.0*aem)*1.0e-03/((mydbl)(1.0*aobs));
	  isw = k0/pfrw - 1.0;
	  fprintf(geodesic, "%16.8Le %16.8Le %16.8Le %16.8Le %16.8Le %16.8Le %16.8Le %16.8e %16.8Le %16.8Le\n", x0, r, k0, kr, energy, v, difft, difftfrw, difference, isw);
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
