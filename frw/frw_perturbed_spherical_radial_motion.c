/*This program evaluates the differential equations for a photon's geodesic in a perturbed spacetime in SPHERICAL coordinates and restricted to photon RADIAL MOTION only.*/

/*This program solves the particular case for flat FRW perturbed spacetime with metric: $g_{ab} = {[g]}_{ab} + h_{ab}$. Where ${[g]}_{ab}$ corresponds to the flat FRW metric and $h_{ab}$ corresponds to the perturbation in the conformal Newtonian gauge. A Plummer potential or a Hernquist potential with adequate parameters can be used to simulate the perturbation.
The equations are written in the form $\frac{d(x or p)^{\alpha}}{d\lambda}=f(x^{\alpha},p^{\alpha})$ and the indice $\alpha$ runs from 0 to 1 since the motion is only radial.
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
#define M 15000.0     //Mass of the perturbation
#define C 299792.458  //Speed of light
#define GAMMA 0.8    //Parameter for halo mass evolution

/* PROGRAM PARAMETERS (ONLY USING NLINES FRW) */
#define NLINES 100000 //Number of lines in geodesic_solution.dat file
#define NSTEPS 150000000 //Number of steps for solving geodesics
#define NLINESFRW 10000 //Number of lines in frw.dat file
#define DLAMBDA 0.01   //Geodesics parameter step

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
void plummer_model(mydbl a, mydbl r, mydbl *potential, mydbl *der_potential)
{
  *potential = -G*mass(a)/(sqrtl(A*A + r*r));
  *der_potential = G*mass(a)*r/(powl(A*A+r*r, 1.5));
}

/*Hernquist model. Provided two variables pointers, this function stores the potential and the derivative of the potential for the Hernquist model in these variables at a given radius.*/
void hernquist_model(mydbl a, mydbl r, mydbl *potential, mydbl *der_potential)
{
  double rsign = copysign(1.0,(double)(1.0*r));
  mydbl absr = fabsl(r);
  *potential = -G*mass(a)/(A + absr);
  if(rsign == -1.0)
    {
      *der_potential = -G*mass(a)/powl(absr+A,2.0);
    }
  else if(rsign == 1.0)
    {
      *der_potential = G*mass(a)/powl(absr+A,2.0);
    }
}

/*Function of the 0th momentum component differential equation for the geodesics.
${p0}^{dot} = f0(x^{\alpha},p^{\alpha})$.*/
mydbl geodesic_equation_0(gsl_spline *spline1, gsl_interp_accel *acc1, gsl_spline *spline2, gsl_interp_accel *acc2, mydbl p0, mydbl pr, mydbl x0, mydbl r, void (*model)(mydbl, mydbl, mydbl *, mydbl *))
{
  double t = (double)(1.0*x0/C);
  mydbl a = (mydbl) 1.0*interpolator(spline1, t, acc1);
  mydbl adot = (mydbl) 1.0*interpolator(spline2, t, acc2);
  mydbl potential, der_potential;
  (*model)(a, r, &potential, &der_potential);
  mydbl f = -2.0*der_potential*p0*pr/(C*C + 2.0*potential) - (1.0 - 2.0*potential/(C*C))*(a*adot)*(pr*pr)/(C + 2.0*potential/C);
  return f;
}

/*Function of the 1th (radial) momentum component differential equation for the geodesics.
${p1}^{dot} = f1(x^{\alpha},p^{\alpha})$.*/
mydbl geodesic_equation_r(gsl_spline *spline1, gsl_interp_accel *acc1, gsl_spline *spline2, gsl_interp_accel *acc2, mydbl p0, mydbl pr, mydbl x0, mydbl r, void (*model)(mydbl, mydbl, mydbl *, mydbl *))
{
  double t = (double)(1.0*x0/C);
  mydbl a = (mydbl) 1.0*interpolator(spline1, t, acc1);
  mydbl adot = (mydbl) 1.0*interpolator(spline2, t, acc2);
  mydbl potential, der_potential;
  (*model)(a, r, &potential, &der_potential);
  mydbl f = - (der_potential*p0*p0)/(a*a*(C*C - 2.0*potential)) - (2.0*adot*p0*pr)/(C*a) + der_potential*(pr*pr)/(C*C - 2.0*potential);
  return f;
}

/*Function for solving the geodesics differential equations using 4th order Runge-Kutta method.
Arguments are pointer so variables in that memory addresses are changed every time this function is called.*/
void runge_kutta_4(gsl_spline *spline1, gsl_interp_accel *acc1, gsl_spline *spline2, gsl_interp_accel *acc2, mydbl *x0, mydbl *x1, mydbl *p0, mydbl *p1, mydbl *lambda, void (*model)(mydbl, mydbl, mydbl *, mydbl *))
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
  k1p0 = DLAMBDA*geodesic_equation_0(spline1, acc1, spline2, acc2, *p0, *p1, *x0, *x1, (*model));
  k1p1 = DLAMBDA*geodesic_equation_r(spline1, acc1, spline2, acc2, *p0, *p1, *x0, *x1, (*model));

  /*This section calculates the k2 quantities*/
  k2x0 = DLAMBDA*(*p0 + 0.5*k1p0); k2x1 = DLAMBDA*(*p1 + 0.5*k1p1);
  k2p0 = DLAMBDA*geodesic_equation_0(spline1, acc1, spline2, acc2, *p0 + 0.5*k1p0, *p1 + 0.5*k1p1, *x0 + 0.5*k1x0, *x1 + 0.5*k1x1, (*model));
  k2p1 = DLAMBDA*geodesic_equation_r(spline1, acc1, spline2, acc2, *p0 + 0.5*k1p0, *p1 + 0.5*k1p1, *x0 + 0.5*k1x0, *x1 + 0.5*k1x1, (*model));

  /*This section calculates the k3 quantities*/
  k3x0 = DLAMBDA*(*p0 + 0.5*k2p0); k3x1 = DLAMBDA*(*p1 + 0.5*k2p1);
  k3p0 = DLAMBDA*geodesic_equation_0(spline1, acc1, spline2, acc2, *p0 + 0.5*k2p0, *p1 + 0.5*k2p1, *x0 + 0.5*k2x0, *x1 + 0.5*k2x1, (*model));
  k3p1 = DLAMBDA*geodesic_equation_r(spline1, acc1, spline2, acc2, *p0 + 0.5*k2p0, *p1 + 0.5*k2p1, *x0 + 0.5*k2x0, *x1 + 0.5*k2x1, (*model));

  /*This section calculates the k4 quantities*/
  k4x0 = DLAMBDA*(*p0 + k3p0); k4x1 = DLAMBDA*(*p1 + k3p1);
  k4p0 = DLAMBDA*geodesic_equation_0(spline1, acc1, spline2, acc2, *p0 + k3p0, *p1 + k3p1, *x0 + k3x0, *x1 + k3x1, (*model));
  k4p1 = DLAMBDA*geodesic_equation_r(spline1, acc1, spline2, acc2, *p0 + k3p0, *p1 + k3p1, *x0 + k3x0, *x1 + k3x1, (*model));

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
mydbl condition_factor(mydbl r, double a, void (*model)(mydbl, mydbl, mydbl *, mydbl *))
{
  mydbl potential, der_potential;
  (*model)((mydbl)(1.0*a), r, &potential, &der_potential);
  return (mydbl)(1.0/a)*sqrtl((C*C + 2.0*potential)/(C*C - 2.0*potential));
}

/*$cp^{0}$ multiplied by this factor allows to obtain the energy for a local inertial observer in this spacetime.*/
mydbl energy_factor(mydbl r, double a, void (*model)(mydbl, mydbl, mydbl *, mydbl *))
{
  mydbl potential, der_potential;
  (*model)((mydbl)(1.0*a), r, &potential, &der_potential);
  mydbl g = sqrtl(1.0 + 2.0*potential/(C*C));
  return g;
}

/*Violation of null geodesics condition $g_{\mu\nu}p^{\mu}p^{\nu} = 0$.*/
mydbl violation(mydbl r, mydbl p0, mydbl pr, double a, void (*model)(mydbl, mydbl, mydbl *, mydbl *))
{
  mydbl potential, der_potential;
  (*model)((mydbl)(1.0*a), r, &potential, &der_potential);
  mydbl f = -(1.0+2.0*potential/(C*C))*p0*p0 + (mydbl)(1.0*a*a)*(1.0-2.0*potential/(C*C))*pr*pr;
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
  void (*model_in_use)(mydbl, mydbl, mydbl *, mydbl *);
  model_in_use = hernquist_model;

  /*Initial conditions*/
  mydbl ti = 7.0 ,x0, r = -1000.0, p0 = 1.0e-3, pr, lambda = 0.0, energy1, energy, v, difft, difference, pfrw, isw;
  double difftfrw, aem, aobs;
  x0 = C*ti;
  aem = interpolator(spline1, (double)(1.0*ti), acc1);
  pr = condition_factor(r, aem, model_in_use)*p0;
  energy1 = C*energy_factor(r, aem, model_in_use)*p0;
  v = violation(r, p0, pr, aem, model_in_use);
  difft = (energy1 - energy1)/energy1;
  difftfrw = (aem/aem) - 1.0;
  difference = difft - (mydbl)(1.0*difftfrw);
  pfrw = ((mydbl)(1.0*aem)*1.0e-03/((mydbl)(1.0*aem)));
  isw = p0/pfrw - 1.0;

  /*Pointer to file where solution of differential equation will be saved.*/
  FILE *geodesic;
  geodesic = fopen("geodesic_solution.dat","w");

  /*Write line of initial values in file*/
  fprintf(geodesic, "%16.8Le %16.8Le %16.8Le %16.8Le %16.8Le %16.8Le %16.8Le %16.8Le %16.8e %16.8Le %16.8Le\n", lambda, x0, r, p0, pr, energy1, v, difft, difftfrw, difference, isw);

  long long int ii;

  /*Solution of the differential equation*/
  for(ii=0; ii< NSTEPS; ii++)
    {
      runge_kutta_4(spline1, acc1, spline2, acc2, &x0, &r, &p0, &pr, &lambda, model_in_use);
      if((ii%(NSTEPS/NLINES)) == 0)
	{
	  ti = x0/C;
	  aobs = interpolator(spline1, (double)(1.0*ti), acc1);
	  energy = C*energy_factor(r, aobs, model_in_use)*p0;
	  v = violation(r, p0, pr, aobs, model_in_use);
	  difft = (energy - energy1)/energy1;
	  difftfrw = (aem/aobs) - 1.0;
	  difference = difft - (mydbl)(1.0*difftfrw);
	  pfrw = (mydbl)(1.0*aem)*1.0e-03/((mydbl)(1.0*aobs));
	  isw = p0/pfrw - 1.0;
	  fprintf(geodesic, "%16.8Le %16.8Le %16.8Le %16.8Le %16.8Le %16.8Le %16.8Le %16.8Le %16.8e %16.8Le %16.8Le\n", lambda, x0, r, p0, pr, energy, v, difft, difftfrw, difference, isw);
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
