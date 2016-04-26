/*This program evaluates the differential equations for a photon's geodesic in a perturbed spacetime.
The equations are written in the form $\frac{d(x or p)^{\alpha}}{d\lambda}=f(x^{\alpha},p^{\alpha})$.
Where $p^{\alpha}={\dot{x}}^{\alpha}$*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_errno.h>          //GSL error management module
#include <gsl/gsl_spline.h>         //GSL interpolation module

#define A 10.0     //Distance parameter of the perturbations
#define G 43007.01     //Gravitational constant
#define M 15000.0     //Mass of the perturbation
#define C 299792.458  //Speed of light
#define NLINES 4000000 //Number of lines in geodesic_solution.dat file
#define NLINESFRW 999 //Number of lines in frw.dat file
#define DLAMBDA 1.0   //Geodesics parameter step

typedef long double mydbl;

/*Interpolation of scale factor at time eta.
Argument *spline is a pointer to a spline object which stores the type of interpolation to be made. eta is value of conformal time to be evaluated. *acc is a pointer to a lookup object for interpolations.*/
double interp_scale_factor(gsl_spline *spline, double eta, gsl_interp_accel *acc)
{
  double a = gsl_spline_eval(spline, eta, acc);  //Interpolates data to abcisa eta using method in spline and acceleration object acc

  return a;  //Return scale factor at a time eta
}

/*Interpolation of derivative of scale factor at time eta.
Argument *spline is a pointer to a spline object which stores the type of interpolation to be made. eta is value of conformal time to be evaluated. *acc is a pointer to a lookup object for interpolations.*/
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
mydbl geodesic_equation_0(gsl_spline *spline1, gsl_interp_accel *acc1, gsl_spline *spline2, gsl_interp_accel *acc2, mydbl p0, mydbl pr, mydbl ptheta, mydbl pphi, mydbl x0, mydbl r, mydbl theta)
{
  double t = (double)(x0/C);
  mydbl a = interp_scale_factor(spline1, t, acc1);
  mydbl adot = interp_der_scale_factor(spline2, t, acc2);
  mydbl f = -(2.0/(C*C))*der_potential(r)*p0*pr - (1.0 - 4.0*potential(r)/(C*C))*((a*adot)/C)*(pr*pr + r*r*(powl(sinl(theta),2.0)*pphi*pphi + ptheta*ptheta));
  return f;
}

/*Function of the 1th (radial) momentum component differential equation for the geodesics.
${p1}^{dot} = f1(x^{\alpha},p^{\alpha})$.*/
mydbl geodesic_equation_r(gsl_spline *spline1, gsl_interp_accel *acc1, gsl_spline *spline2, gsl_interp_accel *acc2, mydbl p0, mydbl pr, mydbl ptheta, mydbl pphi, mydbl x0, mydbl r, mydbl theta)
{
  double t = (double)(x0/C);
  mydbl a = (mydbl)interp_scale_factor(spline1, t, acc1);
  mydbl adot = (mydbl)interp_der_scale_factor(spline2, t, acc2);
  mydbl f = - (2.0*adot*p0*pr)/(C*a) - (der_potential(r)*p0*p0)/(a*a*C*C) + r*(ptheta*ptheta + sinl(theta)*sinl(theta)*pphi*pphi) - (der_potential(r)/(C*C))*(-pr*pr + powl(r,2.0)*(powl(ptheta,2.0) + powl(sinl(theta)*pphi,2.0)));
  return f;
}

/*Function of the 2th (polar) momentum component differential equation for the geodesics.
${p2}^{dot} = f2(x^{\alpha},p^{\alpha})$.*/
mydbl geodesic_equation_theta(gsl_spline *spline1, gsl_interp_accel *acc1, gsl_spline *spline2, gsl_interp_accel *acc2, mydbl p0, mydbl pr, mydbl ptheta, mydbl pphi, mydbl x0, mydbl r, mydbl theta)
{
  double t = (double)(x0/C);
  mydbl a = (mydbl)interp_scale_factor(spline1, t, acc1);
  mydbl adot = (mydbl)interp_der_scale_factor(spline2, t, acc2);
  mydbl f = 0.5*sinl(2.0*theta)*pphi*pphi + 2.0*(der_potential(r)/powl(C,2.0) - 1/r)*pr*ptheta - (2.0*adot*p0*ptheta)/(C*a);
  return f;
}

/*Function of the 3th (azimuthal) momentum component differential equation for the geodesics.
${p3}^{dot} = f3(x^{\alpha},p^{\alpha})$.*/
mydbl geodesic_equation_phi(gsl_spline *spline1, gsl_interp_accel *acc1, gsl_spline *spline2, gsl_interp_accel *acc2, mydbl p0, mydbl pr, mydbl ptheta, mydbl pphi, mydbl x0, mydbl r, mydbl theta)
{
  double t = (double)(x0/C);
  mydbl a = (mydbl)interp_scale_factor(spline1, t, acc1);
  mydbl adot = (mydbl)interp_der_scale_factor(spline2, t, acc2);
  mydbl f = 2.0*( der_potential(r)/powl(C,2.0) - 1.0/r)*pr*pphi - 2.0*(1.0/tanl(theta))*ptheta*pphi - (2.0*adot*p0*pphi)/(a*C);
  return f;
}

/*Function for solving the geodesics differential equations using Euler's method.
Arguments are pointer so variables in that memory addresses are changed every time this function is called.*/
void runge_kutta_4(gsl_spline *spline1, gsl_interp_accel *acc1, gsl_spline *spline2, gsl_interp_accel *acc2, mydbl *x0, mydbl *x1, mydbl *x2, mydbl *x3, mydbl *p0, mydbl *p1, mydbl *p2, mydbl *p3, mydbl *lambda)
{
  /*Increment in the variables of the differential equation we want to solve*/
  mydbl dx0, dx1, dx2, dx3, dp0, dp1, dp2, dp3;

  /*dxi = (k1,j + 2*k2,j + 2*k3,j + k4,j)/6. In this sections the ki,j are declared with i=1,2,3,4.*/
  mydbl k1x0, k1x1, k1x2, k1x3, k1p0, k1p1, k1p2, k1p3;
  mydbl k2x0, k2x1, k2x2, k2x3, k2p0, k2p1, k2p2, k2p3;
  mydbl k3x0, k3x1, k3x2, k3x3, k3p0, k3p1, k3p2, k3p3;
  mydbl k4x0, k4x1, k4x2, k4x3, k4p0, k4p1, k4p2, k4p3;

  /*This section calculates the k1 quantities*/
  k1x0 = *p0*DLAMBDA; k1x1 = *p1*DLAMBDA; k1x2 = *p2*DLAMBDA; k3x0 = *p3*DLAMBDA;
  k1p0 = DLAMBDA*geodesic_equation_0(spline1, acc1, spline2, acc2, *p0, *p1, *p2, *p3, *x0, *x1, *x2);
  k1p1 = DLAMBDA*geodesic_equation_r(spline1, acc1, spline2, acc2, *p0, *p1, *p2, *p3, *x0, *x1, *x2);
  k1p2 = DLAMBDA*geodesic_equation_theta(spline1, acc1, spline2, acc2, *p0, *p1, *p2, *p3, *x0, *x1, *x2);
  k1p3 = DLAMBDA*geodesic_equation_phi(spline1, acc1, spline2, acc2, *p0, *p1, *p2, *p3, *x0, *x1, *x2);

  /*This section calculates the k2 quantities*/
  k2x0 = DLAMBDA*(*p0 + 0.5*k1x0); k2x1 = DLAMBDA*(*p1 + 0.5*k1x1); k2x2 = DLAMBDA*(*p2 + 0.5*k1x2); k2x3 = DLAMBDA*(*p3 + 0.5*k1x3);
  k2p0 = DLAMBDA*geodesic_equation_0(spline1, acc1, spline2, acc2, *p0 + 0.5*k1p0, *p1 + 0.5*k1p1, *p2 + 0.5*k1p2, *p3 + 0.5*k1p3, *x0 + 0.5*k1x0, *x1 + 0.5*k1x1, *x2 + 0.5*k1x2);
  k2p1 = DLAMBDA*geodesic_equation_r(spline1, acc1, spline2, acc2, *p0 + 0.5*k1p0, *p1 + 0.5*k1p1, *p2 + 0.5*k1p2, *p3 + 0.5*k1p3, *x0 + 0.5*k1x0, *x1 + 0.5*k1x1, *x2 + 0.5*k1x2);
  k2p2 = DLAMBDA*geodesic_equation_theta(spline1, acc1, spline2, acc2, *p0 + 0.5*k1p0, *p1 + 0.5*k1p1, *p2 + 0.5*k1p2, *p3 + 0.5*k1p3, *x0 + 0.5*k1x0, *x1 + 0.5*k1x1, *x2 + 0.5*k1x2);
  k2p3 = DLAMBDA*geodesic_equation_phi(spline1, acc1, spline2, acc2, *p0 + 0.5*k1p0, *p1 + 0.5*k1p1, *p2 + 0.5*k1p2, *p3 + 0.5*k1p3, *x0 + 0.5*k1x0, *x1 + 0.5*k1x1, *x2 + 0.5*k1x2);

  /*This section calculates the k3 quantities*/
  k3x0 = DLAMBDA*(*p0 + 0.5*k2x0); k3x1 = DLAMBDA*(*p1 + 0.5*k2x1); k3x2 = DLAMBDA*(*p2 + 0.5*k2x2); k3x3 = DLAMBDA*(*p3 + 0.5*k2x3);
  k3p0 = DLAMBDA*geodesic_equation_0(spline1, acc1, spline2, acc2, *p0 + 0.5*k2p0, *p1 + 0.5*k2p1, *p2 + 0.5*k2p2, *p3 + 0.5*k2p3, *x0 + 0.5*k2x0, *x1 + 0.5*k2x1, *x2 + 0.5*k2x2);
  k3p1 = DLAMBDA*geodesic_equation_r(spline1, acc1, spline2, acc2, *p0 + 0.5*k2p0, *p1 + 0.5*k2p1, *p2 + 0.5*k2p2, *p3 + 0.5*k2p3, *x0 + 0.5*k2x0, *x1 + 0.5*k2x1, *x2 + 0.5*k2x2);
  k3p2 = DLAMBDA*geodesic_equation_theta(spline1, acc1, spline2, acc2, *p0 + 0.5*k2p0, *p1 + 0.5*k2p1, *p2 + 0.5*k2p2, *p3 + 0.5*k2p3, *x0 + 0.5*k2x0, *x1 + 0.5*k2x1, *x2 + 0.5*k2x2);
  k3p3 = DLAMBDA*geodesic_equation_phi(spline1, acc1, spline2, acc2, *p0 + 0.5*k2p0, *p1 + 0.5*k2p1, *p2 + 0.5*k2p2, *p3 + 0.5*k2p3, *x0 + 0.5*k2x0, *x1 + 0.5*k2x1, *x2 + 0.5*k2x2);

  /*This section calculates the k4 quantities*/
  k4x0 = DLAMBDA*(*p0 + k3x0); k4x1 = DLAMBDA*(*p1 + k3x1); k4x2 = DLAMBDA*(*p2 + k3x2); k4x3 = DLAMBDA*(*p3 + k3x3);
  k4p0 = DLAMBDA*geodesic_equation_0(spline1, acc1, spline2, acc2, *p0 + 0.5*k3p0, *p1 + 0.5*k3p1, *p2 + 0.5*k3p2, *p3 + 0.5*k3p3, *x0 + 0.5*k3x0, *x1 + 0.5*k3x1, *x2 + 0.5*k3x2);
  k4p1 = DLAMBDA*geodesic_equation_r(spline1, acc1, spline2, acc2, *p0 + 0.5*k3p0, *p1 + 0.5*k3p1, *p2 + 0.5*k3p2, *p3 + 0.5*k3p3, *x0 + 0.5*k3x0, *x1 + 0.5*k3x1, *x2 + 0.5*k3x2);
  k4p2 = DLAMBDA*geodesic_equation_theta(spline1, acc1, spline2, acc2, *p0 + 0.5*k3p0, *p1 + 0.5*k3p1, *p2 + 0.5*k3p2, *p3 + 0.5*k3p3, *x0 + 0.5*k3x0, *x1 + 0.5*k3x1, *x2 + 0.5*k3x2);
  k4p3 = DLAMBDA*geodesic_equation_phi(spline1, acc1, spline2, acc2, *p0 + 0.5*k3p0, *p1 + 0.5*k3p1, *p2 + 0.5*k3p2, *p3 + 0.5*k3p3, *x0 + 0.5*k3x0, *x1 + 0.5*k3x1, *x2 + 0.5*k3x2);

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

/*To set the initial value of pr, it must hold $g_{\mu\nu}p^{\mu}p^{\nu} = 0$.
This factor multiplies p0 to guarantee that p1 fulfill the null geodesic condition.*/
mydbl condition_factor(mydbl r)
{
  return 1.0+2.0*potential(r)/(C*C);
}

/*$cp^{0}$ multiplied by this factor allows to obtain the energy for a local inertial observer in this spacetime.*/
mydbl energy_factor(mydbl r)
{
  mydbl g = 1.0 + potential(r)/(C*C);
  return g;
}

/*Violation of null geodesics condition $g_{\mu\nu}p^{\mu}p^{\nu} = 0$.*/
mydbl violation(mydbl r, mydbl theta, mydbl phi, mydbl p0, mydbl pr, mydbl ptheta, mydbl pphi)
{
  mydbl f = -(1.0+2.0*potential(r)/(C*C))*p0*p0 + (1.0-2.0*potential(r)/(C*C))*(pr*pr + r*r*(powl(sinl(theta),2.0)*pphi*pphi + ptheta*ptheta));
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
  gsl_spline_init(spline1, cosmictime, scale_factor, NLINESFRW);  //Initializes spline object for data conf_time, scale_factor of size NLINES
  gsl_spline_init(spline2, cosmictime, der_scale_factor, NLINESFRW);  //Initializes spline object for data conf_time, der_scale_factor of size NLINES

  /************************************************************************************/
  
  /***SOLVES GEODESIC EQUATIONS FOR PERTURBED FRW UNIVERSE WITH STATIC PLUMMER POTENTIAL ***/

  /*Initial conditions*/
  mydbl ti = 7.0 ,x0, r = -3500.0, theta = M_PI*0.5, phi = M_PI*0.5, p0 = 1.0e-3, pr, ptheta = 0.0, pphi = 0.0, lambda = 0.0, energy1, energy, v, difft;
  x0 = C*ti;
  pr = condition_factor(r)*p0;
  energy1 = C*energy_factor(r)*p0;
  v = violation(r, theta, phi, p0, pr, ptheta, pphi);
  difft = (energy1 - energy1)/energy1;

  /*Pointer to file where solution of differential equation will be saved.*/
  FILE *geodesic;
  geodesic = fopen("geodesic_solution.dat","w");

  /*Write line of initial values in file*/
  fprintf(geodesic, "%.18Lf %.18Lf %.18Lf %.18Lf %.18Lf %.18Lf %.18Lf %.18Lf %.18Lf %.18Lf %.18Lf %.18Lf\n", lambda, x0, r, theta, phi, p0, pr, ptheta, pphi, energy1, v, difft);

  /*Solution of the differential equation*/
  for(i=0; i<(1+ NLINES); i++)
    {
      runge_kutta_4(spline1, acc1, spline2, acc2, &x0, &r, &theta, &phi, &p0, &pr, &ptheta, &pphi, &lambda);
      energy = C*energy_factor(r)*p0;
      v = violation(r, theta, phi, p0, pr, ptheta, pphi);
      difft = (energy - energy1)/energy1;
      fprintf(geodesic, "%.18Lf %.18Lf %.18Lf %.18Lf %.18Lf %.18Lf %.18Lf %.18Lf %.18Lf %.18Lf %.18Lf %.18Lf\n", lambda, x0, r, theta, phi, p0, pr, ptheta, pphi, energy, v, difft);
    }

  /************************************************************************************/

  /* FILE *graph; */
  /* graph = fopen("graph.dat","w"); */

  /* /\*Graphication of functions depending on time*\/ */
  /* for(eta = 0.12; eta <= 3.29; eta = eta + 0.01) */
  /*   { */
  /*     double a = interp_scale_factor(spline1, eta, acc1); */
  /*     double aprime = interp_der_scale_factor(spline2, eta, acc2); */
  /*     fprintf(graph, "%.12lf %.12lf %.12lf %.12lf\n", eta, a, aprime, aprime/a); */
  /*   } */

  /* fclose(graph); */

  /* /\*Graphication of functions depending on position*\/ */
  /* double x1; */
  /* FILE *graph2; */
  /* graph2 = fopen("graph2.dat", "w"); */
  
  /* for(x1 = 0.0; x1 < 100.0; x1 = x1 + 1.0) */
  /*   { */
  /*     fprintf(graph2, "%.12lf %.12lf %.12lf\n", x1, potential(x1,0.0,0.0), ith_der_potential(x1,0.0,0.0)); */
  /*   } */
  
  /* fclose(graph2); */

  /* /\*GNUPLOT*\/ */

  /* FILE *plot; */
  /* plot = popen("gnuplot -persist","w"); */
  /* fprintf(plot, "set terminal x11 0\n"); */
  /* fprintf(plot, "set multiplot layout 1,3\n"); */
  /* fprintf(plot, "plot 'graph.dat' using 1:2 not\n"); */
  /* fprintf(plot, "plot 'graph.dat' using 1:3 not\n"); */
  /* fprintf(plot, "plot 'graph.dat' using 1:4 not\n"); */
  /* fprintf(plot, "unset multiplot\n"); */
  /* fprintf(plot, "reset\n"); */
  /* fprintf(plot, "set terminal x11 1\n"); */
  /* fprintf(plot, "set multiplot layout 1,2\n"); */
  /* fprintf(plot, "plot 'graph2.dat' using 1:2 not\n"); */
  /* fprintf(plot, "plot 'graph2.dat' using 1:3 not\n"); */
  /* fprintf(plot, "unset multiplot\n"); */
  /* //fprintf(plot, "plot 'graph2.dat' using 1:4 not\n"); */

  /* pclose(plot); */
  /* system("rm graph.dat"); */
  /* system("rm graph2.dat"); */

  

  /*** Releasing all used space in memory ***/
  fclose(geodesic); //Close file storing the results
  gsl_spline_free(spline1);  //Free memory of spline object
  gsl_spline_free(spline2);
  gsl_interp_accel_free(acc1);  //Free memory of accel object
  gsl_interp_accel_free(acc2);
}
