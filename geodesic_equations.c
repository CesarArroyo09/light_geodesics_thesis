/*This program evaluates the differential equations for a photon's geodesic in a perturbed spacetime.
The equations are written in the form $\frac{d(x or p)^{\alpha}}{d\lambda}=f(x^{\alpha},p^{\alpha})$.
Where $p^{\alpha}={\dot{x}}^{\alpha}$*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_errno.h>          //GSL error management module
#include <gsl/gsl_spline.h>         //GSL interpolation module

#define A 1.0     //Distance parameter of the perturbations
#define G 1.0     //Gravitational constant
#define M 1.0     //Mass of the perturbation
#define NLINES 999 //Number of lines in frw.dat file
#define DLAMBDA 0.001   //Geodesics parameter step

/*Interpolation of scale factor at time eta.
Argument conf_time is an array with the conformal times. scale_factor is an array of scale factors which corresponds in array position to the conformal time array.*/
double interp_scale_factor(double conf_time[], double scale_factor[], double eta)
{
  /*Allocate space in memory*/
  gsl_interp_accel *acc = gsl_interp_accel_alloc();  //Acceleration type object (for index lookup)
  gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline, NLINES);  //Spline type object (define interpolation type and space in memory)

  gsl_spline_init(spline, conf_time, scale_factor, NLINES);  //Initializes spline object for data conf_time, scale_factor of size NLINES
  double a = gsl_spline_eval(spline, eta, acc);  //Interpolates data to abcisa eta using method in spline and acceleration object acc
  
  /*Free space in memory*/
  gsl_spline_free(spline);  //Free memory of spline object
  gsl_interp_accel_free(acc);  //Free memory of accel object

  return a;  //Return scale factor at a time eta
}

/*Interpolation of derivative of scale factor at time eta.
This function works the same way as interp_scale_factor.*/
double interp_der_scale_factor(double conf_time[], double der_scale_factor[], double eta)
{
  /*Allocate space in memory*/
  gsl_interp_accel *acc = gsl_interp_accel_alloc();  //Acceleration type object (for index lookup)
  gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline, NLINES);  //Spline type object (define interpolation type and space in memory)

  gsl_spline_init(spline, conf_time, der_scale_factor, NLINES);  //Initializes spline object for data conf_time, der_scale_factor of size NLINES
  double adot = gsl_spline_eval(spline, eta, acc);  //Interpolates data to abcisa eta using method in spline and acceleration object acc
  
  /*Free space in memory*/
  gsl_spline_free(spline);  //Free memory of spline object
  gsl_interp_accel_free(acc);  //Free memory of accel object

  return adot;  //Return derivative of scale factor at a time eta
}

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
double geodesic_equation_0(double aprime, double a, double p0, double p1, double p2, double p3, double eta, double x1, double x2, double x3)
{
  double f = - (aprime/a)*pow(p0,2) - 2*p0*(ith_der_potential(x1,x2,x3)*p1 + ith_der_potential(x2,x1,x3)*p2 + ith_der_potential(x3,x1,x2)*p3) - (aprime/a)*(1-4*potential(x1,x2,x3))*(pow(p1,2) + pow(p2,2) + pow(p3,2));
  return f;
}

/*Function for the ith differential equation for the geodesics, where i=1,2,3.
xi is the ith coordinate, pi is the ith momentum.
The other position and momentum quantities are denoted xj1, xj2, pj1 and pj2. These quantities appear in a symmetric way in the equation so in doesn't matter the order.*/
double geodesic_equation_i(double aprime, double a, double p0, double pi, double pj1, double pj2, double pj3, double eta, double xi, double xj1, double xj2)
{
  double f =  -ith_der_potential(xi,xj1,xj2)*pow(p0,2) - 2*(aprime/a)*p0*pi + ith_der_potential(xi,xj1,xj2)*(pow(pi,2)-pow(pj1,2)-pow(pj2,2)) + 2*ith_der_potential(xj1,xi,xj2)*pi*pj1 + 2*ith_der_potential(xj2,xi,xj1)*pi*pj2;
  return f;
}

void euler1(double a, double aprime, double *eta, double *x1, double *x2, double *x3, double *p0, double *p1, double *p2, double *p3, double *lambda)
{
  /*Increment in the variables of the differential equation we want to solve*/
  double deta, dx1, dx2, dx3, dp0, dp1, dp2, dp3;

  /*Calculation of the increments*/
  deta = *p0*DLAMBDA; dx1 = *p1*DLAMBDA; dx2 = *p2*DLAMBDA; dx3 = *p3*DLAMBDA;
  dp0 = geodesic_equation_0(aprime, a, *p0, *p1, *p2, *p3, *eta, *x1, *x2, *x3, *x4)*DLAMBDA;
  dp1 = geodesic_equation_i(aprime, a, *p0, *p1, *p2 *p3, *eta, *x1, *x2, *x3)*DLAMBDA;
  dp2 = geodesic_equation_i(aprime, a, *p0, *p2, *p1 *p3, *eta, *x2, *x1, *x3)*DLAMBDA;
  dp3 = geodesic_equation_i(aprime, a, *p0, *p3, *p2 *p1, *eta, *x3, *x2, *x1)*DLAMBDA;

  /*New values of the variables of the differential equation. Since we are using pointers, when called the routine the value of variable change.*/
  *eta = *eta + deta; *x1 = *x1 + dx1; *x2 = *x2 + dx2; *x3 = *x3 + dx3;
  *p0 = *p0 + dp0; *p1 = *p1 + dp1; *p2 = *p2 + dp2; *p3 = *p3 + dp3;
  
  /*Increment of parameter of geodesics*/
  *lambda = *lambda + DLAMBDA;
}

int main(void)
{
  int i;            //For array manipulation
  
  /*Pointer to frw.data file*/
  FILE *frw;        
  frw = fopen("frw.dat","r");

  /*Variables and arrays to read the data*/
  double cosmictime, conftime[NLINES], scale_factor[NLINES], der_scale_factor[NLINES];

  /*Reading the data*/
  for(i=0; i<NLINES; i++)
    {
      fscanf(frw,"%lf %lf %lf %lf", &cosmictime, &conftime[i], &scale_factor[i], &der_scale_factor[i]);
    }
  
  /*Initial conditions*/
  double eta = 0, x1 = 0, x2 = 0, x3 = 0, p0 = 0, p1 = 0, p2 = 0, p3 = 0, lambda = 0;

  /*Initial values for scale factor and derivative of scale factor*/
  double a = interp_scale_factor(conftime, scale_factor, eta);
  double aprime = interp_der_scale_factor(conftime, der_scale_factor, eta);

  /*Pointer to file where solution of differential equation will be saved.*/
  FILE *geodesic;
  geodesic = fopen("geodesic_solution.dat","w");

  /*Write line of initial values in file*/
  fprintf(geodesic, "%12lf %12lf %12lf %12lf %12lf %12lf %12lf %12lf %12lf", lambda, eta, x1, x2, x3, p0, p1, p2, p3);

  /*Solution of the differential equation*/
  for(i=0; i<1000; i++)
    {
      euler1(a, aprime, &eta, &x1, &x2, &x3, &p0, &p1, &p2, &p3, &lambda);
      a = interp_scale_factor(conftime, scale_factor, eta);
      aprime = interp_der_scale_factor(conftime, der_scale_factor, eta);
    }
}
