/*This program evaluates the differential equations for a photon's geodesic in a perturbed spacetime.
The equations are written in the form $\frac{d(x or p)^{\alpha}}{d\lambda}=f(x^{\alpha},p^{\alpha})$.
Where $p^{\alpha}={\dot{x}}^{\alpha}$*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_errno.h>          //GSL error management module
#include <gsl/gsl_spline.h>         //GSL interpolation module

#define A 50.0     //Distance parameter of the perturbations
#define G 43007.01     //Gravitational constant
#define M 50.0     //Mass of the perturbation
#define NLINES 999 //Number of lines in frw.dat file
#define DLAMBDA 0.001   //Geodesics parameter step

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

/*Function for the gravitational potential to be used*/
double potential(double x1, double x2, double x3)
{
  return -G*M/(sqrt(A*A + x1*x1 + x2*x2 +x3*x3));
}

/*First partial derivative respect to the ith coordinate.
First argument xi denotes the ith coordinate.
The xj's denote the other coordinates. The order doesn't matter since they appear in the same way in the equation.*/
double ith_der_potential(double xi, double xj1, double xj2)
{
  return G*M*xi/(pow(A*A+xi*xi+xj1*xj1+xj2*xj2, 1.5));
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
double geodesic_equation_i(double aprime, double a, double p0, double pi, double pj1, double pj2, double eta, double xi, double xj1, double xj2)
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
  dp0 = geodesic_equation_0(aprime, a, *p0, *p1, *p2, *p3, *eta, *x1, *x2, *x3)*DLAMBDA;
  dp1 = geodesic_equation_i(aprime, a, *p0, *p1, *p2, *p3, *eta, *x1, *x2, *x3)*DLAMBDA;
  dp2 = geodesic_equation_i(aprime, a, *p0, *p2, *p1, *p3, *eta, *x2, *x1, *x3)*DLAMBDA;
  dp3 = geodesic_equation_i(aprime, a, *p0, *p3, *p2, *p1, *eta, *x3, *x2, *x1)*DLAMBDA;

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

  /*Free space in memory*/
  fclose(frw);
  
  /*** Initializes objects for interpolation. 1 is for interpolation of scale factor, 2 is for interpolation of derivative of scale factor ***/

  /*Allocate space in memory*/
  gsl_interp_accel *acc1 = gsl_interp_accel_alloc();  //Acceleration type object (for index lookup)
  gsl_interp_accel *acc2 = gsl_interp_accel_alloc();
  gsl_spline *spline1 = gsl_spline_alloc(gsl_interp_cspline, NLINES);  //Spline type object (define interpolation type and space in memory, works for both)
  gsl_spline *spline2 = gsl_spline_alloc(gsl_interp_cspline, NLINES);

  /*Initializes objects for interpolation*/
  gsl_spline_init(spline1, conftime, scale_factor, NLINES);  //Initializes spline object for data conf_time, scale_factor of size NLINES
  gsl_spline_init(spline2, conftime, der_scale_factor, NLINES);  //Initializes spline object for data conf_time, der_scale_factor of size NLINES

  /************************************************************************************/

  /*Initial conditions*/
  double eta = 0.12; // x1 = 0, x2 = 0, x3 = 0, p0 = 0, p1 = 0, p2 = 0, p3 = 0, lambda = 0;

  /*Initial values for scale factor and derivative of scale factor*/
  //double a = interp_scale_factor(spline1, eta, acc1);
  //double aprime = interp_der_scale_factor(spline2, eta, acc2);

  FILE *graph;
  graph = fopen("graph.dat","w");

  /*Graphication of functions depending on time*/
  for(eta = 0.12; eta <= 3.29; eta = eta + 0.01)
    {
      double a = interp_scale_factor(spline1, eta, acc1);
      double aprime = interp_der_scale_factor(spline2, eta, acc2);
      fprintf(graph, "%.12lf %.12lf %.12lf %.12lf\n", eta, a, aprime, aprime/a);
    }

  fclose(graph);

  /*Graphication of functions depending on position*/
  double x1;
  FILE *graph2;
  graph2 = fopen("graph2.dat", "w");
  
  for(x1 = 0.0; x1 < 100.0; x1 = x1 + 1.0)
    {
      fprintf(graph2, "%.12lf %.12lf %.12lf\n", x1, potential(x1,0.0,0.0), ith_der_potential(x1,0.0,0.0));
    }
  
  fclose(graph2);

  /*GNUPLOT*/

  FILE *plot;
  plot = popen("gnuplot -persist","w");
  fprintf(plot, "set terminal x11 0\n");
  fprintf(plot, "set multiplot layout 1,3\n");
  fprintf(plot, "plot 'graph.dat' using 1:2 not\n");
  fprintf(plot, "plot 'graph.dat' using 1:3 not\n");
  fprintf(plot, "plot 'graph.dat' using 1:4 not\n");
  fprintf(plot, "unset multiplot\n");
  fprintf(plot, "reset\n");
  fprintf(plot, "set terminal x11 1\n");
  fprintf(plot, "set multiplot layout 1,2\n");
  fprintf(plot, "plot 'graph2.dat' using 1:2 not\n");
  fprintf(plot, "plot 'graph2.dat' using 1:3 not\n");
  fprintf(plot, "unset multiplot\n");
  //fprintf(plot, "plot 'graph2.dat' using 1:4 not\n");

  pclose(plot);
  system("rm graph.dat");
  system("rm graph2.dat");

  

  

  /*Pointer to file where solution of differential equation will be saved.*/
  //FILE *geodesic;
  //geodesic = fopen("geodesic_solution.dat","w");

  /*Write line of initial values in file*/
  //fprintf(geodesic, "%12lf %12lf %12lf %12lf %12lf %12lf %12lf %12lf %12lf", lambda, eta, x1, x2, x3, p0, p1, p2, p3);

  /*Solution of the differential equation*/
  //for(i=0; i<1000; i++)
  //  {
  //    euler1(a, aprime, &eta, &x1, &x2, &x3, &p0, &p1, &p2, &p3, &lambda);
  //    a = interp_scale_factor(conftime, scale_factor, eta);
  //    aprime = interp_der_scale_factor(conftime, der_scale_factor, eta);
  //  }

  /** Releasing all used space in memory **/
  //fclose(geodesic); //Close file storing the results
  gsl_spline_free(spline1);  //Free memory of spline object
  gsl_spline_free(spline2);
  gsl_interp_accel_free(acc1);  //Free memory of accel object
  gsl_interp_accel_free(acc2);
}
