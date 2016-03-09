/*This program evaluates the differential equations for a photon's geodesic in a perturbed spacetime.
This program solves the particular case for the Minkowski perturbed spacetime with metric: g_{ab} = {\eta}_{ab} + h_{ab}. Where h_{ab} i
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
#define C 299792.458  //Speed of light
#define NLINES 999 //Number of lines in frw.dat file
#define DLAMBDA 0.001   //Geodesics parameter step


/*Function for the gravitational potential to be used*/
double potential(double r)
{
  return -G*M/(sqrt(A*A + r*r));
}

/*First partial derivative respect to the ith coordinate.
First argument xi denotes the ith coordinate.
The xj's denote the other coordinates. The order doesn't matter since they appear in the same way in the equation.*/
double der_potential(double r)
{
  return G*M*r/(pow(A*A+r*r, 1.5));
}

/*This is the function of the first differential equation for the geodesics.
l^{dot} = f1(variables of the problem)*/
double geodesic_equation_0(double p0, double pr, double r)
{
  double f = -(2/pow(C,2))*(1-(2*potential(r)/pow(C,2)))*der_potential(r)*p0*pr;
  return f;
}

/*Function for the ith differential equation for the geodesics, where i=1,2,3.
xi is the ith coordinate, pi is the ith momentum.
The other position and momentum quantities are denoted xj1, xj2, pj1 and pj2. These quantities appear in a symmetric way in the equation so in doesn't matter the order.*/
double geodesic_equation_r(double p0, double pr, double ptheta, double pphi, double r, double theta)
{
  double f =  (2/pow(C,2))*(1+(2*potential(r)/pow(C,2)))*der_potential(r)*pr*pr - der_potential(r)*(1+(2*potential(r)/pow(C,2)))*p0*p0 + r*(pow(sin(theta),2)*pphi*pphi + ptheta*ptheta);
  return f;
}

double geodesic_equation_theta(double pr, double ptheta, double pphi, double r, double theta)
{
  double f = 2*(der_potential(r)/pow(C,2) - 1/r)*pr*ptheta + 0.5*sin(2*theta)*pphi*pphi;
  return f;
}

double geodesic_equation_phi(double pr, double ptheta, double pphi, double r, double theta)
{
  double f = 2*(der_potential(r)/pow(C,2) - 1/r)*pr*pphi + 2*(1/tan(theta))*ptheta*pphi;
  return f;
}

void euler1(double *t, double *r, double *theta, double *phi, double *p0, double *pr, double *ptheta, double *pphi, double *lambda)
{
  /*Increment in the variables of the differential equation we want to solve*/
  double dt, dr, dtheta, dphi, dp0, dpr, dptheta, dpphi;

  /*Calculation of the increments*/
  dt = *p0*DLAMBDA; dr = *pr*DLAMBDA; dtheta = *theta*DLAMBDA; dphi = *phi*DLAMBDA;
  dp0 = geodesic_equation_0(*p0, *pr, *r)*DLAMBDA;
  dpr = geodesic_equation_r(*p0, *pr, *ptheta, *pphi, *r, *theta)*DLAMBDA;
  dptheta = geodesic_equation_theta(*pr, *ptheta, *pphi, *r, *theta)*DLAMBDA;
  dpphi = geodesic_equation_phi(*pr, *ptheta, *pphi, *r, *theta)*DLAMBDA;

  /*New values of the variables of the differential equation. Since we are using pointers, when called the routine the value of variable change.*/
  *t = *t + dt; *r = *r + r; *theta = *theta + dtheta; *phi = *phi + dphi;
  *p0 = *p0 + dp0; *pr = *pr + dpr; *ptheta = *ptheta + dptheta; *pphi = *pphi + dpphi;
  
  /*Increment of parameter of geodesics*/
  *lambda = *lambda + DLAMBDA;
}

int main(void)
{
  int i;            //For array manipulation

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
  
  
  /************************************************************************************/

  
  /*Initial conditions*/
  double t = 0.0, r = 500, theta = 0.0, phi = 0.0, p0 = 0.0, pr = 0.0, ptheta = 0.0, pphi = 0.0, lambda = 0.0;

  /*Pointer to file where solution of differential equation will be saved.*/
  FILE *geodesic;
  geodesic = fopen("geodesic_solution.dat","w");

  /*Write line of initial values in file*/
  fprintf(geodesic, "%.12lf %.12lf %.12lf %.12lf %.12lf %.12lf %.12lf %.12lf %.12lf", lambda, t, r, theta, phi, p0, pr, ptheta, pphi);

  /*Solution of the differential equation*/
  for(i=0; i<1000; i++)
    {
      euler1(&t, &r, &theta, &phi, &p0, &pr, &ptheta, &pphi, &lambda);
      fprintf(geodesic, "%.12lf %.12lf %.12lf %.12lf %.12lf %.12lf %.12lf %.12lf %.12lf", lambda, t, r, theta, phi, p0, pr, ptheta, pphi);
    }

  /** Releasing all used space in memory **/
  fclose(geodesic); //Close file storing the results
}
