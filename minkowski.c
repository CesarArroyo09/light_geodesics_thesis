/*This program evaluates the differential equations for a photon's geodesic in Minkowski spacetime.
The equations are written in the form $\frac{d(x or p)^{\alpha}}{d\lambda}=f(x^{\alpha},p^{\alpha})$.
Where $p^{\alpha}={\dot{x}}^{\alpha}$*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_errno.h>          //GSL error management module

#define C 299792.458     //Speed of light in velocity units
#define DLAMBDA 0.001   //Geodesics parameter step


/*Function for solving the differential equations of Minkowski geodesics.*/
void euler1(double *x0, double *x1, double *x2, double *x3, double *p0, double *p1, double *p2, double *p3, double *lambda)
{
  /*Increment in the variables of the differential equation we want to solve*/
  double dx0, dx1, dx2, dx3, dp0, dp1, dp2, dp3;

  /*Calculation of the increments*/
  dx0 = *p0*DLAMBDA; dx1 = *p1*DLAMBDA; dx2 = *p2*DLAMBDA; dx3 = *p3*DLAMBDA;
  dp0 = 0.0; dp1 = 0.0; dp2 = 0.0; dp3 = 0.0;

  /*New values of the variables of the differential equation. Since we are using pointers, when called the routine the value of variable change.*/
  *x0 = *x0 + dx0; *x1 = *x1 + dx1; *x2 = *x2 + dx2; *x3 = *x3 + dx3;
  *p0 = *p0 + dp0; *p1 = *p1 + dp1; *p2 = *p2 + dp2; *p3 = *p3 + dp3;
  
  /*Increment of parameter of geodesics*/
  *lambda = *lambda + DLAMBDA;
}

int main(void)
{
  int i;            //For array manipulation

  /*Initial conditions*/
  double x0 = 0.0, x1 = 0.0, x2 = 0.0, x3 = 0.0, p1 = 0.1*C, p2 = 0.9*C, p3 = 0.0, lambda = 0.0;
  double p0 = sqrt(pow(p1,2) + pow(p2,2) + pow(p3,2))/C;

  /*Pointer to file where solution of differential equation will be saved.*/
  FILE *geodesic;
  geodesic = fopen("geodesic_solution.dat","w");
  
  /*Write line of initial values in file*/
  fprintf(geodesic, "%.12lf %.12lf %.12lf %.12lf %.12lf %.12lf %.12lf %.12lf %.12lf\n", lambda, x0, x1, x2, x3, p0, p1, p2, p3);

  /*Solution of the differential equation*/
  for(i=0; i<1000; i++)
    {
      euler1(&x0, &x1, &x2, &x3, &p0, &p1, &p2, &p3, &lambda);
      fprintf(geodesic, "%12lf %12lf %12lf %12lf %12lf %12lf %12lf %12lf %12lf\n", lambda, x0, x1, x2, x3, p0, p1, p2, p3);
    }

  /** Releasing all used space in memory **/
  fclose(geodesic); //Close file storing the results

  /*GNUPLOT*/

  FILE *plot;
  plot = popen("gnuplot -persist","w");
  fprintf(plot, "set terminal x11 0\n");
  fprintf(plot, "set multiplot layout 1,3\n");
  fprintf(plot, "plot 'geodesic_solution.dat' using 1:2 not\n");
  fprintf(plot, "plot 'geodesic_solution.dat' using 1:3 not\n");
  fprintf(plot, "plot 'geodesic_solution.dat' using 1:4 not\n");
  fprintf(plot, "unset multiplot\n");
  fprintf(plot, "reset\n");
  fprintf(plot, "set terminal x11 1\n");
  fprintf(plot, "set multiplot layout 1,3\n");
  fprintf(plot, "plot 'geodesic_solution.dat' using 1:6 not\n");
  fprintf(plot, "plot 'geodesic_solution.dat' using 1:7 not\n");
  fprintf(plot, "plot 'geodesic_solution.dat' using 1:8 not\n");
  fprintf(plot, "unset multiplot\n");
  fprintf(plot, "set terminal x11 2\n");
  fprintf(plot, "splot 'geodesic_solution.dat' using 2:3:1\n");

  pclose(plot);
  //system("rm geodesic_solution.dat");
}
