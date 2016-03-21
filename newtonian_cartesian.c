/*This program evaluates the differential equations for a photon's geodesic in a perturbed spacetime.
This program solves the particular case for the Minkowski perturbed spacetime with metric: g_{ab} = {\eta}_{ab} + h_{ab}. Where h_{ab} i
The equations are written in the form $\frac{d(x or p)^{\alpha}}{d\lambda}=f(x^{\alpha},p^{\alpha})$.
Where $p^{\alpha}={\dot{x}}^{\alpha}$*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define A 1.0     //Distance parameter of the perturbations
#define G 43007.01     //Gravitational constant
#define M 0.0     //Mass of the perturbation
#define C 300000.0  //Speed of light
#define NLINES 999 //Number of lines in frw.dat file
#define DLAMBDA 0.000001   //Geodesics parameter step

typedef long double mydbl;

/*Function for the gravitational potential to be used*/
mydbl potential(mydbl x1, mydbl x2, mydbl x3)
{
  return -G*M/(sqrtl(A*A + x1*x1 + x2*x2 +x3*x3));
}

/*First partial derivative respect to the ith coordinate.
First argument xi denotes the ith coordinate.
The xj's denote the other coordinates. The order doesn't matter since they appear in the same way in the equation.*/
mydbl ith_der_potential(mydbl xi, mydbl xj1, mydbl xj2)
{
  return G*M*xi/(powl(A*A+xi*xi+xj1*xj1+xj2*xj2, 1.5));
}

/*This is the function of the first differential equation for the geodesics.
l^{dot} = f1(variables of the problem)*/
mydbl geodesic_equation_0(mydbl p0, mydbl p1, mydbl p2, mydbl p3, mydbl x0, mydbl x1, mydbl x2, mydbl x3)
{
  mydbl f = - p0*(p1*ith_der_potential(x1, x2, x3) + p2*ith_der_potential(x2, x1, x3) + p3*ith_der_potential(x3, x2, x1));
  return f;
}

/*Function for the ith differential equation for the geodesics, where i=1,2,3.
xi is the ith coordinate, pi is the ith momentum.
The other position and momentum quantities are denoted xj1, xj2, pj1 and pj2. These quantities appear in a symmetric way in the equation so in doesn't matter the order.*/
mydbl geodesic_equation_i(mydbl p0, mydbl pi, mydbl pj1, mydbl pj2, mydbl x0, mydbl xi, mydbl xj1, mydbl xj2)
{
  mydbl f =  (1/powl(C,2)) * ( 2*pi * (pi*ith_der_potential(xi,xj1,xj2) + pj1*ith_der_potential(xj1,xi,xj2) + pj2*ith_der_potential(xj2,xj1,xi)) - ith_der_potential(xi,xj1,xj2) * (powl(p0,2) + powl(pi,2) + powl(pj1,2) + powl(pj2,2)) ) ;
  return f;
}

void euler1(mydbl *x0, mydbl *x1, mydbl *x2, mydbl *x3, mydbl *p0, mydbl *p1, mydbl *p2, mydbl *p3, mydbl *lambda)
{
  /*Increment in the variables of the differential equation we want to solve*/
  mydbl dx0, dx1, dx2, dx3, dp0, dp1, dp2, dp3;

  /*Calculation of the increments*/
  dx0 = *p0*DLAMBDA; dx1 = *p1*DLAMBDA; dx2 = *p2*DLAMBDA; dx3 = *p3*DLAMBDA;
  dp0 = geodesic_equation_0(*p0, *p1, *p2, *p3, *x0, *x1, *x2, *x3)*DLAMBDA;
  dp1 = geodesic_equation_i(*p0, *p1, *p2, *p3, *x0, *x1, *x2, *x3)*DLAMBDA;
  dp2 = geodesic_equation_i(*p0, *p2, *p1, *p3, *x0, *x2, *x1, *x3)*DLAMBDA;
  dp3 = geodesic_equation_i(*p0, *p3, *p2, *p1, *x0, *x3, *x2, *x1)*DLAMBDA;

  /*New values of the variables of the differential equation. Since we are using pointers, when called the routine the value of variable change.*/
  *x0 = *x0 + dx0; *x1 = *x1 + dx1; *x2 = *x2 + dx2; *x3 = *x3 + dx3;
  *p0 = *p0 + dp0; *p1 = *p1 + dp1; *p2 = *p2 + dp2; *p3 = *p3 + dp3;
  
  /*Increment of parameter of geodesics*/
  *lambda = *lambda + DLAMBDA;
}

mydbl g00(mydbl x1, mydbl x2, mydbl x3)
{
  mydbl g = 1.0 + 2*potential(x1,x2,x3)/(C*C);
  return g;
}

mydbl gii(mydbl x1, mydbl x2, mydbl x3)
{
  mydbl g = 1.0 - 2*potential(x1,x2,x3)/(C*C);
  return g;
}

int main(void)
{
  int i;            //For array manipulation

  /* /\*Graphication of functions depending on position*\/ */
  /* mydbl x1; */
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
  mydbl t = 0.0, x1 = 0.0, x2 = 0.0, x3 = 0.0, p0 = 0.001, p1 = p0, p2 = 0.0, p3 = 0.0, lambda = 0.0, nu = 0.0;
  //p1 = -sqrtl(g00(x1,x2,x3)/gii(x1,x2,x3))*p0;

  /*Pointer to file where solution of differential equation will be saved.*/
  FILE *geodesic;
  geodesic = fopen("geodesic_solution.dat","w");

  /*Write line of initial values in file*/
  fprintf(geodesic,"%.18Lf %.18Lf %.18Lf %.18Lf %.18Lf %.18Lf %.18Lf %.18Lf %.18Lf\n", lambda, t, x1, x2, x3, p0, p1, p2, p3);

  /*Solution of the differential equation*/
  for(i=0; i<100000; i++)
    {
      euler1(&t, &x1, &x2, &x3, &p0, &p1, &p2, &p3, &lambda);
      fprintf(geodesic,"%.18Lf %.18Lf %.18Lf %.18Lf %.18Lf %.18Lf %.18Lf %.18Lf %.18Lf\n", lambda, t, x1, x2, x3, p0, p1, p2, p3);
    }

  /** Releasing all used space in memory **/
  fclose(geodesic); //Close file storing the results
}
