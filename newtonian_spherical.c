/*This program evaluates the differential equations for a photon's geodesic in a perturbed spacetime in SPHERICAL coordinates.
This program solves the particular case for the Minkowski perturbed spacetime with metric: $g_{ab} = {\eta}_{ab} + h_{ab}$. Where $h_{ab}$ corresponds to the perturbation in the Newtonian weak field limit. A Plummer potential with adequate parameters have been used to simulate the perturbation.
The equations are written in the form $\frac{d(x or p)^{\alpha}}{d\lambda}=f(x^{\alpha},p^{\alpha})$.
Where $p^{\alpha}={\dot{x}}^{\alpha}$ and the indice $\alpha$ runs from 0 to 3.
The coordinates for the photon's geodesics are: (ct,r,\theta,\phi) = (x0,x1,x2,x3).*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define A 1.0     //Distance parameter of the perturbations
#define G 43007.01     //Gravitational constant
#define M 0.0     //Mass of the perturbation
#define C 299792.458  //Speed of light
#define NLINES 99999 //Number of lines in frw.dat file
#define DLAMBDA 0.00001   //Geodesics parameter step

typedef long double mydbl;

/*Function for the gravitational potential to be used. Potential for Plummer model.*/
mydbl potential(mydbl r)
{
  return -G*M/(sqrt(A*A + r*r));
}

/*Derivative of potential respecto to radial coordinate.*/
mydbl der_potential(mydbl r)
{
  return G*M*r/(powl(A*A+r*r, 1.5));
}

/*Function of the 0th momentum component differential equation for the geodesics.
${p0}^{dot} = f0(x^{\alpha},p^{\alpha})$.*/
mydbl geodesic_equation_0(mydbl p0, mydbl pr, mydbl r)
{
  mydbl f = -(2/powl(C,2))*der_potential(r)*p0*pr;
  return f;
}

/*Function of the 1th (radial) momentum component differential equation for the geodesics.
${p1}^{dot} = f1(x^{\alpha},p^{\alpha})$.*/
mydbl geodesic_equation_r(mydbl p0, mydbl pr, mydbl ptheta, mydbl pphi, mydbl r, mydbl theta)
{
  mydbl f =  r*(powl(ptheta,2) + powl(sinl(theta),2)*powl(pphi,2)) - (der_potential(r)/powl(C,2))*(powl(p0,2) -powl(pr,2) + powl(r,2)*(powl(ptheta,2) + powl(sinl(theta),2)*powl(pphi,2)));
  return f;
}

/*Function of the 2th (polar) momentum component differential equation for the geodesics.
${p2}^{dot} = f2(x^{\alpha},p^{\alpha})$.*/
mydbl geodesic_equation_theta(mydbl pr, mydbl ptheta, mydbl pphi, mydbl r, mydbl theta)
{
  mydbl f = 0.5*sinl(2*theta)*pphi*pphi - 2*(1/r - der_potential(r)/powl(C,2))*pr*ptheta;
  return f;
}

/*Function of the 3th (azimuthal) momentum component differential equation for the geodesics.
${p3}^{dot} = f3(x^{\alpha},p^{\alpha})$.*/
mydbl geodesic_equation_phi(mydbl pr, mydbl ptheta, mydbl pphi, mydbl r, mydbl theta)
{
  mydbl f = 2*(1/r - der_potential(r)/powl(C,2))*pr*pphi - 2*(1/tanl(theta))*ptheta*pphi;
  return f;
}

/*Function for solving the geodesics differential equations using Euler's method.
Arguments are pointer so variables in that memory addresses are changed every time this function is called.*/
void euler1(mydbl *x0, mydbl *r, mydbl *theta, mydbl *phi, mydbl *p0, mydbl *pr, mydbl *ptheta, mydbl *pphi, mydbl *lambda)
{
  /*Increment in the variables of the differential equation we want to solve*/
  mydbl dx0, dr, dtheta, dphi, dp0, dpr, dptheta, dpphi;

  /*Calculation of the increments*/
  dx0 = *p0*DLAMBDA; dr = *pr*DLAMBDA; dtheta = *ptheta*DLAMBDA; dphi = *phi*DLAMBDA;
  dp0 = geodesic_equation_0(*p0, *pr, *r)*DLAMBDA;
  dpr = geodesic_equation_r(*p0, *pr, *ptheta, *pphi, *r, *theta)*DLAMBDA;
  dptheta = geodesic_equation_theta(*pr, *ptheta, *pphi, *r, *theta)*DLAMBDA;
  dpphi = geodesic_equation_phi(*pr, *ptheta, *pphi, *r, *theta)*DLAMBDA;

  /*New values of the variables of the differential equation. Since we are using pointers, when called the routine the value of variable change.*/
  *x0 = *x0 + dx0; *r = *r + dr; *theta = *theta + dtheta; *phi = *phi + dphi;
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
  mydbl x0 = 0.0, r = 500, theta = 0.0, phi = 0.0, p0 = 0.001, pr = -p0, ptheta = 0.0, pphi = 0.0, lambda = 0.0;

  /*Pointer to file where solution of differential equation will be saved.*/
  FILE *geodesic;
  geodesic = fopen("geodesic_solution.dat","w");

  /*Write line of initial values in file*/
  fprintf(geodesic, "%.18Lf %.18Lf %.18Lf %.18Lf %.18Lf %.18Lf %.18Lf %.18Lf %.18Lf\n", lambda, x0, r, theta, phi, p0, pr, ptheta, pphi);

  /*Solution of the differential equation*/
  for(i=0; i<(1+ NLINES); i++)
    {
      euler1(&x0, &r, &theta, &phi, &p0, &pr, &ptheta, &pphi, &lambda);
      fprintf(geodesic, "%.18Lf %.18Lf %.18Lf %.18Lf %.18Lf %.18Lf %.18Lf %.18Lf %.18Lf\n", lambda, x0, r, theta, phi, p0, pr, ptheta, pphi);
    }

  /** Releasing all used space in memory **/
  fclose(geodesic); //Close file storing the results
}
