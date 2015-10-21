/*This program is the first step in solving the Merlin-Salgado equation.
It produces a .dat file with the scale factor in terms of the conformal time and the cosmic time.*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

#define NLINES 11

/*Structure with parameters needed for the function to integrate
hubble: The Hubble parameter
omega_m: Density parameter for matter (dark matter and baryonic)
omega_l: Density parameter for dark energy*/
struct Iparams
{
  double hubble;
  double omega_m;
  double omega_l;
};

/*This function multiplied by the differential of the scale factor is the differential of cosmic time*/
double integrando_cosmictime(double a, void *parameters)
{
  struct Iparams *p = (struct Iparams *) parameters; //Structure of parameters
  
  double h0 = p->hubble;
  double om = p->omega_m;
  double ol = p->omega_l;
  double f = (1/h0)/(sqrt(om*pow(a,-1) + ol*pow(a,2))); //Integrate to cosmic time
  
  return f;
}

/*This function multiplied by the differential of the scale factor is the differential of conformal time*/
double integrando_conformaltime(double a, void *parameters)
{
  struct Iparams *p = (struct Iparams *) parameters; //Structure of parameters
  
  double h0 = p->hubble;
  double om = p->omega_m;
  double ol = p->omega_l;
  double f = (1/h0)/(sqrt(om*a + ol*pow(a,4))); //Integrate to conformal time
  
  return f;
}

int main (void)
{
  struct Iparams parameters; // Structure of parameters for functions to integrate
  double om = 0.3, ol = 1.0-om; //Density parameters
  double h0 = 100; //Hubble parameter
  /*Asign values to the structure of parameters*/
  parameters.hubble = h0;
  parameters.omega_m = om;
  parameters.omega_l = ol;

  /*File to be written with data*/
  FILE *pf;
  pf = fopen("frw.dat","w");

  /*GSL object. Necessary object to separate enough space in memory for integration. 1 and 2 refer to cosmic and conformal time respectively*/
  gsl_integration_workspace *w1;
  gsl_integration_workspace *w2;
  
  /*result: Result of integration. error: Maximum relative error. 1 and 2 refer to cosmic and conformal time respectively*/
  double result1, result2, error1, error2, epsrel = 1.0e-11;
  size_t limit = 1000; //Maximum number of subintervals for integration
  /*a: Inferior integration limit. b: Superior integration limit*/
  double a = 0.0, b = 0.0;

  /*GSL object. Contain all information about the integrand.*/
  gsl_function Fcosmic, Fconformal;
  /*Defining Fcosmic*/
  Fcosmic.function = &integrando_cosmictime;
  Fcosmic.params = &parameters;
  /*Defining Fconformal*/
  Fconformal.function = &integrando_conformaltime;
  Fconformal.params = &parameters;

  int j=1; //Index for arrays
  /*Arrays for cosmic time, conformal time and scale factor*/
  double cos_tim[NLINES], con_tim[NLINES], sca_fac[NLINES];
  
  /*Initial values*/
  cos_tim[0] = 0.0; con_tim[0]=0.0; sca_fac[0] = 0.0;

  /*Allocate space in memory*/
  w1 = gsl_integration_workspace_alloc(limit);
  w2 = gsl_integration_workspace_alloc(limit);

  for(b = 0.001;b <= 1.0;b += 0.001)
    {
      /*Numerical integration of functions*/
      gsl_integration_qags(&Fcosmic, a, b, 0, epsrel, limit, w1, &result1, &error1);
      gsl_integration_qags(&Fconformal, a, b, 0, epsrel, limit, w2, &result2, &error2);
      fprintf(pf, "%.12f %.12f %.12f\n", result1, result2, b);
      /*Save values in an array*/
      //cos_tim[j] = result1; con_tim[j] = result2; sca_fac[j] = b;
      //j = j+1;
    }

  /*Close streams and free space in memory*/
  fclose(pf);
  gsl_integration_workspace_free(w1);
  gsl_integration_workspace_free(w2);
}
