# light_geodesics_thesis
This repository contains the calculations of light shift using the formula obtained by Merlin and Salgado.

The first stage consists in calculating the scale factor and the derivative of scale factor as a function of cosmic and conformal time. This is done in the program scale_factor.c which produces a data file with the desired output.

To run the program write in the terminal:
make scale_factor
./scale_factor.x

This produces a data file described above. Cosmic time corresponds to the first column, conformal time to the second column and the third column corresponds to the scale factor.

The second stage consists in calculate the geodesics equations. The first part is done in geodesic_equations.c, this contains the evaluation of the differential equations which are written in the for  $\frac{d(x or p)^{\alpha}}{d\lambda}=f(x^{\alpha},p^{\alpha})$.
Up to now I must be introducing the interpolation of scale factor and a simple Euler method to solve the differential equations.