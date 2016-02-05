# light_geodesics_thesis

## Description
This repository contains the calculations of light shift using the formula obtained by Merlin and Salgado.
The use of these programas asume you have GNU Scientific Libraries installed in your computer.

## Stages of the work
### Calculating scale factor
The first stage consists in calculating the scale factor and the derivative of scale factor as a function of cosmic and conformal time. This is done in the program scale_factor.c which produces a data file with the desired output.
By default the program is setup to solve scale factor for ${\Omega}_{matter}=0.3$ and ${\Omega}_{\Lambda}=0.7$, this can be change in the first part of the program where constants are defined.
To run the program write in the terminal:
```
make scale_factor
./scale_factor.x
```
This produces the data file described above. Cosmic time corresponds to the first column, conformal time to the second column, third column corresponds to the scale factor and the fourth column corresponds to the derivative of scale factor.

### Solving geodesics equations
The second stage consists in calculate the geodesics equations. This is done in geodesic_equations.c, which contains the evaluation of the differential equations for FRLW model with an scalar perturbation corresponding to a Plummer potential.
Since geodesic equations are second order differential equations for the coordinates $x^{\alpha}$ (in the code I denote $x^{0}=eta$), the system must be converted to a first order set of differential equations. This is easily done by setting $p^{\alpha} = {\dot{x}}^{\alpha}$.
Up to now interpolation of scale factor and the derivative of scale factor is done and the functions for the potential and the potential derivative have been made to work.
When writing the next commands on the terminal you can obtain some graphs of scale factor, derivative of it and the potential an its derivative:
```
make geodesic_equations
./geodesic_equations.x
```