#!/bin/bash

make clean
make geodesic_solution_spherical
sleep 3
./geodesic_solution_spherical.x
