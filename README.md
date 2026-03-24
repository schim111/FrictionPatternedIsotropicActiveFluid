# FrictionPatternedIsotropicActiveFluid
MATLAB code that solves dynamical equations for an isotropic active fluid in the presence of anisotropic friction.

The script ActiveFrictionStartUp.m is a startup script for the solver and sets model parameters, mesh geometry, friction pattern, boundary conditions, and initial conditions. Once these have been set, call this script to run the solver.
BJGSolver2D.m is the actual solver given the parameters and initial condition from the startup file. This file does not need to be edited.

The solver uses the MATLAB/C++ finite element package FELICITY: (https://github.com/walkersw/felicity-finite-element-toolbox) which must be downloaded and installed. The scripts MatAssemCAniso2D.m and MatAssemVAniso2D.m must be compiled using FELICITY before running the solver. They do not need to be edited.
