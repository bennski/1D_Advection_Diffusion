"README.txt"

File: 1D_Advection_Diffusion.py
Author: Ben K. D. Pearce
Created: July 14, 2016
Updated: November 11, 2018

***Prerequisite packages***
1. python 2.7
2. ffmpeg

Description:
-The following code was written to simulate the mixing of a local concentration of molecules in a pond.
-The advection cell is driven by a temperature difference from the top to the bottom of the pond.
-The code produces an animation of the mixing in the moving frame of the convective fluid.

Boundary conditions:
-Cyclic (periodic) boundaries, representing a convection cell in a pond.

Initial condition:
-Sharp Guassian function at the base of the pond (x=0)

Note:
-The upwind method was chosen due to its lack of spurious oscillations, lowest error accumulation in mass conservation tests, 
and convergence for increasing levels of refinement. This method has mass losses of up to ~1% for 9 levels of refinement. 
This error decreases for increasing levels of refinement.


