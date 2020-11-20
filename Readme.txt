A super-vague-and-not-informative-enough readme file:

Will Riedel wrote this code (pretty badly and haphazardly) in 2017, with some updates later in 2018. It wasn’t well planned out, and has plenty of gaps and poor form. This readme file is also very incomplete, but feel free to follow up at wriedel@stanford.edu.

Gas can be initialized as homogeneously distributed in a rectangular simulation volume, or split between two regions A and B with different conditions (used for benchmarking), or with no gas at time 0 and an input source with a finite pulse length. If an input source is used, plasma passes through a surface into the gun volume from a reservoir assumed to be at a constant pressure (flux density is n*C/4).

Ultimately the code was used to estimate mass flow into a cylindrical gun volume. The valve opens from the left, and gas travels into the gun volume. An image of the gun volume is shown in simplified_wall_geometry_2018-9-13.pdf. 

Because the code is limited in various ways (e.g. only single-core processing), it is too slow to handle atmospheric densities over large areas, and the estimates were made assuming free stream collisionless flow. However, the collisional behavior has been validated for homogenous density and for flux between split regions at various conditions, so collisions can be included at lower densities.

Input files
Input_Parameters:   Many general input parameters, see Read_Input_Data.f90 for slightly more description.
x_inlet.txt:      If using an exponentially varying grid (I probably wouldn’t), sets where the exponential variation begins.
y_inlet.txt:      The upper and lower edges of the source reservoir
num_walls.txt:    The number of walls included in the simulation
x_walls.txt:      Set endpoints of the wall boundaries that are included in the simulation. Each row are endpoints of a flat wall (x1,y1,x2,y2)

Main.f90 has the main loop. The general procedure at each time step is:
1. Particles are advanced in time through collisionless motion.
2. If applicable, particle velocities are adjusted based on a collisional sampling procedure.
3. Boundary conditions are applied, including particle reflections at walls and inflow/outflow of particles at sources and exit planes.

Position, velocity, and number of particles in each cell for each saved time step is saved in plain text in a directory defined in Input_Parameters. Also saved are vectors of various values (number particles added to simulation, number of particles in simulation, candidate/accepted collisional pairs) at each time step.

Anything referring to SPLIT: I modified this code to run benchmarking simulations of two regions at different conditions separated by a slit opening (sides A and B). If I were doing this over again from scratch I would be much more general about how initial conditions could be supplied.