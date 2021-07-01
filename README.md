# Grease Rheology Modelling
This repository contains datafiles and inputfiles that can be used with the molecular dynamics package LAMMPS (https://lammps.sandia.gov/).

### REQUIREMENTS
To run the simulations in this repository you need to have LAMMPS installed. Different ways to install LAMMPS are described in the LAMMPS documentation: https://lammps.sandia.gov/doc/Install.html

### DATAFILES
Most datafiles include systems of randomly oriented chains of particles. The folder "OrthogonalBox" contains systems with an orthogonal simulation box, which are useful for equilibrium simulations. The folder "TriclinicBox" contains simulation boxes with a small amount of tilt, which are useful for simulations in which the simulation box is deformed (Alternatively, you can use an orthogonal box combined with the change box command in your input file for this). A folder "Ordered" was later added. This folder contains systems with chains that are aligned along the z-axis of the simulation box. The simulation boxes in this folder are all orthogonal.

All folders will include the subfolders "ChainLength" and "VolumeDensity". All the datafiles in the "ChainLength" folder have a volume density of chains equal to 0.15 (15%) while chain lengths are varied. The length of the chains (in σ) are indicated in the filename, e.g. data.CL*20*_1.gz corresponds to a chain length of 20σ. Similarly, the folder "VolumeDensity" contains datafiles with chains of length 20σ while the volume density of chains is varied. This is also indicated in the filename in a similar way, e.g. data.D*20*_1.gz corresponds to a volume density of 0.20 (20%).

### INPUTFILES
This repository contains three types of input file. Firstly, the in.GK file which uses the Green-Kubo method to calculate the viscosity of the system. Secondly, the in.DC file which calculates the diffusion coefficient based on the mean-squared displacement of the fibres. And finally, the in.Wiggle file which uses LAMMPS' wiggle command to subject the simulation box to oscillatory shear.
