# Overview of PIMC

This program was created as a sandbox program to learn about path integral
simulations. Since this is a learning tool, many of the options are defined
at compile-time rather than appearing in the input. Most of these options 
should only be changed if you have a thorough understanding of the code.

Currently PIMC.exe can simulate atoms, molecules, electrons, and positrons.
The available potentials are Lennard-Jones, Coulomb, harmonic bond, harmonic
angle, and eFF (electrons/positrons, Phys. Rev. Lett.,99,185003,2007).

The program runs in parallel through OpenMP. Some loops are parallel over
the number of atoms, while others are over the number of beads. The number of
atoms and beads limits the number of CPUs that the program can use efficiently.
PIMC.exe is most efficient if both the number of atoms and the number of beads
are greater than or equal to the number of CPUs.

When PIMC.exe is executed, the program looks for two files (pimc.param and
pimc.pot). Descriptions of the files and commands are given below.

## Simulation Input Description (pimc.param)

Note: The format of this file is fixed by the code. The order of the options
cannot be changed and the spaces should be replaced with underscores. PIMC.exe
reads a single string before reading the value(s).

### Options

Ensemble: NVT or NPT simulations. If the NPT ensemble is chosen, the system
should be large enough to ensure that the simulation box is twice the size of
the cutoff of the Lennard-Jones and Coulomb potentials.

Temperature: The temperature of the system must be greater than zero. Since
quantum effects are more important at low temperatures, lowering the
temperature may require more beads.

Pressure: External pressure applied to the system. The pressure can be positive
or negative, but is only important in NPT simulations.

Number of Eq steps: Equilibration steps before the trajectory and energies are
recorded. Simulations containing a large number of atoms and beads will require
a large number of steps to reach equilibrium. If the system has a large number
of atoms, it may be beneficial to run a classical (Num. beads == 1) simulation
first and use the last trajectory frame as the starting structure of the 
quantum simulation (Num. beads > 1). Only accepted moves are counted.

Number of MC steps: Number of Monte Carlo steps for the production run. Only
accepted moves are counted.

Number of beads: Number of beads or "time slices" to include in the simulation.
The number must be greater or equal to 1, and 1 gives a classical simulation.

Acceptance ratio: Ratio of accepted moves to the total moves attempted. The
maximum size of the MC moves are adjusted during the equilibration to force
the acceptance ratio to be the specified value. A small acceptance ratio allows
the atoms and beads to move further, however, the simulations will be more 
expensive. The value should typically be between 0.2 and 0.5 (20-50% of the
moves are accepted).

Traj. steps before printing: Number of MC steps between trajectory frames.

Print mode: The trajectory can be printed for the center-of-mass positions of
the atoms (COM/com) or for all beads (All/ALL). 

XYZ file: File name for the starting coordinates. The number of atoms in the
simulation is determined from the XYZ file as well.

Box size: The size of the simulation box in the x, y, and z directions. The
sizes must be at least twice the cutoff.

Cutoff: Cutoff for the Lennard-Jones and Coulomb potentials. Force fields are
typically parametrized for cutoffs between 9-12 angstroms. The force field
will not give accurate results if the cutoffs are changed, however, PIMC.exe
only allows for a single cutoff.

Connectivity file: A file describing the harmonic bond and harmonic angle
potentials. If there are harmonic bonds, the file can be listed as None/none.
An example of the format is show below.
Bonds: Num. bonds
<atom num. 1> <atom num. 2> <force constant (eV/A^2)> <Eq. length (A)>
...
Angles: Num. angles
<cent. atom> <atom 1> <atom 2> <force constant (eV/rad^2)> <Eq. angle (deg.)>
...

### Options for free electrons or positrons using the eFF model

Note: The eFF model represents electrons as spherical Gaussian charge
distributions which have 5 degrees of freedom (x, y, z, radius, and spin). This
model was designed to produce semi-classical simulations of electron dynamics,
and hence, does not have the correct behavior in PIMC simulations. Within
PIMC.exe positrons can be simulated as positively charged electrons.

Spin mode: Since each eFF particle has a spin, it greatly speeds up the 
equilibration if particles are allowed to change their spin. The total spin 
can be restricted to a low spin configuration (LPIMC/lpimc/low/LOW) or the 
total spin can change (high/HIGH/HPIMC/hpimc).

Number of pre-Eq. steps: Since the electrons are placed randomly, the 
initial configurations often high energy states. Setting the number of
pre-equilibration steps to be greater than zero allows the electrons to relax
before the nuclei and atoms are allowed to move.

Electrons: Number of free eFF electrons.

Positrons: Number of free eFF positrons. Note that electron-positron
interactions are purely Coulombic.

## Non-Bonded Potential Input (pimc.pot)

Note: This file contains information for the Lennard-Jones and Coulomb
potentials, as well as the masses of the atoms. The properties are assigned
to the atoms using the atom types, therefore extra types can be defined in the
file without causing errors. Make sure to check the atom types in the XYZ files
since an atom type that is not in this file will have no mass and no 
interactions.

Number of types: Number of atom types listed in the file. The number of entries
in the remaining sections must match this number.

Types: A list of names for the atom types. All remaining sections much have the
same order.

Masses: Masses of the atoms in amu.

Charges: Atomic charges in units of the electron charge.

Lennard-Jones: Epsilon (eV) and sigma (A) for the atom types. The potential is
given by
E(r) = 4*<epsilon>*((<sigma>/r)^12-(<sigma>/r)^6)
and the parameters for different atom types are determined through combination
rules.
