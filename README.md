# CHEM279

This repo will contain all code for CHEM279 Fall 2024. 
A previous student suggested making different code bits re-usable across homeworks, so I tried to create a system where different functionality would be contained in submodules. For example, code to parse xyz input files is separated from code to calculate energies. 

## Homeworks
Each homework is given a folder at the top level. In the HWX folder the executables for each specific homeork will be written, compiled, and run. These executables will depend on code in the `src/` directory for different components. 

**Each HW folder will have its own README with compilation and run instructions!!**

Below I give an overview of what comes from what homework:
### HW1
 - `src/Simulation/Simulation.cpp`: basic code to read xyz file into a `Simulation` object to store elements and coordinates
 - `src/Energy/Leanard_Jones.cpp`: Defines the Leanard_Jones energy calculation and analytical force
 - `src/Gradients/Finite_Difference.cpp`: Methods for calculating numerical forces with finite difference (forward and central differences)
 - `src/Optimization/Steepest_Descent.cpp`: Local optimization methods

 - `HW1/problem1.cpp`: program to report energy
 - `HW1/problem2.cpp`: program to show truncation error in finite difference
 - `HW1/problem3.cpp`: program to optimize a cluster of atoms with steepest descent
   
### HW2
 - `src/Integrands/Gaussian.cpp`: Defines `Gauss` and `Shell` objects, as well as `ShellOverlap` function for problem2. 
 - `src/Integrals/Numerical.cpp`: Defines the `QSIMP` numerical integration method.

 - `HW2/problem1.cpp`: program to calculate 1D `Gauss` overlap integral.
 - `HW2/problem2.cpp`: program to calculate 3D `Shell` overlap integral.

### HW3
This homework required a big re-write to better handle the gaussians as actual AOs
 - `src/Gaussian/Gaussian.cpp`: Defines `AO` objects, `read_basis_set` function to parse a basis set of contracted gaussians, and `AO_overlap` function to compute their overlaps
 - `src/Simulation/Simulation.cpp`: was heavily modified and now has ability to store AOs and get more important information about the system, such as number of electrons and valence electrons

 - `HW3/problem.cpp`: program to read xyz file, read basis set file, perform Extended Huckel calculation
 
 ### HW4
 This homework I kept most of the changes from HW3 and added the CNDO code on top
 - `src/CONDO_2/CNDO_2.cpp`: this is where all of the CNDO related functions are
 - `src/Simulation/Simulation.cpp`: a few small tweaks to make the element information and setting P and Q electrons.
 - `HW4/problem.cpp`: Computes the CNDO_2 energy of the input molecule given a basis set file 
 - `HW4/scan_bond_distance.cpp`: Uses CNDO_2 to scan the bond distance of N2.
 
 ### HW5
 This homework builds on HW4
 - `HW5/problem.cpp`: Computes the CNDO_2 energy and gradient 
