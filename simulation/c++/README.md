# Description
This is the reference code for HW5. It reads the molecular input composed exclusively of H, C, N, O, and F atoms, along with the corresponding basis set, and calculates the numerical gradient of CNDO/2 energy with respect to the atomic coordinates.

To compile the code on Datahub, enter the command `make all` under this directory.

# Usage
```
./HW5 <filename>
```
`<filename>`: The path to the input file containing information of the molecule.

# MD Simulation Compilation

g++ -o md_simulation_cndo md_simulation_cndo.cpp AO.cpp CNDO.cpp util.cpp -larmadillo -std=c++11
./md_simulation_cndo


