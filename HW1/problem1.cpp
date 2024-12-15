#include <iostream>
#include <vector>

////////////////////////////////////////////////////////////////////////////

#include "Simulation/Simulation.hpp"
#include "Energy/Leanard_Jones.hpp"

////////////////////////////

int main(int argc, char** argv) {

  if (argc != 2) {
    std::cout << "You must provide an input file!" << std::endl;
    std::cout << "Usage: $ hw1 <input xyz file>" << std::endl;
    exit(1);
  }

  // Read input xyz file
  Simulation sim(argv[1]);
  std::cout << "Input:" << std::endl;
  std::cout << sim << std::endl;

  std::vector<int> bad_atoms;

  for (int i=0; i<sim.n_atoms; i++) {

    if (sim.elements[i] != 79) {
      bad_atoms.push_back(i);
    }
  }

  if (bad_atoms.size() > 0) {

      std::cout << "ERROR: Only Gold atoms accepted..." << std::endl;
      std::cout << "Atoms ";
      for (int atom_idx : bad_atoms) {
        std::cout << atom_idx << " ";
      }
      std::cout << "are bad!" << std::endl;
      exit(1);
  }

  // Compute and report energy
  std::cout << "Energy: " << LJ_energy(sim) << std::endl;

  return 0;
}
