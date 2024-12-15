#include <iostream>
#include <vector>
#include <fstream>

////////////////////////////////////////////////////////////////////////////

#include "Simulation/Simulation.hpp"
#include "Energy/Leanard_Jones.hpp"
#include "Gradients/Finite_Differences.hpp"
#include "Optimization/Steepest_Descent.hpp"

////////////////////////////

int main(int argc, char** argv) {

  if (argc != 2) {
    std::cout << "You must provide an input!" << std::endl;
    std::cout << "Usage: $ hw3 <input xyz file>" << std::endl;
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

  // Optimize
  GradFunc grad = [](Simulation& sim) -> MatrixXd {return central_difference(sim, LJ_energy, 0.0001);};
  Simulation opt = steepest_descent_line_search(sim, LJ_energy, grad, 500, 0.01);

  return 0;
}
