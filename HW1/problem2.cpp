#include <iostream>
#include <vector>
#include <fstream>

////////////////////////////////////////////////////////////////////////////

#include "Simulation/Simulation.hpp"
#include "Energy/Leanard_Jones.hpp"
#include "Gradients/Finite_Differences.hpp"

////////////////////////////

int main(int argc, char** argv) {

  if (argc != 3) {
    std::cout << "You must provide an input and output file!" << std::endl;
    std::cout << "Usage: $ hw2 <input xyz file> <output csv file>" << std::endl;
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
  std::cout << "Initial Energy: " << LJ_energy(sim) << std::endl;

  // Gradients
  std::cout << "Analytical Gradient: " << std::endl;
  MatrixXd analytical = LJ_gradient_analytical(sim);
  std::cout << analytical << std::endl;
  std::cout << std::endl;


  std::ofstream csv;
  csv.open(argv[2]);

  // Check output file exists
  if (!csv.is_open()) {
    std::cout << "ERROR: could not open output file " << argv[2] << std::endl;
    exit(1);
  }

  // Header for csv file
  csv << "atom_idx,analytical_x,analytical_y,analytical_z,h,fd_x,fd_y,fd_z,cd_x,cd_y,cd_z" << std::endl;


  std::vector<double> hs = {0.1, 0.01, 0.001, 0.0001};

  for (double h : hs) {

    std::cout << "Forward Difference Gradient (h=" << h << "):" << std::endl;
    MatrixXd fd = forward_difference(sim, LJ_energy, h);
    std::cout << fd << std::endl;

    std::cout << "Central Difference Gradient (h=" << h << "):" << std::endl;
    MatrixXd cd = central_difference(sim, LJ_energy, h);
    std::cout << cd << std::endl;

    std::cout << std::endl;

    for (int i=0; i<sim.n_atoms; i++) {

      csv << i << ","

          << analytical(0,i) << ","
          << analytical(1,i) << ","
          << analytical(2,i) << ","

          << h << ","

          << fd(0,i) << ","
          << fd(1,i) << ","
          << fd(2,i) << ","

          << cd(0,i) << ","
          << cd(1,i) << ","
          << cd(2,i)

          << std::endl;
    }
  }



  return 0;
}
