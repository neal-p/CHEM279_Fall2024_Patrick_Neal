#include "Gradients/Finite_Differences.hpp"

MatrixXd central_difference(Simulation& sim, EnergyFunc E, double h) {

  MatrixXd gradients(3, sim.n_atoms);

  // For each atom, for each dimension x,y,z
  for (int atom_idx=0; atom_idx < sim.n_atoms; atom_idx++) {
    for (int d=0; d < 3; d++) {

      // perterb atom
      // eval energy difference
      sim.xyz(atom_idx, d) -= h;
      double E0 = E(sim); 

      sim.xyz(atom_idx, d) += (2*h);
      double E1 = E(sim); 

      gradients(d, atom_idx) = E1 - E0;
      sim.xyz(atom_idx, d) -= h;
    }
  }

  return gradients.array() / (-2 * h);
}


MatrixXd forward_difference(Simulation& sim, EnergyFunc E, double h) {

  double E0 = E(sim); 
  MatrixXd gradients(3, sim.n_atoms);
  gradients.fill(-1 * E0);

  // For each atom, for each dimension x,y,z
  for (int atom_idx=0; atom_idx < sim.n_atoms; atom_idx++) {
    for (int d=0; d < 3; d++) {

      // perterb atom
      // eval energy difference
      sim.xyz(atom_idx, d) += h;
      double E1 = E(sim); 

      gradients(d, atom_idx) += E1;
      sim.xyz(atom_idx, d) -= h;
    }
  }

  return gradients.array() / (-h);
}


