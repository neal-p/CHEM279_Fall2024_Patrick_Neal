#include "Energy/Leanard_Jones.hpp"

double LJ_energy(Simulation& sim) {

  int pairwise_indices = sim.n_atoms * (sim.n_atoms - 1) / 2;

  // compute pairwise distances
  ArrayXd e_s(pairwise_indices);
  ArrayXd o_s(pairwise_indices);
  ArrayXd distances(pairwise_indices);

  int _idx = 0;
  for (int i=0; i < sim.n_atoms; i++) {
    for (int j=i+1; j < sim.n_atoms; j++) {

      double distance = (sim.xyz.row(i) - sim.xyz.row(j)).norm();

      // Currently always the same for gold atoms
      // but may need to lookup based on element
      // in the future
      double e = e_au;
      double o = o_au;

      e_s(_idx) = e;
      o_s(_idx) = o;
      distances(_idx) = distance;

      _idx++;
    }
  }

  // Calculate LJ potential terms
  ArrayXd div = o_s / distances;
  ArrayXd t6  = 2 * div.pow(6);
  ArrayXd t12 = div.pow(12);

  return (e_s * (t12 - t6)).sum();
}


MatrixXd LJ_gradient_analytical(Simulation& sim) {

  MatrixXd gradients(3, sim.n_atoms);

  for (int i=0; i < sim.n_atoms; i++) {
    for (int j=i+1; j < sim.n_atoms; j++) {

      // distance rij
      double rij = (sim.xyz.row(i) - sim.xyz.row(j)).norm();

      double t12 = 12 * pow(o_au, 12) * pow((1.0 / rij), 13);
      double t6 =  12 * pow(o_au, 6)  * pow((1.0 / rij), 7);

      // taking de with respect to each coordinate -> operate on vector between atoms
      MatrixXd de = e_au * (t12 - t6) * (1.0 / rij) * (sim.xyz.row(i).array() - sim.xyz.row(j).array());

      // Force is equal and opposite on each atom
      gradients.col(i).array() += de.transpose().array();
      gradients.col(j).array() -= de.transpose().array();

    }
  }

  return gradients;
}

