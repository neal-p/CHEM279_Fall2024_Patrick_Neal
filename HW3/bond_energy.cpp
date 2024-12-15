#include <iostream>
#include "Simulation/Simulation.hpp"

#include "../external/eigen-3.4.0/Eigen/Dense"
#include "../external/eigen-3.4.0/Eigen/Eigenvalues"
using namespace Eigen;


double total_energy(Simulation& sim) {

  MatrixXd S(sim.AOs.size(), sim.AOs.size());

  for (int i=0; i < sim.AOs.size(); i++) {
      for (int j=0; j < sim.AOs.size(); j++) {

          S(i,j) = AO_overlap(sim.AOs[i], sim.AOs[j]);

      }
  }

  MatrixXd H(sim.AOs.size(), sim.AOs.size());

  int i=0;
  for (AO& ao : sim.AOs) {

      int element = sim.elements[ao.owning_atom_idx];
      int total_angular_momentum = ao.l + ao.m + ao.n;

      if (element == 1) {

          // H 1S orbital
          if (total_angular_momentum == 0) {
              H(i,i) = -13.6;
          } else {
              std::cout << "cant handle H with angular momentum > 0" << std::endl;
          }

      } else if (element == 6) {

          // C 2S orbital
          if (total_angular_momentum == 0) {
              H(i,i) = -21.4;
          } else if (total_angular_momentum == 1) {
              H(i,i) = -11.4;
          } else {
              std::cout << "cant handle C with angular momentum > 1" << std::endl;
          }
      }

      i++;
  }

  double K = 1.75;

  for (int i=0; i<sim.AOs.size(); i++) {
    for (int j=0; j<sim.AOs.size(); j++) {

        // Only off-diagonal
        if (i != j) {
            H(i,j) = K * 0.5 * (H(i,i) + H(j,j)) * S(i,j);
        }
    }
  }

    // Solve the inversion thingy to get X
  SelfAdjointEigenSolver<Eigen::MatrixXd> overlap_solver(S);
  MatrixXd U = overlap_solver.eigenvectors(); 
  VectorXd s = overlap_solver.eigenvalues(); 

  // Get X matrix
  //     Lecture 9 - slide 14
  //     set inverse square root of the eigenvalues VECTOR as the diagonal 
  //     so that we can multiply it by U and U^t
  MatrixXd X = U * s.cwiseSqrt().cwiseInverse().asDiagonal() * U.transpose();

  // Check orthogonalization 
  //    also Lecture 9 - slide 14
  //    lecture 10 - slide 19 good example of the diagonal bit
  MatrixXd should_be_orthogonal = X.transpose() * S * X;
  bool is_identity = (bool)should_be_orthogonal.isIdentity(1e-8);

  if (!is_identity) {
      std::cout << "    HEY ORTHOGONALIZATION DIDNT WORK RIGHT, ENERGIES PROP NOT RIGHT EITHER " << std::endl;
  }
  // Now get *fancy* H ('H prime')
  //    lecture 10 - slide 119 example
  //    S^-1/2 = X from above
  MatrixXd H_p = X.transpose() * H * X;

  // Solve for H_p V = V E
  SelfAdjointEigenSolver<MatrixXd> hamiltonain_solver(H_p);
  MatrixXd V = hamiltonain_solver.eigenvectors();
  VectorXd E = hamiltonain_solver.eigenvalues();

  // Un-orthogonalize V to get the original MO coefficients C
  MatrixXd C = X * V;

  // Calculate total energy
  //    lecture 
  double total_energy = 0;

  for (int i=0; i < sim.n_v_e / 2; i++) {
      total_energy += 2 * E[i];
  }

  return total_energy;
}



int main(int argc, char** argv) {

  if (argc != 3) {
    std::cout << "You must provide an input files!" << std::endl;
    std::cout << "Usage: $ bond_energy <h2 xyz> <input basis file>" << std::endl;
    exit(1);
  }

  Simulation h2(argv[1]);
  h2.load_basis(argv[2]);

  double h2_energy = total_energy(h2);
  double h_energy = -13.6;


  std::cout << "Computed energy of H2: " << h2_energy << " eV" << std::endl;
  std::cout << "Given energy of H: " << h_energy << " eV" << std::endl;

  std::cout << "H2 bond energy: " << h2_energy - (2*h_energy) << " eV" << std::endl;



  return 0;
}
