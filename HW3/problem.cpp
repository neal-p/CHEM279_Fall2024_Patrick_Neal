#include <iostream>
#include "Simulation/Simulation.hpp"

#include "../external/eigen-3.4.0/Eigen/Dense"
#include "../external/eigen-3.4.0/Eigen/Eigenvalues"
using namespace Eigen;


// Make small values 0
// helpful for debugging these bigger matrices
void zero_out_small(MatrixXd& inout)
{
    inout = (1e-10 < inout.array().abs()).select(inout, 0.0f);
}
void zero_out_small(VectorXd& inout)
{
    inout = (1e-10 < inout.array().abs()).select(inout, 0.0f);
}




int main(int argc, char** argv) {

  if (argc != 3) {
    std::cout << "You must provide an input files!" << std::endl;
    std::cout << "Usage: $ problem1 <input xyz file> <input basis file>" << std::endl;
    exit(1);
  }

  Simulation sim(argv[1]);
  sim.load_basis(argv[2]);

  std::cout << "Read in the following:" << std::endl;
  std::cout << sim << std::endl;

  if (sim.n_e % 2 != 0) {
      std::cout << "Simulation must have EVEN number of electrons!! Exiting..." << std::endl;
      exit(1);
  }


  MatrixXd S(sim.AOs.size(), sim.AOs.size());

  for (int i=0; i < sim.AOs.size(); i++) {
      for (int j=0; j < sim.AOs.size(); j++) {

          S(i,j) = AO_overlap(sim.AOs[i], sim.AOs[j]);

      }
  }
  zero_out_small(S);

  std::cout << "Overlap Matrix S:" << std::endl;
  std::cout << S << std::endl;
  std::cout << std::endl;


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

  zero_out_small(H);

  std::cout << "Huckel Hamiltonian Matrix H:" << std::endl;
  std::cout << H << std::endl;
  std::cout << std::endl;

    // Solve the inversion thingy to get X
  SelfAdjointEigenSolver<Eigen::MatrixXd> overlap_solver(S);
  MatrixXd U = overlap_solver.eigenvectors(); 
  VectorXd s = overlap_solver.eigenvalues(); 

  // Get X matrix
  //     Lecture 9 - slide 14
  //     set inverse square root of the eigenvalues VECTOR as the diagonal 
  //     so that we can multiply it by U and U^t
  MatrixXd X = U * s.cwiseSqrt().cwiseInverse().asDiagonal() * U.transpose();
  zero_out_small(X);

  std::cout << "Matrix X:" << std::endl;
  std::cout << X << std::endl;
  std::cout << std::endl;

  // Check orthogonalization 
  //    also Lecture 9 - slide 14
  //    lecture 10 - slide 19 good example of the diagonal bit
  MatrixXd should_be_orthogonal = X.transpose() * S * X;
  zero_out_small(should_be_orthogonal);
  bool is_identity = (bool)should_be_orthogonal.isIdentity(1e-8);

  std::cout << "Check orthogonalizatoin X^t * S * X:" << std::endl;
  std::cout << should_be_orthogonal << std::endl;

  std::cout << "    is identity check: " << is_identity << std::endl;
  if (!is_identity) {
      std::cout << "    HEY ORTHOGONALIZATION DIDNT WORK RIGHT, ENERGIES PROP NOT RIGHT EITHER " << std::endl;
  }
  std::cout << std::endl;

  // Now get *fancy* H ('H prime')
  //    lecture 10 - slide 119 example
  //    S^-1/2 = X from above
  MatrixXd H_p = X.transpose() * H * X;

  // Solve for H_p V = V E
  SelfAdjointEigenSolver<MatrixXd> hamiltonain_solver(H_p);
  MatrixXd V = hamiltonain_solver.eigenvectors();
  VectorXd E = hamiltonain_solver.eigenvalues();
  zero_out_small(V);
  zero_out_small(E);

  std::cout << "Matrix V: " << std::endl;
  std::cout << V << std::endl;
  std::cout << std::endl;

  // Un-orthogonalize V to get the original MO coefficients C
  MatrixXd C = X * V;
  zero_out_small(C);
  std::cout << "Matrix C: " << std::endl;
  std::cout << C << std::endl;
  std::cout << std::endl;

  // Orbital energies
  std::cout << "Orbital energies Vector E:" << std::endl;
  std::cout << E << std::endl;
  std::cout << std::endl;


  // Calculate total energy
  //    lecture 
  double total_energy = 0;

  for (int i=0; i < sim.n_v_e / 2; i++) {
      total_energy += 2 * E[i];
  }

  std::cout << "TOTAL ENERGY: " << std::endl;
  std::cout << total_energy << std::endl;

  return 0;
}
