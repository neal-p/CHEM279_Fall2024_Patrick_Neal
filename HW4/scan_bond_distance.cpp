#include "Simulation/Simulation.hpp"
#include <unordered_map>

#include "CNDO_2/CNDO_2.hpp"


double get_CNDO_2_total_energy(Simulation& sim) {

  PeriodicTable pt;

  // Based on nuclear coordinates,
  // these dont change for each scf cycle
  MatrixXd overlap = sim.Overlap();
  MatrixXd gamma = CalculateGamma(sim);

  MatrixXd Hcore = CalculateHcore(sim, gamma, overlap);

  // SCF 
  MatrixXd p_a(sim.AOs.size(), sim.AOs.size());
  MatrixXd p_b(sim.AOs.size(), sim.AOs.size());

  MatrixXd p_a_old(sim.AOs.size(), sim.AOs.size());
  MatrixXd p_b_old(sim.AOs.size(), sim.AOs.size());

  // Set initial guesss to 0 density
  p_a.fill(0.0);
  p_b.fill(0.0);

  // Initial Fock matrices
  MatrixXd Fock_a = CalculateFock(sim, gamma, p_a, p_b, overlap);
  MatrixXd Fock_b = CalculateFock(sim, gamma, p_b, p_a, overlap);


  // Iterate
  double density_change = 1.0 + 1e-6; // ensure at least one cycle

  int SCF_cycle = 0; 
  double energy;
  while (density_change > 1e-6) {

    SCF_cycle++;

   // Update Fock matrices
    overlap = sim.Overlap();
    gamma = CalculateGamma(sim);
    Fock_a = CalculateFock(sim, gamma, p_a, p_b, overlap);
    Fock_b = CalculateFock(sim, gamma, p_b, p_a, overlap);

    // Eq. 2.3/4
    SelfAdjointEigenSolver<MatrixXd> fock_a_solver(Fock_a);
    MatrixXd c_a = fock_a_solver.eigenvectors();
    VectorXd e_a = fock_a_solver.eigenvalues();
    

    SelfAdjointEigenSolver<MatrixXd> fock_b_solver(Fock_b);
    MatrixXd c_b = fock_b_solver.eigenvectors();
    VectorXd e_b = fock_b_solver.eigenvalues();
 
    // Copy old densities
    p_a_old = p_a;
    p_b_old = p_b;

    // Calculate updated densities
    p_a.fill(0.0);
    p_b.fill(0.0);

    int P = 4;
    int Q = 4;

    for (int i=0; i < sim.AOs.size(); i++) {
      for (int j=0; j < sim.AOs.size(); j++) {

        for (int p=0; p < sim.p; p++) {
          p_a(i,j) += c_a(i,p) * c_a(j,p);
        }

        for (int q=0; q < sim.q; q++) {
          p_b(i,j) += c_b(i,q) * c_b(j,q);
        }

      }
    }

    MatrixXd density_change_a = (p_a_old - p_a).array().abs();
    MatrixXd density_change_b = (p_b_old - p_b).array().abs();

    density_change = density_change_a.norm() + density_change_b.norm();

    // double total_energy = 

    overlap = sim.Overlap();
    gamma = CalculateGamma(sim);
    Hcore = CalculateHcore(sim, gamma, overlap);
    Fock_a = CalculateFock(sim, gamma, p_a, p_b, overlap);
    Fock_b = CalculateFock(sim, gamma, p_b, p_a, overlap);


    // Calculate total energy
    double first_sum_term = 0.0;
    double second_sum_term = 0.0;



    for (int u=0; u < sim.AOs.size(); u++) {
      for (int v=0; v < sim.AOs.size(); v++) {

        first_sum_term += p_a(u, v) * (Hcore(u,v) + Fock_a(u, v));
        second_sum_term += p_b(u, v) * (Hcore(u,v) + Fock_b(u, v));

      }
    }

    first_sum_term *= 0.5;
    second_sum_term *= 0.5;

    double third_sum_term = 0.0;
    for (int A=0; A < sim.n_atoms; A++) {

      for (int B=A + 1; B < sim.n_atoms; B++) {

        if (A != B) {
          Element A_element = pt.get_element(sim.elements[A]);
          Element B_element = pt.get_element(sim.elements[B]);

          VectorXd A_center = sim.xyz.row(A);
          VectorXd B_center = sim.xyz.row(B);

          double distance = (A_center - B_center).norm();

          double A_Z = A_element.z - (A_element.n_e - A_element.n_v_e);
          double B_Z = B_element.z - (B_element.n_e - B_element.n_v_e);

          third_sum_term += (A_Z * B_Z) / distance;

        }
      }
    }

  energy = first_sum_term + second_sum_term + third_sum_term;

  if (SCF_cycle > 100) {
    break;
  }
  }

  return energy;
}


int main(int argc, char** argv) {

  if (argc < 1) {
    std::cout << "You must provide an input files!" << std::endl;
    std::cout << "Usage: $ bond_distance_scan <input basis file> <geometry 1> <geometry 2> ... " << std::endl;
    exit(1);
  }

  std::cout << "file,energy" << std::endl;

  for (int i=2; i < argc; i++) {

    Simulation sim(argv[i]);
    sim.load_basis(argv[1]);

    double energy = get_CNDO_2_total_energy(sim);

    std::cout << argv[i] << "," << energy << std::endl;
  }

  return 0;
}


