#include "Simulation/Simulation.hpp"
#include <unordered_map>

#include "CNDO_2/CNDO_2.hpp"



int main(int argc, char** argv) {

  if (argc != 3) {
    std::cout << "You must provide an input files!" << std::endl;
    std::cout << "Usage: $ problem1 <input xyz file> <input basis file>" << std::endl;
    exit(1);
  }

  Simulation sim(argv[1]);
  sim.load_basis(argv[2]);


  PeriodicTable pt;

  // Based on nuclear coordinates,
  // these dont change for each scf cycle
  MatrixXd overlap = sim.Overlap();
  MatrixXd gamma = CalculateGamma(sim);

  std::cout << "p=" << sim.p << ", q=" << sim.q << std::endl;
  std::cout << std::endl;
  std::cout << "gamma:" << std::endl;
  std::cout << gamma * eV << std::endl;
  std::cout << std::endl;

  std::cout << "overlap:" << std::endl;
  std::cout << overlap << std::endl;
  std::cout << std::endl;

  MatrixXd Hcore = CalculateHcore(sim, gamma, overlap);

  std::cout << "Hcore:" << std::endl;
  std::cout << Hcore * eV << std::endl;
  std::cout << std::endl;


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
  while (density_change > 1e-6) {

    SCF_cycle++;

    std::cout << "----------------" << std::endl;
    std::cout << "SCF Cycle: " << SCF_cycle << std::endl;


    std::cout << "Starting Densities: " << std::endl;
    std::cout << p_a << std::endl;
    std::cout << std::endl;
    std::cout << p_b << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;

    // Update Fock matrices
    overlap = sim.Overlap();
    gamma = CalculateGamma(sim);
    Fock_a = CalculateFock(sim, gamma, p_a, p_b, overlap);
    Fock_b = CalculateFock(sim, gamma, p_b, p_a, overlap);

    std::cout << "Fock Matrices: " << std::endl;
    std::cout << Fock_a * eV << std::endl;
    std::cout << std::endl;
    std::cout << Fock_b * eV << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;

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


    std::cout << "End of SCF Cycle: " << SCF_cycle << std::endl;
    std::cout << "Change in density: " << density_change << std::endl;
    std::cout << std::endl;


    std::cout << "Ca:" << std::endl;
    std::cout << c_a << std::endl;
    std::cout << std::endl;

    std::cout << "Cb:" << std::endl;
    std::cout << c_b << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;


    std::cout << "Ea:" << std::endl;
    std::cout << e_a * eV << std::endl;
    std::cout << std::endl;

    std::cout << "Eb:" << std::endl;
    std::cout << e_b * eV << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;

    std::cout << "Ending Densities:" << std::endl;
    std::cout << p_a << std::endl;
    std::cout << std::endl;
    std::cout << p_b << std::endl;
    std::cout << std::endl;
    std::cout << "----------------" << std::endl;


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

  double energy = first_sum_term + second_sum_term + third_sum_term;

  std::cout << "Nuclear Repulsion Energy is " << third_sum_term  * eV << std::endl;
  std::cout << "Electron Energy is " << (first_sum_term + second_sum_term) * eV << std::endl;

  std::cout << "Total Energy=" << energy * eV << std::endl;

    if (SCF_cycle > 100) {
      break;
    }
  }



  return 0;
}


