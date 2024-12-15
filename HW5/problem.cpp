#include "Simulation/Simulation.hpp"
#include <unordered_map>
#include "Gaussian/Gaussian.hpp"

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

  //std::cout << "p=" << sim.p << ", q=" << sim.q << std::endl;
  //std::cout << std::endl;
  //std::cout << "gamma:" << std::endl;
  //std::cout << gamma * eV << std::endl;
  //std::cout << std::endl;

  //std::cout << "overlap:" << std::endl;
  //std::cout << overlap << std::endl;
  //std::cout << std::endl;

  MatrixXd Hcore = CalculateHcore(sim, gamma, overlap);

  //std::cout << "Hcore:" << std::endl;
  //std::cout << Hcore * eV << std::endl;
  //std::cout << std::endl;


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

  double energy, nuclear_repulsion, electron_energy;

  int SCF_cycle = 0; 
  while (density_change > 1e-6) {

    SCF_cycle++;

    //std::cout << "----------------" << std::endl;
    //std::cout << "SCF Cycle: " << SCF_cycle << std::endl;


    // std::cout << "Starting Densities: " << std::endl;
    // std::cout << p_a << std::endl;
    // std::cout << std::endl;
    // std::cout << p_b << std::endl;
    // std::cout << std::endl;
    // std::cout << std::endl;

    // Update Fock matrices
    overlap = sim.Overlap();
    gamma = CalculateGamma(sim);
    Fock_a = CalculateFock(sim, gamma, p_a, p_b, overlap);
    Fock_b = CalculateFock(sim, gamma, p_b, p_a, overlap);

    // std::cout << "Fock Matrices: " << std::endl;
    // std::cout << Fock_a * eV << std::endl;
    // std::cout << std::endl;
    // std::cout << Fock_b * eV << std::endl;
    // std::cout << std::endl;
    // std::cout << std::endl;

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


    // std::cout << "End of SCF Cycle: " << SCF_cycle << std::endl;
    // std::cout << "Change in density: " << density_change << std::endl;
    // std::cout << std::endl;


    // std::cout << "Ca:" << std::endl;
    // std::cout << c_a << std::endl;
    // std::cout << std::endl;

    // std::cout << "Cb:" << std::endl;
    // std::cout << c_b << std::endl;
    // std::cout << std::endl;
    // std::cout << std::endl;


    // std::cout << "Ea:" << std::endl;
    // std::cout << e_a * eV << std::endl;
    // std::cout << std::endl;

    // std::cout << "Eb:" << std::endl;
    // std::cout << e_b * eV << std::endl;
    // std::cout << std::endl;
    // std::cout << std::endl;

    // std::cout << "Ending Densities:" << std::endl;
    // std::cout << p_a << std::endl;
    // std::cout << std::endl;
    // std::cout << p_b << std::endl;
    // std::cout << std::endl;
    // std::cout << "----------------" << std::endl;


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
  nuclear_repulsion = third_sum_term;
  electron_energy = first_sum_term + second_sum_term;

  //std::cout << "Nuclear Repulsion Energy is " << third_sum_term  * eV << std::endl;
  //std::cout << "Electron Energy is " << (first_sum_term + second_sum_term) * eV << std::endl;


    if (SCF_cycle > 100) {
      break;
    }
  }

  std::cout << "Nuclear Repulsion Energy=" << nuclear_repulsion * eV << std::endl;
  std::cout << "Electron Energy=" << electron_energy * eV << std::endl;
  std::cout << "Total Energy=" << energy * eV << std::endl;

  /////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////

  // GRADIENT!!!!!!! HW5

  std::unordered_map<int, double> CNDO_2_BETA = get_CNDO_2_BETA();

  //Overlap Term
  MatrixXd overlap_term(3, sim.AOs.size() * sim.AOs.size());
  // MatrixXd overlap_term_coef_X(3, sim.AOs.size() * sim.AOs.size());
  MatrixXd overlap_term_coef_X(sim.AOs.size(), sim.AOs.size());
  overlap_term.fill(0.0);
  overlap_term_coef_X.fill(0.0);

  for (int u=0; u < sim.AOs.size(); u++) {
    for (int v=0; v < sim.AOs.size(); v++) {

      AO& AO_u = sim.AOs[u];
      AO& AO_v = sim.AOs[v];

      int atom_u = AO_u.owning_atom_idx;
      int atom_v = AO_v.owning_atom_idx;

      double Beta_u = CNDO_2_BETA[sim.elements[atom_u]];
      double Beta_v = CNDO_2_BETA[sim.elements[atom_v]];

      int col = u * sim.AOs.size() + v;
      // std::cout << "u=" << u << ", v=" << v << ", col=" << col << std::endl;

      overlap_term_coef_X(u,v) = (Beta_u + Beta_v) * (p_a(u,v) +  p_b(u,v));
      overlap_term_coef_X(u,v) = (Beta_u + Beta_v) * (p_a(u,v) +  p_b(u,v));
      overlap_term_coef_X(u,v) = (Beta_u + Beta_v) * (p_a(u,v) +  p_b(u,v));


      if (atom_u != atom_v) {

        double outer_xterm = 0;
        double outer_yterm = 0;
        double outer_zterm = 0;

        for (int k=0; k < AO_u.primatives.size(); k++) {
          for (int l=0; l < AO_v.primatives.size(); l++) {

            _Primative& p_u = AO_u.primatives[k];
            _Primative& p_v = AO_v.primatives[l];


            // std::cout << "u x coord=" << AO_u.center(0) << std::endl;
            // std::cout << "v x coord=" << AO_v.center(0) << std::endl;

            // std::cout << "u alpha=" << p_u.alpha << std::endl;
            // std::cout << "v alpha=" << p_v.alpha << std::endl;

            // std::cout << "u lmn=" << AO_u.l << std::endl;
            // std::cout << "v lmn=" << AO_v.l << std::endl;


            double x_overlap = __inner_overlap(p_u.alpha, p_v.alpha, AO_u.center(0), AO_v.center(0), AO_u.l, AO_v.l);
            double y_overlap = __inner_overlap(p_u.alpha, p_v.alpha, AO_u.center(1), AO_v.center(1), AO_u.l, AO_v.l);
            double z_overlap = __inner_overlap(p_u.alpha, p_v.alpha, AO_u.center(2), AO_v.center(2), AO_u.l, AO_v.l);

            double dx = __inner_overlap_derivative(p_u.alpha, p_v.alpha, AO_u.center(0), AO_v.center(0), AO_u.l, AO_v.l);
            double dy = __inner_overlap_derivative(p_u.alpha, p_v.alpha, AO_u.center(1), AO_v.center(1), AO_u.l, AO_v.l);
            double dz = __inner_overlap_derivative(p_u.alpha, p_v.alpha, AO_u.center(2), AO_v.center(2), AO_u.l, AO_v.l);

            // std::cout << "Ix=" << x_overlap << std::endl;
            // std::cout << "Iy=" << y_overlap << std::endl;
            // std::cout << "Iz=" << z_overlap << std::endl;

            // std::cout << "dIx=" << dx << std::endl;
            // std::cout << "dIy=" << dy << std::endl;
            // std::cout << "dIz=" << dz << std::endl;


            double xterm = dx * y_overlap * z_overlap;
            double yterm = x_overlap * dy * z_overlap;
            double zterm = x_overlap * y_overlap * dz;

            // std::cout << "x y z terms: " << xterm << " " << yterm << " " << zterm << std::endl;

            outer_xterm += p_u.N * p_u.coef * p_v.N * p_v.coef * xterm;
            outer_yterm += p_u.N * p_u.coef * p_v.N * p_v.coef * yterm;
            outer_zterm += p_u.N * p_u.coef * p_v.N * p_v.coef * zterm;

            // std::cout << "coef " << p_u.N * p_u.coef << " and " << p_v.N * p_v.coef << std::endl;

          }
        }

        overlap_term(0, col) =  outer_xterm;
        overlap_term(1, col) =  outer_yterm;
        overlap_term(2, col) =  outer_zterm;

      }
    }
  }

  std::cout << std::endl;
  std::cout << "Suv_RA" << std::endl;
  std::cout << overlap_term << std::endl;
  std::cout << std::endl;
  std::cout << "coef Xuv" << std::endl;
  std::cout << overlap_term_coef_X * eV << std::endl;
  std::cout << std::endl;


  // Gamma Term
  // FOUR NESTED LOOPS - SEE LAB SLIDES
  //
  // NEED TO BE CAREFUL TO ONLY HAVE S ORBTIALS AGAIN
  int n_s_orbitals = 0;
  for (AO& ao : sim.AOs) {
    if (ao.L() == 0) {
      n_s_orbitals++;
    }
  }

  MatrixXd gamma_term(3, n_s_orbitals * n_s_orbitals);
  gamma_term.fill(0.0);

  int i=0;
  for (AO& ao1 : sim.AOs) {

    if (ao1.L() == 0) {

      int j=0;
      for (AO& ao2 : sim.AOs) {

        if (ao2.L() == 0) {


          int col = i * n_s_orbitals + j;

          //gamma(i, j) = _Eval_2ei_sAO(ao1, ao2);
          // DERIVATIVE HERE

          for (int lp=0; lp < ao1.primatives.size(); lp++) {
            for (int l=0; l < ao1.primatives.size(); l++) {
              for (int kp=0; kp < ao2.primatives.size(); kp++) {
                for(int k=0; k < ao2.primatives.size(); k++) {

                  double lp_fact = ao1.primatives[lp].coef * ao1.primatives[lp].N;
                  double l_fact = ao1.primatives[l].coef * ao1.primatives[l].N;
                  double kp_fact = ao2.primatives[kp].coef * ao2.primatives[kp].N;
                  double k_fact = ao2.primatives[k].coef * ao2.primatives[k].N;

                  double sigma_a = 1.0 / (ao1.primatives[lp].alpha + ao1.primatives[l].alpha);
                  double sigma_b = 1.0 / (ao2.primatives[kp].alpha + ao2.primatives[k].alpha);

                  double Ua = pow((M_PI * sigma_a), 1.5);
                  double Ub = pow((M_PI * sigma_b), 1.5);

                  double Vsq = 1.0 / (sigma_a + sigma_b);

                  VectorXd displacement = (ao1.center - ao2.center);
                  // I think I need to keep separate for the x,y,z terms?? 
                  // like dont take the norm to get distance?

                  // but do need distance to tell if same atom
                  double distance = displacement.norm();
                  // std::cout << "distance=" << distance << std::endl;
                  double Rab_sq = displacement.dot(displacement);
                  // std::cout << "Rab_sq=" << Rab_sq << std::endl;

                  // IF ON SAME ATOM
                  if (distance == 0.0) {
                    // ALL ZERO
                    gamma_term(0, col) = 0;
                    gamma_term(1, col) = 0;
                    gamma_term(2, col) = 0;
                  } else {

                  //VectorXd T = Vsq * displacement.array().pow(2.0); // keep each dim
                  double T = Vsq * Rab_sq;

                  double x_term = ((Ua * Ub) * (displacement(0))) / pow(distance, 2);
                  double y_term = ((Ua * Ub) * (displacement(1))) / pow(distance, 2);
                  double z_term = ((Ua * Ub) * (displacement(2))) / pow(distance, 2);

                  double second_term = -(erf(sqrt(T))/ distance) + ((2 * sqrt(Vsq)) / sqrt(M_PI)) * exp(-T);
                  // std::cout << "second term=" << second_term << std::endl;


                  gamma_term(0, col) += l_fact * lp_fact * k_fact * kp_fact * x_term * second_term;
                  gamma_term(1, col) +=  l_fact * lp_fact * k_fact * kp_fact * y_term * second_term;
                  gamma_term(2, col) +=  l_fact * lp_fact * k_fact * kp_fact * z_term * second_term;


                  // std::cout << "l_fact=" << l_fact << std::endl;
                  // std::cout << "lp_fact=" << lp_fact << std::endl;
                  // std::cout << "k_fact=" << k_fact << std::endl;
                  // std::cout << "kp_fact=" << kp_fact << std::endl;

                  // std::cout << "sigma_a=" << sigma_a << std::endl;
                  // std::cout << "sigma_b=" << sigma_b << std::endl;

                  // std::cout << "Ua=" << Ua << std::endl;
                  // std::cout << "Ub=" << Ub << std::endl;

                  // std::cout << "Vsq=" << Vsq << std::endl;
                  // std::cout << "T=" << T << std::endl;

                  // std::cout << "xterm=" << x_term << std::endl;
                  // std::cout << "yterm=" << y_term << std::endl;
                  // std::cout << "zterm=" << z_term << std::endl;

                  // std::cout << "xterm_err=" << x_term_with_err_fnc << std::endl;
                  // std::cout << "yterm_err=" << y_term_with_err_fnc << std::endl;
                  // std::cout << "zterm_err=" << z_term_with_err_fnc << std::endl;

                  // std::cout << "[0] x=" << x_term * second_term << std::endl;
                  // std::cout << "[0] y=" << y_term * second_term << std::endl;
                  // std::cout << "[0] z=" << z_term * second_term << std::endl;


                  }
                }
              }
            }
          }

          j++;
        }
      }

      i++;
    }
  }

  MatrixXd gamma_term_coef_Y(sim.n_atoms, sim.n_atoms);
  gamma_term_coef_Y.fill(0.0);

  for (int i=0; i < sim.n_atoms; i++) {
    for (int j=0; j < sim.n_atoms; j++) {

      // same code used to get total e density in regular gamma 
      double p_aa_tot = 0.0;
      for (int ao_idx=0; ao_idx < sim.AOs.size(); ao_idx++) {
        if (sim.AOs[ao_idx].owning_atom_idx == i) {
          p_aa_tot += p_a(ao_idx, ao_idx);
          p_aa_tot += p_b(ao_idx, ao_idx);
        }
      }

      // same code used to get total e density in regular gamma 
      double p_bb_tot = 0.0;
      for (int ao_idx=0; ao_idx < sim.AOs.size(); ao_idx++) {
        if (sim.AOs[ao_idx].owning_atom_idx == j) {
          p_bb_tot += p_a(ao_idx, ao_idx);
          p_bb_tot += p_b(ao_idx, ao_idx);
        }
      }

      Element A_element = pt.get_element(sim.elements[i]);
      Element B_element = pt.get_element(sim.elements[j]);

      double A_Z = A_element.z - (A_element.n_e - A_element.n_v_e);
      double B_Z = B_element.z - (B_element.n_e - B_element.n_v_e);

      double sum=0;

      for (int k=0; k < sim.AOs.size(); k++) {
        for (int l=0; l < sim.AOs.size(); l++) {

          if (sim.AOs[k].owning_atom_idx == i && sim.AOs[l].owning_atom_idx == j) {
            sum += pow(p_a(k,l),2) + pow(p_b(k,l),2);

          }
        }
      }

      gamma_term_coef_Y(i,j) = p_aa_tot * p_bb_tot - B_Z * p_aa_tot - A_Z * p_bb_tot - sum;

    }
  }


  std::cout << std::endl;
  std::cout << "gradient (gamma part)" << std::endl;
  std::cout << gamma_term * eV << std::endl;
  std::cout << std::endl;

  std::cout << "Yab" << std::endl;
  std::cout << gamma_term_coef_Y << std::endl;
  std::cout << std::endl;


  MatrixXd grad_term1(3, sim.n_atoms);
  MatrixXd grad_term2(3, sim.n_atoms);

  grad_term1.fill(0);
  grad_term2.fill(0);

  for (int i=0; i < sim.AOs.size(); i++) {
    for (int j=0; j < sim.AOs.size(); j++) {


      if (i != j) {

        grad_term1(0, sim.AOs[i].owning_atom_idx) += overlap_term_coef_X(i,j) * overlap_term(0, i*sim.AOs.size() + j);
        grad_term1(1, sim.AOs[i].owning_atom_idx) += overlap_term_coef_X(i,j) * overlap_term(1, i*sim.AOs.size() + j);
        grad_term1(2, sim.AOs[i].owning_atom_idx) += overlap_term_coef_X(i,j) * overlap_term(2, i*sim.AOs.size() + j);

      }
    }
  }

  for (int i=0; i < sim.n_atoms; i++) {
    for (int j=0; j < sim.n_atoms; j++) {

      if (i != j) {

        grad_term2(0, i) += gamma_term_coef_Y(i,j) * gamma_term(0, i*sim.n_atoms + j);
        grad_term2(1, i) += gamma_term_coef_Y(i,j) * gamma_term(1, i*sim.n_atoms + j);
        grad_term2(2, i) += gamma_term_coef_Y(i,j) * gamma_term(2, i*sim.n_atoms + j);

      }
    }
  }

  MatrixXd grad = grad_term1 + grad_term2;

  std::cout << "e term 1" << std::endl;
  std::cout << grad_term1 * eV << std::endl;
  std::cout << std::endl;
  std::cout << "e term 2" << std::endl;
  std::cout << grad_term2 * eV << std::endl;
  std::cout << std::endl;

  std::cout << "e grad total" << std::endl;
  std::cout << grad * eV << std::endl;
  std::cout << std::endl;


  // Nuclear Term
  MatrixXd nuclear_term(3, sim.n_atoms);
  nuclear_term.fill(0.0);

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

        VectorXd diff = A_center - B_center;
        VectorXd denom = diff.array().abs().pow(3.0);
        denom = denom.array() + 1e-16; // small eps so we dont div by zero
        VectorXd num = -1 * (A_Z * B_Z) * diff;

        nuclear_term(0, A) += num(0) / denom(0);
        nuclear_term(1, A) += num(1) / denom(1);
        nuclear_term(2, A) += num(2) / denom(2);

        nuclear_term(0, B) -= num(0) / denom(0);
        nuclear_term(1, B) -= num(1) / denom(1);
        nuclear_term(2, B) -= num(2) / denom(2);


      }
    }
  }

  std::cout << std::endl;
  std::cout << "gradient (Nuclear part)" << std::endl;
  std::cout << nuclear_term * eV << std::endl;
  std::cout << std::endl;

  std::cout << "TOTAL GRAD" << std::endl;
  std::cout << eV * (grad + nuclear_term) << std::endl;
  std::cout << std::endl;



  return 0;
}


