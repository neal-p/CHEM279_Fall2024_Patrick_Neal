#include "CNDO_2/CNDO_2.hpp"

std::unordered_map<int, double> get_CNDO_2_BETA() {

  std::unordered_map<int, double> CNDO_2_BETA;
  CNDO_2_BETA[1] = -9.0 / eV;
  CNDO_2_BETA[6] = -21.0 / eV;
  CNDO_2_BETA[7] = -25.0 / eV;
  CNDO_2_BETA[8] = -31.0 / eV;
  CNDO_2_BETA[9] = -39.0 / eV;

  return CNDO_2_BETA;
}

std::unordered_map<int, double> get_CNDO_2_s() {
  std::unordered_map<int, double> CNDO_2_s;
  CNDO_2_s[1] = 7.176 / eV;
  CNDO_2_s[6] = 14.051 / eV;
  CNDO_2_s[7] = 19.316 / eV;
  CNDO_2_s[8] = 25.390 / eV;
  CNDO_2_s[9] = 32.272 / eV;

  return CNDO_2_s;
}

std::unordered_map<int, double> get_CNDO_2_p() {
  std::unordered_map<int, double> CNDO_2_p;
  CNDO_2_p[6] = 5.572 / eV;
  CNDO_2_p[7] = 7.275 / eV;
  CNDO_2_p[8] = 9.111 / eV;
  CNDO_2_p[9] = 11.080 / eV;

  return CNDO_2_p;
}





 
double _I2e_pG(Array3d& Ra, Array3d& Rb, double sigmaA, double sigmaB) {

  double Ua = pow(M_PI * sigmaA, (3.0 / 2.0)); // 3.11
  double Ub = pow(M_PI * sigmaB, (3.0 / 2.0)); // 3.11

  double U =  Ua * Ub; // 3.8 and 3.11
  double V2 = 1.0 / (sigmaA + sigmaB); // 3.9
  
  Vector3d Rab = Ra - Rb;

  double Rd = Rab.norm();; // norm of vector difference - this is just a distance
  
  if (Rd == 0.0) {
    // Centers are the same
    return U * 2.0 * sqrt(V2 / M_PI); // 3.15
  }

  double srT = sqrt(V2 * pow(Rd, 2.0)); // 3.7 sqrt

  // double erf_srT = (2.0 / sqrt(M_PI)) * 
  
  double result = U * sqrt(1.0 / pow(Rd, 2.0)) * std::erf(srT); // 3.14
  
  return result;
}

double _Eval_2ei_sAO(AO& ao1, AO& ao2) {

  int len = ao1.primatives.size();
  if (ao1.primatives.size() != ao2.primatives.size()) {
   std::cout << "ao1 and ao2 have different number of primatives!!" << std::endl;
   exit(1);
  }


  // std::cout << "AO1:" << std::endl;
  // std::cout << "    alpha=" << ao1.primatives[0].alpha << std::endl;
  // std::cout << "    alpha=" << ao1.primatives[1].alpha << std::endl;
  // std::cout << "    alpha=" << ao1.primatives[2].alpha << std::endl;

  // std::cout << std::endl;

  // std::cout << "    coef=" << ao1.primatives[0].coef << std::endl;
  // std::cout << "    coef=" << ao1.primatives[1].coef << std::endl;
  // std::cout << "    coef=" << ao1.primatives[2].coef << std::endl;

  // std::cout << std::endl;

  // std::cout << "    N=" << ao1.primatives[0].N << std::endl;
  // std::cout << "    N=" << ao1.primatives[1].N << std::endl;
  // std::cout << "    N=" << ao1.primatives[2].N << std::endl;

  // std::cout << std::endl;
  // std::cout << std::endl;

  // std::cout << "AO2:" << std::endl;
  // std::cout << "    alpha=" << ao2.primatives[0].alpha << std::endl;
  // std::cout << "    alpha=" << ao2.primatives[1].alpha << std::endl;
  // std::cout << "    alpha=" << ao2.primatives[2].alpha << std::endl;

  // std::cout << std::endl;

  // std::cout << "    coef=" << ao2.primatives[0].coef << std::endl;
  // std::cout << "    coef=" << ao2.primatives[1].coef << std::endl;
  // std::cout << "    coef=" << ao2.primatives[2].coef << std::endl;

  // std::cout << std::endl;

  // std::cout << "    N=" << ao2.primatives[0].N << std::endl;
  // std::cout << "    N=" << ao2.primatives[1].N << std::endl;
  // std::cout << "    N=" << ao2.primatives[2].N << std::endl;



  double gamma = 0.0;


  // Double loop over K and l - 3.3
  for (int k1=0; k1 < len; k1++) {
    for (int k2=0; k2 < len; k2++) {

      double alpha_k_1 = ao1.primatives[k1].alpha;
      double alpha_k_2 = ao1.primatives[k2].alpha;
      double coef_k_1 = ao1.primatives[k1].coef;
      double coef_k_2 = ao1.primatives[k2].coef;
      double N_k_1 = ao1.primatives[k1].N;
      double N_k_2 = ao1.primatives[k2].N;

      double sigmaA =  1.0 / (alpha_k_1 + alpha_k_2); // 3.10
      
      
      for (int j1=0; j1 < len; j1++) {
        for (int j2=0; j2 < len; j2++) {

          double alpha_j_1 = ao2.primatives[j1].alpha;
          double alpha_j_2 = ao2.primatives[j2].alpha;
          double coef_j_1 = ao2.primatives[j1].coef;
          double coef_j_2 = ao2.primatives[j2].coef;
          double N_j_1 = ao2.primatives[j1].N;
          double N_j_2 = ao2.primatives[j2].N;


          double sigmaB =  1.0 / (alpha_j_1 + alpha_j_2); // 3.10
 
          double I2e = _I2e_pG(ao1.center, ao2.center, sigmaA, sigmaB);
          gamma += coef_k_1 * N_k_1 * coef_k_2 * N_k_2 * coef_j_1 * N_j_1 * coef_j_2 * N_j_2 * I2e; // 3.3


        }
      }


    }
  }

  return gamma;

}



MatrixXd CalculateGamma(Simulation& sim) {

  int n_s_orbitals = 0;
  for (AO& ao : sim.AOs) {
    if (ao.L() == 0) {
      n_s_orbitals++;
    }
  }

  MatrixXd gamma(n_s_orbitals, n_s_orbitals);

  int i=0;
  for (AO& ao1 : sim.AOs) {

    if (ao1.L() == 0) {

      int j=0;
      for (AO& ao2 : sim.AOs) {

        if (ao2.L() == 0) {
          gamma(i, j) = _Eval_2ei_sAO(ao1, ao2);
          j++;
        }
      }

      i++;
    }
  }

  return gamma;

}


MatrixXd CalculateFock(Simulation& sim, MatrixXd& gamma, MatrixXd& p, MatrixXd& p_other, MatrixXd& overlap) {

  std::unordered_map<int, double> CNDO_2_s = get_CNDO_2_s();
  std::unordered_map<int, double> CNDO_2_p = get_CNDO_2_p();
  std::unordered_map<int, double> CNDO_2_BETA = get_CNDO_2_BETA();

  PeriodicTable pt;

  // Initial guess is all 0s
  MatrixXd Fock(sim.AOs.size(), sim.AOs.size());

  for (int i=0; i < sim.AOs.size(); i++) {
    for (int j=0; j < sim.AOs.size(); j++) {

      // 1.4 for on-diagonal
      // 1.5 for off-diagonal
      
      if (i == j) {

        int a_atom_idx = sim.AOs[i].owning_atom_idx;
        int a_element = sim.elements[a_atom_idx];

        Element a_element_obj = pt.get_element(a_element);
        int a_Z_eff = a_element_obj.z - (a_element_obj.n_e - a_element_obj.n_v_e);

        double param;
        if (sim.AOs[i].L() == 0) {
          param = CNDO_2_s[a_element];
        } else if (sim.AOs[i].L() == 1) {
          param = CNDO_2_p[a_element];
        } else {
          std::cout << "Dont have param for this orbital" << std::endl;
          exit(1);
        }

        double p_aa_tot = 0.0;
        for (int ao_idx=0; ao_idx < sim.AOs.size(); ao_idx++) {
          if (sim.AOs[ao_idx].owning_atom_idx == a_atom_idx) {
            p_aa_tot += p(ao_idx, ao_idx);
            p_aa_tot += p_other(ao_idx, ao_idx);
          }
        }

        double sum = 0;

        for (int b_atom_idx=0; b_atom_idx < sim.n_atoms; b_atom_idx++) {
          if (b_atom_idx != a_atom_idx) {

            double p_bb_tot = 0.0;
            for (int ao_idx=0; ao_idx < sim.AOs.size(); ao_idx++) {
              if (sim.AOs[ao_idx].owning_atom_idx == b_atom_idx) {
                p_bb_tot += p(ao_idx, ao_idx);
                p_bb_tot += p_other(ao_idx, ao_idx);

              }
            }

            int b_element = sim.elements[b_atom_idx];
            double gamma_a_b = gamma(a_atom_idx, b_atom_idx);

            Element b_element_obj = pt.get_element(b_element);
            int b_Z_eff = b_element_obj.z - (b_element_obj.n_e - b_element_obj.n_v_e);

            sum += (p_bb_tot - b_Z_eff) * gamma_a_b;
          }
        }

        Fock(i,j) = -param + ((p_aa_tot - a_Z_eff) - (p(i, j) - 0.5)) * gamma(a_atom_idx, a_atom_idx) + sum;


      } else {

        int atom_idx_u = sim.AOs[i].owning_atom_idx;
        int element_u = sim.elements[atom_idx_u];

        int atom_idx_v = sim.AOs[j].owning_atom_idx;
        int element_v = sim.elements[atom_idx_v];

        double Beta_a = CNDO_2_BETA[element_u];
        double Beta_b = CNDO_2_BETA[element_v];

        double gamma_a_b = gamma(atom_idx_u, atom_idx_v);

        Fock(i,j) = 0.5 * (Beta_a + Beta_b) * overlap(i,j) - p(i,j) * gamma_a_b;

      }
    }
  }


  return Fock;
}



MatrixXd CalculateHcore(Simulation& sim, MatrixXd& gamma, MatrixXd& overlap) {

  std::unordered_map<int, double> CNDO_2_s = get_CNDO_2_s();
  std::unordered_map<int, double> CNDO_2_p = get_CNDO_2_p();
  std::unordered_map<int, double> CNDO_2_BETA = get_CNDO_2_BETA();

  PeriodicTable pt;

  // Initial guess is all 0s
  MatrixXd Hcore(sim.AOs.size(), sim.AOs.size());

  for (int i=0; i < sim.AOs.size(); i++) {
    for (int j=0; j < sim.AOs.size(); j++) {

      // 1.4 for on-diagonal
      // 1.5 for off-diagonal
      
      if (i == j) {

        int a_atom_idx = sim.AOs[i].owning_atom_idx;
        int a_element = sim.elements[a_atom_idx];

        Element a_element_obj = pt.get_element(a_element);
        int a_Z_eff = a_element_obj.z - (a_element_obj.n_e - a_element_obj.n_v_e);

        double param;
        if (sim.AOs[i].L() == 0) {
          param = CNDO_2_s[a_element];
        } else if (sim.AOs[i].L() == 1) {
          param = CNDO_2_p[a_element];
        } else {
          std::cout << "Dont have param for this orbital" << std::endl;
          exit(1);
        }

       double sum = 0;

        for (int b_atom_idx=0; b_atom_idx < sim.n_atoms; b_atom_idx++) {
          if (b_atom_idx != a_atom_idx) {

            int b_element = sim.elements[b_atom_idx];
            double gamma_a_b = gamma(a_atom_idx, b_atom_idx);

            Element b_element_obj = pt.get_element(b_element);
            int b_Z_eff = b_element_obj.z - (b_element_obj.n_e - b_element_obj.n_v_e);

            sum += b_Z_eff * gamma_a_b;
          }
        }

        Hcore(i,j) = -param - (a_Z_eff - 0.5) * gamma(a_atom_idx, a_atom_idx) - sum;


      } else {

        int atom_idx_u = sim.AOs[i].owning_atom_idx;
        int element_u = sim.elements[atom_idx_u];

        int atom_idx_v = sim.AOs[j].owning_atom_idx;
        int element_v = sim.elements[atom_idx_v];

        double Beta_a = CNDO_2_BETA[element_u];
        double Beta_b = CNDO_2_BETA[element_v];

        double gamma_a_b = gamma(atom_idx_u, atom_idx_v);

        Hcore(i,j) = 0.5 * (Beta_a + Beta_b) * overlap(i,j);

      }
    }
  }


  return Hcore;
}






