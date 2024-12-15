#include "Gaussian/Gaussian.hpp"

// Copy constructor bc we want
// to COPY an AO from the basis set
// any time we need it
AO::AO(const AO& ao) {
    this->center = ao.center;
    this->owning_atom_idx = ao.owning_atom_idx;
    this->l = ao.l;
    this->m = ao.m;
    this->n = ao.n;
    this->primatives = ao.primatives;
}

// Same with primatives!
_Primative::_Primative(const _Primative& p) {
    this->coef = p.coef;
    this->alpha = p.alpha;
    this->N = p.N;
}


void AO::add_primative(double coef, double alpha) {

    this->primatives.emplace_back(coef, alpha);
    double overlap = PrimativeOverlap(this->primatives.back(), this->primatives.back(), this->center, this->center, this->l, this->l, this->m, this->m, this->n, this->n);
    this->primatives.back().N = 1.0 / (sqrt(overlap));
}

BasisSet read_basis_set(std::string file) {
    BasisSet bs;
    PeriodicTable PT; // Convert from atomic symbol to z

    std::ifstream infile(file);
    if (!infile) {
        std::cout << "Cant read file" << std::endl;
        return bs;
    }

    std::string line;
    std::string element;
    int z;
    int num_primatives;
    bool element_block_start = true;

    while (std::getline(infile, line)) {

        // Skip empty or comment lines
        if (line.empty() || line.substr(0, 1) == "!") {
            continue;

        // If next line is '****' then start a new element
        } else if ( line.substr(0, 4) == "****") {
            element_block_start = true;
            continue;
        }

        if (element_block_start) {

            // Element block starts with Element STRING
            std::istringstream iss(line);
            int number_idk_what_means;
            iss >> element >> number_idk_what_means;
            element_block_start = false;

            z = PT.get_element(element).z;
            bs[z];
            continue;
        }

        // Not empty, comment, new element indicator
        // Not element start block
        
        // So next line should be start of orbital
        // definitions
        // Next line is <type> <number> <scale>

        std::istringstream iss(line);
        std::string type;
        double scale;
        iss >> type >> num_primatives >> scale;

        // Handle each type separately
        if (type == "S") {

            // S orbital - 0 angular momentum
            AO ao = AO(0, 0, 0);

            // Get line for each primative
            for (int i=0; i < num_primatives; i++) {
                std::getline(infile, line);
                std::istringstream line_stream(line);

                double alpha, coef;

                line_stream >> alpha >> coef;
                ao.add_primative(coef, alpha);
            }
            bs[z].push_back(ao);

        } else if (type == "P") {

            // P orbital - 1 angular momentum
            AO px(1, 0, 0);
            AO py(0, 1, 0);
            AO pz(0, 0, 1);


            // Get line for each primative
            for (int i=0; i < num_primatives; i++) {
                std::getline(infile, line);
                std::istringstream line_stream(line);

                double alpha, coef;

                line_stream >> alpha >> coef;

                px.add_primative(coef, alpha);
                py.add_primative(coef, alpha);
                pz.add_primative(coef, alpha);
            }

            bs[z].push_back(px);
            bs[z].push_back(py);
            bs[z].push_back(pz);

        } else if (type == "SP") {

            // S orbital - 0 angular momentum
            AO s2(0, 0, 0);

            // P orbital - 1 angular momentum
            AO px(1, 0, 0);
            AO py(0, 1, 0);
            AO pz(0, 0, 1);


            // Get line for each primative
            for (int i=0; i < num_primatives; i++) {
                std::getline(infile, line);
                std::istringstream line_stream(line);

                double alpha, S_coef, P_coef;
                line_stream >> alpha >> S_coef >> P_coef;

                s2.add_primative(S_coef, alpha);
                px.add_primative(P_coef, alpha);
                py.add_primative(P_coef, alpha);
                pz.add_primative(P_coef, alpha);
                
            }

            bs[z].push_back(s2);
            bs[z].push_back(px);
            bs[z].push_back(py);
            bs[z].push_back(pz);


        } else {
            std::cout << "type " << type << " not implemented" << std::endl;
            // Get line for each primative
            for (int i=0; i < num_primatives; i++) {
                std::getline(infile, line);
            }
        }
    }

    return bs;
}


// double __inner_overlap(double alpha, double beta, double Ra, double Rb, double Rp, int la, int lb) {
double __inner_overlap(double alpha, double beta, double Ra, double Rb, int la, int lb) {

    double sum = 0.0;
    double new_coord = (alpha * Ra + beta * Rb) / (alpha + beta);

    for (int i=0; i <= la; i++) {
      for (int j=0; j <= lb; j++) {

          if ((i + j) % 2 == 0) {

            int binom_A = binom(la, i);
            int binom_B = binom(lb, j);
            sum += binom_A * binom_B * double_factorial(i + j -1) * pow(new_coord - Ra, la-i) * pow(new_coord - Rb, lb-j) / pow(2.0 * (alpha + beta), double(i+j) / 2.0);
          }
        }
      }

  double root_term = sqrt(M_PI / (alpha + beta));
  //double expo_term = exp(-(alpha * beta / (alpha + beta)) * (Ra - Rb));
  double expo_term = exp(-alpha * beta * (Ra - Rb) * (Ra - Rb) / (alpha + beta));

  return sum * root_term * expo_term;
}

double __inner_overlap_derivative(double alpha, double beta, double Ra, double Rb, int la, int lb) {

  double first_term = -la * __inner_overlap(alpha, beta, Ra, Rb, la -1, lb);
  double second_term = 2 * alpha * __inner_overlap(alpha, beta, Ra, Rb, la+1, lb);
  return first_term + second_term;

}



double AO_overlap(AO& A, AO& B) {

    double sum = 0.0;

    // All A by All B
    int A_idx=0;
    for (_Primative& A_p : A.primatives) {

        int B_idx=0;
        for (_Primative& B_p : B.primatives) {

            double overlap = PrimativeOverlap(A_p, B_p, A.center, B.center, A.l, B.l, A.m, B.m, A.n, B.n);
            sum += A_p.coef * B_p.coef * A_p.N * B_p.N * overlap ;
        }
    }

    return sum;
}

double PrimativeOverlap(_Primative& A, _Primative& B, Array3d A_center, Array3d B_center, int A_l, int B_l, int A_m, int B_m, int A_n, int B_n) {

    // Exponential terms
    // Array3d new_xyz = (A.alpha * A_center + B.alpha * B_center) / (A.alpha + B.alpha);
    // double root_term = sqrt(M_PI / (A.alpha + B.alpha));
    // Array3d expo_term = exp(-(A.alpha * B.alpha / (A.alpha + B.alpha)) * (A_center - B_center).pow(2.0));

    // Angular momentum terms
    // double x_overlap = __inner_overlap(A.alpha, B.alpha, A_center(0), B_center(0), new_xyz(0), A_l, B_l);
    // double y_overlap = __inner_overlap(A.alpha, B.alpha, A_center(1), B_center(1), new_xyz(1), A_m, B_m);
    // double z_overlap = __inner_overlap(A.alpha, B.alpha, A_center(2), B_center(2), new_xyz(2), A_n, B_n);
  
    double x_overlap = __inner_overlap(A.alpha, B.alpha, A_center(0), B_center(0), A_l, B_l);
    double y_overlap = __inner_overlap(A.alpha, B.alpha, A_center(1), B_center(1), A_m, B_m);
    double z_overlap = __inner_overlap(A.alpha, B.alpha, A_center(2), B_center(2), A_n, B_n);

    double overlap = x_overlap * y_overlap * z_overlap;

	//overlap *= expo_term.prod();
    //overlap *= pow(root_term, 3.0);

	return overlap;
}


