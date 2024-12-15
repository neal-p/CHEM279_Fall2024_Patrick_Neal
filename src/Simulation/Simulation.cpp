#include "Simulation/Simulation.hpp"

Simulation::Simulation(std::string filename) {

  std::ifstream f;
  f.open(filename);

  // Check input file exists
  if (!f.is_open()) {
    throw InputFileException("Input file " + filename + " does not exist!");
  }

  // Read input format

  /* n_atoms charge
     ELEMENT x y z
     ELEMENT x y z
     ...
  */

  PeriodicTable PT;

  std::string line;
  std::getline(f, line);
  std::istringstream iss(line);

  int n_atoms, charge;

  if (iss >> n_atoms >> charge) {
      this->n_atoms = n_atoms;
      this->charge = charge;

  } else {
      throw InputFileException("Input file " + filename + 
                               " must start with '<n_atoms(int)> <charge(int)>'"
                               " on first line. n_atoms must be >=1."
                               " However, read in: " + line);
  }

  if (this->n_atoms < 1) {
    throw InputFileException("Input file " + filename + 
                               " must start with '<n_atoms(int)> <charge(int)>'"
                               " on first line. n_atoms must be >=1."
                               " However, read in: " + line);
  }

  // Take charge into account for n_e and n_e_v (assume we haven't lost/gained a core e-"
  //    + charge means loosing e-
  //    - charge means adding e-
  this->n_e = -1 * this->charge;
  this->n_v_e = -1 * this->charge;

  // Read Atoms into elements and xyz arrays
  this->xyz = MatrixXd::Zero(this->n_atoms, 3);

  for (int i=0; i < this->n_atoms; i++) {

    std::getline(f, line);
    std::istringstream iss(line);
    int element;
    double x, y, z;

    if (iss >> element >> x >> y >> z) {

      Element e = PT.get_element(element);
      this->n_e += e.n_e;
      this->n_v_e += e.n_v_e;

      this->elements.push_back(element);
      this->xyz(i,0) = x;
      this->xyz(i,1) = y;
      this->xyz(i,2) = z;

      
    } else {
      throw InputFileException("Input file " + filename + 
                               " must provide atom coordinates starting on the second line\n" + 
                               "The format must be '<atomic_number> <x> <y> <z>' separated by whitespace\n\n" +
                               "The following line is incorrect: " + line);
    }
  }

  // We have no way of inupting a spin multiplicity
  // So for now, we assume either we have all paired electrons
  // or that we have 1 extra
  // We will never have two unpaired electrons
  this->p = this->n_v_e / 2;
  this->q = this->n_v_e / 2;

  if (this->p + this->q < this->n_v_e) {
    this->p++;
  }
}


void Simulation::load_basis(std::string basis_file) {

  BasisSet bs = read_basis_set(basis_file);

  for (int i=0; i < this->n_atoms; i++) {

      int z = this->elements[i];
      Array3d center = {this->xyz(i,0), this->xyz(i,1), this->xyz(i, 2)};

      auto it = bs.find(z);
      if (it == bs.end()) {
          std::cout << "Element " << z << " not defined in basis set!" << std::endl;
          exit(1);
      } else {

          // Store a copy of each basis function with proper center etc
          // as an AO
          for (AO new_ao : it->second) {
              new_ao.center = center;
              new_ao.owning_atom_idx = i;

              this->AOs.push_back(new_ao);
          }
      }
  }
}



Simulation::Simulation(const Simulation& sim) {
  this->n_atoms = sim.n_atoms;
  this->charge = sim.charge;
  this->n_e = sim.n_e;
  this->n_v_e = sim.n_v_e;
  this->elements = sim.elements;
  this->xyz = sim.xyz;
  this->AOs = sim.AOs;
}


MatrixXd Simulation::Overlap() {

  MatrixXd S(this->AOs.size(), this->AOs.size());

  for (int i=0; i < this->AOs.size(); i++) {
      for (int j=0; j < this->AOs.size(); j++) {

          S(i,j) = AO_overlap(this->AOs[i], this->AOs[j]);

      }
  }

  return S;
}



std::ostream& operator<<(std::ostream& os, const Simulation& sim) {

  os << "Simulation: " << std::endl;
  os << "    " << "n_atoms=" << sim.n_atoms << ", charge=" << sim.charge << ", n_AOs= " << sim.AOs.size() << ", n_e-=" << sim.n_e << ", n_valence_e=" << sim.n_v_e << std::endl;
  for (int i=0; i<sim.n_atoms; i++) {
    os << "    " << std::left << std::setw(ATOMIC_NUMBER_CHARS) << sim.elements[i];

    for (int j=0; j<3; j++) {
      os << std::right << std::setw(COORD_CHARS) << std::fixed << std::setprecision(COORD_PRECISION) << sim.xyz(i,j);
    }

    os << std::endl;
  }

  return os;
}




