#pragma once

#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <iomanip>

#include "Gaussian/Gaussian.hpp"
#include "PeriodicTable/PeriodicTable.hpp"

#include "eigen-3.4.0/Eigen/Dense"
using namespace Eigen;


////////////////////////////////////////////////////////////

#define ATOMIC_NUMBER_CHARS 5
#define COORD_CHARS 10
#define COORD_PRECISION 5

////////////////////////////////////////////////////////////


class InputFileException {

public:
  std::string message;

  InputFileException(std::string message) : message(message) {}

};


class Simulation {

public:

  int n_atoms, charge, n_e, n_v_e, p, q;
  std::vector<int> elements;
  MatrixXd xyz;
  std::vector<AO> AOs;

  Simulation(std::string filename);
  Simulation(const Simulation& sim);

  void load_basis(std::string basis_file);
  MatrixXd Overlap();

  friend std::ostream& operator<<(std::ostream& os, const Simulation& sim);
};





