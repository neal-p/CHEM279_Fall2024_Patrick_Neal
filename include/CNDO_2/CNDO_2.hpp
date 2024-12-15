#pragma once

#include <iostream>
#include <cmath>
#include <unordered_map>
#include "Simulation/Simulation.hpp"
#include "PeriodicTable/PeriodicTable.hpp"

#include "../../external/eigen-3.4.0/Eigen/Dense"
using namespace Eigen;

#define eV 27.211

std::unordered_map<int, double> get_CNDO_2_BETA();
std::unordered_map<int, double> get_CNDO_2_s();
std::unordered_map<int, double> get_CNDO_2_p();


MatrixXd CalculateGamma(Simulation& sim);
MatrixXd CalculateFock(Simulation& sim, MatrixXd& gamma, MatrixXd& p, MatrixXd& p_other, MatrixXd& overlap);
MatrixXd CalculateHcore(Simulation& sim, MatrixXd& gamma, MatrixXd& overlap);
