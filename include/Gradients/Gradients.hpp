#pragma once

#include "eigen-3.4.0/Eigen/Dense"
#include "Energy/Energy.hpp"
#include "Simulation/Simulation.hpp"

typedef MatrixXd(*GradFunc)(Simulation&);

