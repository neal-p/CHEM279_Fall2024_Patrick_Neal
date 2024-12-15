#pragma once

#include "Gradients/Gradients.hpp"
#include "Simulation/Simulation.hpp"
#include "Energy/Energy.hpp"

MatrixXd central_difference(Simulation& sim, EnergyFunc E, double h);
MatrixXd forward_difference(Simulation& sim, EnergyFunc E, double h);
