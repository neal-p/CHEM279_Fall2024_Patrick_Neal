#pragma once

#include "Energy/Energy.hpp"
#include "Simulation/Simulation.hpp"

/////////////////////

#define e_au 5.29
#define o_au 2.951

/////////////////////


double LJ_energy(Simulation& sim);

MatrixXd LJ_gradient_analytical(Simulation& sim);

