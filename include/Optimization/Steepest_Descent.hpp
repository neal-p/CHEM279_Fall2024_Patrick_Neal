#pragma once
#include <iostream>

/////////////////////

#include "Optimization/Optimization.hpp"

Simulation steepest_descent(Simulation& sim, 
                            EnergyFunc E, GradFunc dE,
                            double step_size, int max_steps, double tol);


Simulation steepest_descent_line_search(Simulation& sim,
                                        EnergyFunc E, GradFunc dE,
                                        int max_steps, double tol); 


