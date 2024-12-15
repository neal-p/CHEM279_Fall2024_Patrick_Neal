#include "Optimization/Steepest_Descent.hpp"

Simulation steepest_descent(Simulation& sim, 
                            EnergyFunc E, GradFunc dE,
                            double step_size, int max_steps, double tol) {

  int steps = 0;
  MatrixXd grad_i, grad_f;
  double convergence, e_i, e_f;

  // make a deep copy
  Simulation updated_sim = sim;

  // Initial Grad
  grad_i = dE(updated_sim);
  convergence = grad_i.norm();
  e_i = E(updated_sim);

  std::cout << "Initial enerrgy: " << e_i << std::endl;
  std::cout << "Central Difference Force" << std::endl;
  std::cout << grad_i << std::endl;

  // For as long as the gradient is 
  //    larger than tol
  //    and
  //    we have taken less than or equal to max_steps
  while (convergence > tol && steps <= max_steps) {

    updated_sim.xyz.array() += (step_size * grad_i.transpose().array());
    e_f = E(updated_sim);

    if (e_f <= e_i + 1e-5) {
      // Increase step size bc we stepped in good direction
      // and now gradient is smaller
      step_size *= 1.2;

      e_i = e_f;
      grad_i = dE(updated_sim);
      convergence = grad_i.norm();

    } else {
      // Reset to previous position
      updated_sim.xyz.array() -= (step_size * grad_i.transpose().array());

      // Reduce step size bc we stepped to far
      step_size *= 0.5;
    }

    steps++;
  }

  if (convergence > tol && steps >= max_steps) {
    std::cout << "WARNING! optimization not finished" << std::endl;
  } else {
    std::cout << "Optimization finished in " << steps << " steps" << std::endl;
  }

  return updated_sim;
}



Simulation steepest_descent_line_search(Simulation& sim,
                                        EnergyFunc E, GradFunc dE,
                                        int max_steps, double tol) {

  int steps = 0;
  int gr_steps, max_index, nearest_index;
  MatrixXd grad;
  double convergence, e_i, e_f, shortest_bond, distance_thresh, largest_atom_grad, nearest_atom_grad, a, b, d, x1, x2, e_x1, e_x2, step_size, max_atom_grad, nearest_atom_distance, max_displacement, e_b;


  double gr = (pow(5.0, 0.5) - 1.0) / 2.0;

  // make a deep copy
  Simulation updated_sim = sim;

  // Will use atomic distances heuristic
  // to inform step size
  int pairwise_indices = sim.n_atoms * (sim.n_atoms - 1) / 2;
  ArrayXd distances(pairwise_indices);

  // Initial Grad
  grad = dE(updated_sim);
  convergence = grad.norm();
  e_i = E(updated_sim);

  std::cout << "Initial Energy: " << e_i << std::endl;
  std::cout << "Initial Gradient:" << std::endl;
  std::cout << grad << std::endl;

// For as long as the gradient is 
  //    larger than tol
  //    and
  //    we have taken less than or equal to max_steps
  while ((steps < 3) || (convergence > tol && steps <= max_steps)) {

    std::cout << "Optimization Step: " << steps << std::endl;


    VectorXd per_atom_grad_norm = grad.colwise().norm();
    max_atom_grad = per_atom_grad_norm.maxCoeff(&max_index);

    distances = (updated_sim.xyz.rowwise() - updated_sim.xyz.row(max_index)).rowwise().norm();
    distances(max_index) = distances.maxCoeff(); // ignore distance to self when finding minimum
    nearest_atom_distance = distances.minCoeff(&nearest_index);

    nearest_atom_grad = per_atom_grad_norm(nearest_index);

    a = 0.0;  // Lower bound is to not move at all
    b = (3.0 * nearest_atom_distance) / (4.0 * (max_atom_grad + nearest_atom_grad)); // upper bound is 
                                                                                     // what it would take to shorten the
                                                                                     // rij between atom with largest grad
                                                                                     // and its nearest neighbor
                                                                                     // to 1/4 the current rij
    
    // std::cout << "a: " << a << " b: " << b << std::endl;

    int gs_steps = 0;
    while ((abs(b - a) > 1e-8) && gs_steps < 1000000) {
        // std::cout << "    GS: " << gs_steps << std::endl;

      d = gr * (b - a);

      x2 = b - d;
      x1 = a + d;

      // std::cout << "        " << "a:" << a << " x2:" << x2 << " x1:" << x1 << " b:" << b << std::endl;

      updated_sim.xyz.array() += (x2 * grad.transpose().array());
      e_x2 = E(updated_sim);

      // std::cout << "Coords at x2:" << std::endl;
      // std::cout << updated_sim.xyz << std::endl;

      updated_sim.xyz.array() -= (x2 * grad.transpose().array());

      updated_sim.xyz.array() += (x1 * grad.transpose().array());
      e_x1 = E(updated_sim);
      updated_sim.xyz.array() -= (x1 * grad.transpose().array());

      // std::cout << "            e_x2:" << e_x2 << std::endl;
      // std::cout << "            e_x1:" << e_x1 << std::endl;


      if (e_x2 > e_x1) {
        a = x2;
      } else {
        b = x1;
      }

      gs_steps++;
    }

    step_size = (a + b) / 2.0;
    updated_sim.xyz.array() += (step_size * grad.transpose().array());
    e_f = E(updated_sim);


    if (e_f > e_i) {
      std::cout << "SOMETHING IS WRONG! energy increased after gr: " << e_i << " -> " << e_f << std::endl;
    }

    grad = dE(updated_sim);
    convergence = grad.norm();

    // If our energy changes significantly
    // even if our gradient is small
    // we should keep trying to optimize
    // if (abs(e_f - e_i) > 1e-4) {
    //   std::cout << "Gradient is small, but energy change was significant: " << e_i << " -> " << e_f << std::endl;
    //   std::cout << "Continuing optimization..." << std::endl;
    //   convergence += tol;
    // }

    e_i = e_f;

    std::cout << "Updated Coordiantes:" << std::endl;
    std::cout << updated_sim << std::endl;
    std::cout << "Updated Energy: " << e_i << std::endl;

    steps++;
  }


  if (convergence > tol && steps >= max_steps) {
    std::cout << "WARNING! optimization not finished" << std::endl;
  } else {
    std::cout << "Optimization finished in " << steps << " steps" << std::endl;
  }

  std::cout << std::endl << std::endl;
  std::cout << "Final Energy:" << e_i << std::endl;
  std::cout << "Final Gradient: " << std::endl;
  std::cout << grad << std::endl;
  std::cout << "Final Coordinates:" << std::endl;
  std::cout << updated_sim << std::endl;
 
  std::cout << "Final atom-atom distances:" << std::endl;
  std::cout << distances << std::endl;

  return updated_sim;
}
