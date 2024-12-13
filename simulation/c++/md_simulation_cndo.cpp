#include <iostream>
#include <vector>
#include <armadillo>
#include "CNDO.h"
#include "AO.h"
#include "util.h"

struct MDSystem {
    Molecule_basis molecule;
    arma::mat velocities; // Velocities of atoms (3 x num_atoms)
    arma::mat forces;     // Forces on atoms (3 x num_atoms)
    double time_step;     // Time step for integration
    int num_steps;        // Number of integration steps

    MDSystem(Molecule_basis &mol, double dt, int steps) :
        molecule(mol), time_step(dt), num_steps(steps) {
        int num_atoms = molecule.getnum_atoms();
        velocities = arma::zeros(3, num_atoms);
        forces = arma::zeros(3, num_atoms);
    }
};

void compute_forces(MDSystem &system, CNDO &cndo) {
    // Compute the forces using CNDO gradient
    arma::mat gradient = cndo.getGradient();
    system.forces = -gradient; // Force is the negative gradient of the energy
}

void velocity_verlet(MDSystem &system, CNDO &cndo) {
    int num_atoms = system.molecule.getnum_atoms();
    arma::mat positions(3, num_atoms);

    // Extract initial positions from the molecule
    for (size_t i = 0; i < num_atoms; ++i) {
        positions.col(i) = system.molecule.mAtoms[i].mAOs[0].get_R0();
    }

    for (int step = 0; step < system.num_steps; ++step) {
        std::cout << "Step " << step << ":" << std::endl;

        // Update positions: R(t + dt) = R(t) + v(t)*dt + 0.5*F(t)/m*dt^2
        positions += system.velocities * system.time_step + 0.5 * system.forces * std::pow(system.time_step, 2);

        // Update molecule's atomic positions
        for (size_t i = 0; i < num_atoms; ++i) {
            arma::vec temp_position = positions.col(i);
            system.molecule.mAtoms[i].set_R0(temp_position);
        }

        // Recompute molecular data
        // system.molecule.eval_OVmat(system.forces); // Example of recomputation if needed
        // cndo.init(); // Reinitialize the CNDO object with updated molecular data

        // Recompute forces at the updated positions
        compute_forces(system, cndo);

        // Update velocities: v(t + dt) = v(t) + 0.5*[F(t) + F(t + dt)]/m*dt
        system.velocities += 0.5 * (system.forces) * system.time_step;

        // Output current positions and velocities
        std::cout << "Positions:\n" << positions << std::endl;
        std::cout << "Velocities:\n" << system.velocities << std::endl;
    }
}

void geometry_optimization(Molecule_basis &molecule, CNDO &cndo, double threshold = 1e-6, double step_size = 0.01) {
    int max_iterations = 100;
    int num_atoms = molecule.getnum_atoms();
    arma::mat gradient;

    for (int iter = 0; iter < max_iterations; ++iter) {
        // Compute the gradient
        gradient = cndo.getGradient();

        // Compute the norm of the gradient to check convergence
        double grad_norm = arma::norm(gradient, "fro"); // Frobenius norm
        std::cout << "Iteration " << iter << ", Gradient Norm: " << grad_norm << std::endl;

        if (grad_norm < threshold) {
            std::cout << "Converged to equilibrium geometry!" << std::endl;
            break;
        }

        // Update atomic positions using the gradient
        for (int i = 0; i < num_atoms; ++i) {
            arma::vec position = molecule.mAtoms[i].mAOs[0].get_R0(); // Get current position
            position -= step_size * gradient.col(i); // Update position
            molecule.mAtoms[i].set_R0(position);    // Update in molecule
        }

        // Reinitialize CNDO to account for updated positions
        cndo.init();
        cndo.run();
    }
}

// int main() {

//     arma::mat H_basis, C_basis;
//     std::string molecule_file = "H2O.txt";

//     // Load STO-3G basis sets for hydrogen and carbon (example)
//     H_basis.load("H_STO3G.txt");
//     C_basis.load("C_STO3G.txt");

//     // Initialize molecule and CNDO
//     Molecule_basis molecule(molecule_file, H_basis, C_basis);
//     CNDO cndo(molecule, 100, 1e-6);

//     if (cndo.init() != 0) {
//         std::cerr << "Failed to initialize CNDO." << std::endl;
//         return 1;
//     }

//     // Perform geometry optimization
//     geometry_optimization(molecule, cndo);

//     // Output final bond length
// // Output final bond length
//     for (auto &atom : molecule.mAtoms) {
//         std::cout << "Final Position: " << atom.mAOs[0].get_R0().t() << std::endl;
//     }

//     return 0;
// }




int main(int argc, char*argv[]) {
    // Define input files for molecule initialization and basis sets
    // std::string molecule_file = "H2O.txt";
    std::string molecule_file = argv[1];
    arma::mat H_basis, C_basis;

    // Load STO-3G basis sets for hydrogen and carbon (example)
    H_basis.load("H_STO3G.txt");
    C_basis.load("C_STO3G.txt");

    // Initialize the molecule
    Molecule_basis molecule(molecule_file, H_basis, C_basis);

    // Set up the CNDO/2 calculation
    CNDO cndo(molecule, 100, 1e-6); // Max iterations and tolerance

    if (cndo.init() != 0) {
        std::cerr << "Failed to initialize CNDO." << std::endl;
        return 1;
    }

    // std::cout << "Before running SCF" << std::endl;

    // double Energy = cndo.getEnergy();
    // arma::mat gradient = cndo.getGradient();
    // gradient.print("gradient");

    if (cndo.run() != 0) {
        std::cerr << "SCF calculation did not converge." << std::endl;
        return 1;
    }

    // Energy = cndo.getEnergy();
    // gradient = cndo.getGradient();
    // gradient.print("gradient");

    // Set up the MD system
    double time_step = 0.5; // femtoseconds
    int num_steps = atoi(argv[2]);    // Number of steps
    MDSystem system(molecule, time_step, num_steps);

    // Initial force computation
    compute_forces(system, cndo);

    // Perform the MD simulation
    velocity_verlet(system, cndo);

    return 0;
}
