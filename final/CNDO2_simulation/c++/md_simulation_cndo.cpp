#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>
#include <armadillo>
#include "CNDO.h"
#include "AO.h"
#include "util.h"
#include <chrono>

struct MDSystem
{
    Molecule_basis &mol;
    arma::mat velocities;                       // Velocities of atoms (3 x num_atoms)
    arma::mat forces;                           // Forces on atoms (3 x num_atoms)
    double time_step;                           // Time step for integration
    int num_steps;                              // Number of integration steps
    std::vector<double> energies;               // Energies at each step
    std::vector<arma::mat> coordinates_history; // Coordinates at each step

    MDSystem(Molecule_basis &mol_i, int num_atoms, double dt, int steps) : mol(mol_i), time_step(dt), num_steps(steps)
    {
        velocities = arma::zeros(3, num_atoms);
        forces = arma::zeros(3, num_atoms);
    }

    void log_state(double energy, const arma::mat &positions)
    {
        energies.push_back(energy);
        coordinates_history.push_back(positions);
    }

    void write_log_to_file(const std::string &filename)
    {
        std::ofstream outfile(filename);

        if (!outfile.is_open())
        {
            std::cerr << "Error opening file: " << filename << std::endl;
            return;
        }

        // Write header
        outfile << "Step,Energy,Coordinates\n";

        for (size_t step = 0; step < energies.size(); ++step)
        {
            outfile << step << "," << std::fixed << std::setprecision(6) << energies[step] << ",";

            // Flatten and write the coordinates as a single Numpy-parsable array
            const arma::mat &coords = coordinates_history[step];
            for (size_t i = 0; i < coords.n_elem; ++i)
            {
                outfile << coords(i);
                if (i < coords.n_elem - 1)
                    outfile << ",";
            }
            outfile << "\n";
        }

        outfile.close();
    }

    void compute_forces(CNDO cndo)
    {
        // Compute the forces using CNDO gradient
        arma::mat gradient = cndo.getGradient();
        forces = -gradient; // Force is the negative gradient of the energy
    }

    void simulation()
    {
        CNDO cndo(mol, 100, 1e-6);
        compute_forces(cndo);

        int num_atoms = mol.getnum_atoms();
        arma::mat positions(3, num_atoms);

        // Initialize positions
        for (size_t i = 0; i < num_atoms; ++i)
        {
            positions.col(i) = mol.mAtoms[i].mAOs[0].get_R0();
        }

        std::vector<double> energies;
        energies.reserve(num_steps); // Optional optimization

        for (int step = 0; step < num_steps; ++step)
        {
            std::cout << "Step " << step << ":" << std::endl;

            double forces_norm = arma::norm(forces);

            for (size_t i = 0; i < num_atoms; ++i)
            {
                arma::vec randomness = arma::randn<arma::vec>(3);
                double randomness_norm = arma::norm(randomness);

                double scaler = 0.5 * forces_norm / randomness_norm;
                randomness *= scaler;

                arma::vec H1 = mol.mAtoms[1].mAOs[0].get_R0();
                arma::vec H2 = mol.mAtoms[2].mAOs[0].get_R0();
                arma::vec O = mol.mAtoms[0].mAOs[0].get_R0();

                // Bond length corrections
                arma::vec H1_O_vec = H1 - O;
                arma::vec H2_O_vec = H2 - O;

                double bd1 = arma::norm(H1_O_vec);
                double bd2 = arma::norm(H2_O_vec);

                if (bd1 > 2.0)
                {
                    arma::vec correction = (bd1 - 2.0) * (H1_O_vec / bd1);
                    arma::vec updated_H1 = H1 - correction / 2.0;
                    arma::vec updated_O = O + correction / 2.0;
                    mol.mAtoms[1].mAOs[0].set_R0(updated_H1);
                    mol.mAtoms[0].mAOs[0].set_R0(updated_O);
                }

                if (bd2 > 2.0)
                {
                    arma::vec correction = (bd2 - 2.0) * (H2_O_vec / bd2);
                    arma::vec updated_H2 = H2 - correction / 2.0;
                    arma::vec updated_O = O + correction / 2.0;
                    mol.mAtoms[2].mAOs[0].set_R0(updated_H2);
                    mol.mAtoms[0].mAOs[0].set_R0(updated_O);
                }

                // Update position using forces and randomness
                arma::vec position = positions.col(i);
                arma::vec force = forces.col(i);

                position += 0.01 * (force + randomness);
                positions.col(i) = position;
            }

            // Update molecule positions in CNDO
            for (size_t i = 0; i < num_atoms; ++i)
            {
                arma::vec temp_position = positions.col(i);
                mol.mAtoms[i].set_R0(temp_position);
            }

            CNDO cndo(mol, 100, 1e-6);

            cndo.init();
            cndo.run();

            compute_forces(cndo);

            velocities += 0.5 * (forces)*time_step;

            double energy = cndo.getEnergy();
            energies.push_back(energy);
            log_state(energy, positions);

            std::cout << "Positions:\n"
                      << positions << std::endl;
            std::cout << "Velocities:\n"
                      << velocities << std::endl;

            // Check for stagnation
            if (energies.size() >= 3)
            {
                double e1 = energies[energies.size() - 1];
                double e2 = energies[energies.size() - 2];

                if (std::abs(e1 - e2) < 1e-6 < 1e-6)
                {
                    for (size_t i = 0; i < num_atoms; ++i)
                    {
                        arma::vec random_displacement = arma::randn<arma::vec>(3) * 0.1; // Small random displacement
                        arma::vec new_position = positions.col(i) + random_displacement;
                        mol.mAtoms[i].set_R0(new_position);
                        positions.col(i) = new_position;
                    }

                    // Reinitialize CNDO with randomized positions
                    CNDO cndo = CNDO(mol, 100, 1e-6);
                    cndo.init();
                    cndo.run();
                    compute_forces(cndo);
                }
            }
        }
    }
};

int main(int argc, char *argv[])
{
    // Define input files for molecule initialization and basis sets
    std::string molecule_file = argv[1];
    arma::mat H_basis, C_basis;

    // Load STO-3G basis sets for hydrogen and carbon (example)
    H_basis.load("H_STO3G.txt");
    C_basis.load("C_STO3G.txt");

    // Initialize the molecule
    Molecule_basis molecule(molecule_file, H_basis, C_basis);

    auto start = std::chrono::high_resolution_clock::now();

    // // Set up the CNDO/2 calculation
    // CNDO cndo(molecule, 100, 1e-6); // Max iterations and tolerance

    // if (cndo.init() != 0)
    // {
    //     std::cerr << "Failed to initialize CNDO." << std::endl;
    //     return 1;
    // }

    // if (cndo.run() != 0)
    // {
    //     std::cerr << "SCF calculation did not converge." << std::endl;
    //     return 1;
    // }

    // Set up the MD system
    double time_step = 0.5;        // femtoseconds
    int num_steps = atoi(argv[2]); // Number of steps
    MDSystem system(molecule, molecule.getnum_atoms(), time_step, num_steps);

    // Initial force computation
    // compute_forces(system, cndo);

    // Perform the MD simulation
    system.simulation();

    // Write results to an output file
    std::string output_file = "md_simulation_output.csv";
    system.write_log_to_file(output_file);

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Time taken for " << num_steps << " steps: " << elapsed.count() << " seconds." << std::endl;
    std::cout << "Results written to " << output_file << std::endl;

    return 0;
}
