# Machine Learning Accelerated Dynamics

Final project for Chem279. We extrapolated from the CNDO/2 code we created through the course of this class to design a molecular dynamics simulation of a simple water molecule. Additionally, we designed an experiment to evaluate the feasibility and computational performance of a machine learning accelearted MD simulation, implemented in Python, where key gradient and energy calculations from the CNDO/2 simulation are replaced with predictions from a neural network. This neural network is trained on gradients and energies originally calculated from the CNDO/2 method. For more information, refer to the project report.

## C++ CNDO/2 Simulation Compilation Instructions

Code for the C++ CNDO/2 simulation is stored in the CNDO2_simulation/c++ directory. Run the following commands to run the C++ CNDO/2 MD simulation with 1000 iterations:

```
cd c++
make
```

## Python Simulation Compilation Instructions

Run the following commands to build the Python simulation Docker container, generate training data, train the model, and run simulations:

```
./build_container.sh
./generate_training_data.sh
./train_model.sh
./run_ML_simulation.sh
```