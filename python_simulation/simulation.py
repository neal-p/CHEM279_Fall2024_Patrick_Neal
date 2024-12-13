from utils import utils
import prediction
from sklearn.preprocessing import StandardScaler
import numpy as np
import pickle
import time


def simulation(coords, num_steps, e_predictor, grad_predictor):

    #all_coords = [coords.copy()]
    all_energies = []

    energy = e_predictor.predict(coords)
    all_energies.append(energy)

    forces = -1 * grad_predictor.predict(coords)
    forces_norm = np.linalg.norm(forces)

    for step in range(10):

        forces = -1 * grad_predictor.predict(coords)

        randomness = np.random.normal(0, 1, (3, 3))
        randomness_norm = np.linalg.norm(randomness)
        scaler = 0.5 * forces_norm / randomness_norm
        randomness *= scaler
        
        H1 = coords[1]
        H2 = coords[2]
        O = coords[0]
        
        bd1 = np.linalg.norm((O - H1))
        bd2 = np.linalg.norm((O - H2))

        if bd1 > 2.0:
            vec = (H1 - O) / bd1  # Unit vector along bond
            correction = (bd1 - 2.0) * vec  # Adjust distance to 2.0
            coords[1] -= correction / 2  # Move H1 slightly
            coords[0] += correction / 2  # Move O slightly to conserve COM

        if bd2 > 2.0:
            vec = (H2 - O) / bd2  # Unit vector along bond
            correction = (bd2 - 2.0) * vec  # Adjust distance to 2.0
            coords[2] -= correction / 2  # Move H2 slightly
            coords[0] += correction / 2  # Move O slightly to conserve COM


        coords = coords + 0.01 * (forces + randomness)
        
            
        energy = e_predictor.predict(coords)
        #all_coords.append(coords.copy())
        all_energies.append(energy)

    return #all_coords, all_energies





if __name__ == "__main__":


    molecule_file = "H2O.txt"
    elements, coords = utils.read_class_xyz(molecule_file)
    masses = np.array([utils.Element(e).getMass() for e in elements])  # in AMU

    time_step = 0.5 # femtoseconds
    num_steps = 100 # number of steps

    e_sizes = [
                (3, 100),

                (100, 100),

                (100, 1)
            ]

    grad_sizes = [
                (3, 100),

                (100, 100),

                (100, 9)
            ]

    with open("energy_scaler.pkl", "rb") as pkl:
        energy_scaler = pickle.load(pkl)

    e_predictor = prediction.EnergyPredictor("energy_fold_0_best_model.pt", energy_scaler, e_sizes)
    grad_predictor = prediction.GradPredictor("gradient_fold_2_best_model.pt", grad_sizes)

    start = time.time()
    simulation(coords, num_steps, e_predictor, grad_predictor)
    end = time.time()

    # convert to angstrom for visualization
    # traj = [c / 1.8897 for c in traj]
    # utils.write_multi_xyz(elements, traj, f"results.xyz")


    print(end - start)


