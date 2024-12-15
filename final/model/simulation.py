import sys
from pathlib import Path
sys.path.insert(0, str((Path(".") / "utils").parent.resolve()))
sys.path.insert(0, str((Path(".") / "model").parent.resolve()))

from utils import parsers 
from model import prediction
import numpy as np
import pickle
import time


def simulation(coords, num_steps, predictor):

    all_coords = [coords.copy()]
    all_energies = []

    energy, grad = predictor.predict(coords)
    all_energies.append(energy)

    forces = -1 * grad 
    forces_norm = np.linalg.norm(forces)

    for step in range(5000):

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

        energy, grad = predictor.predict(coords)
        forces = -1 * grad

        all_energies.append(energy)
        all_coords.append(coords.copy())

    return all_coords, all_energies, [], []





if __name__ == "__main__":

    if len(sys.argv) < 2:
        print("ERROR")
        print("Provide starting structure")
        exit()

    trained_model_dir = Path("model_dir")
    initial_geometry = Path(sys.argv[1])

    elements, coords = parsers.read_xyz(initial_geometry)
    coords *= 1.8897 # convert to Bohr
    masses = np.array([parsers.Element(e).getMass() for e in elements])  # in AMU

    time_step = 0.5 # femtoseconds
    num_steps = 100 # number of steps

    predictor = prediction.Predictor(trained_model_dir)

    start = time.time()
    traj, energies, distances, angles = simulation(coords, num_steps, predictor)
    end = time.time()

    traj = [c / 1.8897 for c in traj]  # convert back to angstrom
    parsers.write_multi_xyz(elements, traj, f"traj.xyz")


    print(end - start)


