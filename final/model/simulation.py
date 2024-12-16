import sys
from pathlib import Path
sys.path.insert(0, str((Path(".") / "utils").parent.resolve()))
sys.path.insert(0, str((Path(".") / "model").parent.resolve()))

from utils import parsers 
from model import prediction, featurize
import numpy as np
import pickle
import time
from matplotlib import pyplot as plt
import pandas as pd


print(str(Path(".").resolve()))
print(list(Path(".").glob("*")))
print(list(Path("results").glob("*")))


def simulation(coords, num_steps, predictor):

    all_coords = [coords.copy()]
    all_energies = []
    h1_distances = []
    h2_distances = []
    angles = []

    energy, grad = predictor.predict(coords)
    all_energies.append(energy)

    forces = -1 * grad 
    forces_norm = np.linalg.norm(forces)

    for step in range(500):

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
        h1_distances.append(bd1)
        h2_distances.append(bd2)
        angles.append(np.degrees(featurize.angle(H1, O, H2)))

    return all_coords, all_energies, h1_distances, h2_distances, angles





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
    traj, energies, h1_distances, h2_distances, angles = simulation(coords, num_steps, predictor)
    end = time.time()

    traj = [c / 1.8897 for c in traj]  # convert back to angstrom
    parsers.write_multi_xyz(elements, traj, f"results/traj.xyz")

    plt.plot(list(range(len(energies))), energies)
    plt.title("Simulation Energy")
    plt.xlabel("Step")
    plt.ylabel("Energy (eV)")
    plt.savefig(f"results/traj_energy.png")
    plt.clf()

    plt.plot(list(range(len(h1_distances))), h1_distances, label="H1-O")
    plt.plot(list(range(len(h2_distances))), h2_distances, label="H2-O")
    plt.title("Simulation Distances")
    plt.xlabel("Step")
    plt.ylabel("H-O Distance (Bohr)")
    plt.legend()
    plt.savefig(f"results/traj_distances.png")
    plt.clf()

    plt.plot(list(range(len(angles))), angles)
    plt.title("Simulation Angles")
    plt.xlabel("Step")
    plt.ylabel("Angles (degrees)")
    plt.savefig(f"results/traj_angles.png")
    plt.clf()

    
    df = pd.DataFrame(list(zip(energies, h1_distances, h2_distances, angles)), columns=["energy", "H1-O_distance", "H2-O_distance", "H1-O-H2_angle"])
    df.to_csv("traj_info.csv")

    print(end - start)


