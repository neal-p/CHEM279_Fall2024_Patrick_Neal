from pathlib import Path
from utils import utils
import pandas as pd


if __name__ == "__main__":

    df = []

    for file in Path(".").glob("s*out"):
        print("Reading", file, "...")

        try:
            energy = utils.read_energy(file)
            grad = utils.read_grad(file)
            elements, coords = utils.read_xyz(str(file).replace(".out", ".xyz"))

            df.append((file, coords, energy, grad))
        except:
            pass

    df = pd.DataFrame(df, columns=["file", "coords", "energy", "grad"])
    df.to_pickle("data.pkl")

    print(df)

