import sys
from pathlib import Path
sys.path.insert(0, str(Path(".").resolve()))

from utils import parsers
import pandas as pd


if __name__ == "__main__":

    df = []

    for file in Path(".").glob("s*out"):
        print("Reading", file, "...")

        energy, grad = parsers.read_our_output(file)
        elements, coords = parsers.read_class_xyz(str(file).replace(".out", ".xyz"))

        df.append((file, coords, energy, grad))

    df = pd.DataFrame(df, columns=["file", "coords", "energy", "grad"])
    df.to_pickle("data.pkl")

    print(df)

