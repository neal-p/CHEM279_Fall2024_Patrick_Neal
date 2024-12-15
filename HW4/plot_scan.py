import pandas as pd
from matplotlib import pyplot as plt
import numpy as np



if __name__ == "__main__":

    distances = list(np.linspace(0, 10, 100) / 1.88973)
    df = pd.read_csv("scan_output.csv")

    df["bond_distance"] = df["file"].str.split("_").str[-1].str.strip(".xyz").astype(int).map(lambda x: distances[x])
    df["bond_distance"] = df["bond_distance"] * 0.529177
    df["energy"] = df["energy"] * 27.211
    df = df.sort_values("energy").reset_index()
    minimum = df.head(1)["bond_distance"].tolist()[0]
    energy = df.head(1)["energy"].tolist()[0]
    
    df.plot.scatter(x="bond_distance", y="energy")

    plt.axvline(x=1.098, color="red", linestyle="--", linewidth=1.5, label=f"Experimental Minimum energy at r={1.098:.3f}A")
    plt.scatter(minimum, energy, color="gold", s=200, marker="*", linewidth=1.5, label=f"CNDO/2 Minimum energy at r={minimum:.3f}A",)


    plt.xlabel("N2 Bond Distance (A)")
    plt.ylabel("Energy (eV)")

    plt.legend()
    plt.tight_layout()


    plt.savefig("scan_bond_distance.png")

    df.to_csv("scan_output.csv")
