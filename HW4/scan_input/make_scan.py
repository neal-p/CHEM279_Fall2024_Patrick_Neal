import numpy as np


if __name__ == "__main__":

    #                                                   vvv

    distances = list(np.linspace(0, 10, 100) / 1.88973)
 

    for idx, d in enumerate(distances):
        with open(f"N2_{idx}.xyz", "w") as f:
            f.write("2 0\n")
            f.write("7 0 0 0\n")
            f.write(f"7 {d} 0 0\n")
