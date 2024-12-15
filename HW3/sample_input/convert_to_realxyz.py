import sys

BOHR_TO_ANGSTROM = 0.529177

if __name__ == "__main__":

    if len(sys.argv) < 2:
        print("HEY GIVE ME A FILE TO CONVERT")
        exit()

    with open(sys.argv[1], "r") as f:
        lines = f.readlines()

    output_file = sys.argv[1].split(".")[0] + ".xyz"
    output = ""

    n_atoms, charge = lines[0].strip().split()
    output += f"{n_atoms}\n"
    output += f"from {sys.argv[1]} with charge {charge}\n"

    for line in lines[1:]:
        element, x, y, z = line.strip().split()
        x, y, z = float(x), float(y), float(z)

        x *= BOHR_TO_ANGSTROM
        y *= BOHR_TO_ANGSTROM
        z *= BOHR_TO_ANGSTROM

        output += f"{element} {x} {y} {z}\n"

    with open(output_file, "w") as f:
        f.write(output)

    print("DONE")







