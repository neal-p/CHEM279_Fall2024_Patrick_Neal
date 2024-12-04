import numpy as np
import pandas as pd
from pathlib import Path


SENTINEL = -10000

def get_gradient(file):

    with open(file, "r") as f:
        lines = f.readlines()

    grad_section = 0
    for idx, line in enumerate(lines):
        line = line.strip()
        if line == "FINAL  POINT  AND  DERIVATIVES":
            grad_section = idx

    if grad_section == 0:
        raise ValueError(f"Could not read gradient from {file}")


    idx = grad_section + 3
    info = []

    while "KCAL/ANGSTROM" in lines[idx]:
        #print(lines[idx])
        _, atom_idx, element, _, dim, _, grad, _ = lines[idx].strip().split()

        info.append((int(atom_idx), dim, float(grad)))
        idx += 1

    natoms = max((int(i[0]) for i in info))
    grad_matrix = np.zeros((natoms, 3)) + SENTINEL

    for atom_idx, dim, grad in info:

        if dim == "X":
            col = 0
        elif dim == "Y":
            col = 1
        elif dim == "Z":
            col = 2
        else:
            raise ValueError(f"unknown dim {dim}")


        grad_matrix[atom_idx-1, col] = grad

    assert (grad_matrix == SENTINEL).sum() == 0, f"grad not complete"
    return grad_matrix


def get_atomic_coordinates(file):

    with open(file, "r") as f:
        lines = f.readlines()

    coord_section = 0
    for idx, line in enumerate(lines):
        line = line.strip()
        if line == "CARTESIAN COORDINATES":
            coord_section = idx

    if coord_section == 0:
        raise ValueError(f"Could not read atomic coordinates from {file}")

    idx = coord_section + 2

    info = []
    
    while lines[idx].strip() != "":
        atom_idx, element, x, y, z = lines[idx].strip().split()
        if atom_idx == "NO.":
            idx += 2
            continue

        info.append((atom_idx, element, x, y, z))
        idx += 1

    natoms = max((int(i[0]) for i in info))
    coords = np.zeros((natoms, 3)) + SENTINEL
    elements = [None for _ in range(natoms)]

    for atom_idx, element, x, y, z in info:
        coords[int(atom_idx)-1, :] = (x ,y ,z)
        elements[int(atom_idx)-1] = element

    assert (coords == SENTINEL).sum() == 0, f"cords cnot compete"
    assert all((e is not None for e in elements)), f"cords not complete"

    return elements, coords


def get_energies(file):

    with open(file, "r") as f:
        lines = f.readlines()

    energy_section = 0
    for idx, line in enumerate(lines):
        line = line.strip()
        if line == "***  SUMMARY OF ENERGY PARTITION  ***":
            energy_section = idx

    if energy_section == 0:
        raise ValueError(f"Could not read energy decompostion from {file}")

    idx = energy_section + 2
    line = lines[idx].strip()

    nuclear_nuclear_repulsion = SENTINEL
    etot = SENTINEL
    
    while not line.startswith("ETOT (EONE + ETWO"):

        if line.startswith("NUCLEAR-NUCLEAR REPULSION"):
            _, _, nuclear_nuclear_repulsion, _ = line.split()
            nuclear_nuclear_repulsion = float(nuclear_nuclear_repulsion)

        idx += 1
        line = lines[idx].strip()

    if line.startswith("ETOT (EONE + ETWO"):
        _, _, _, _, etot, _ = line.split()
        etot = float(etot)

    assert nuclear_nuclear_repulsion != SENTINEL, "could not read nuclear nuclear repulsion"
    assert etot != SENTINEL, "coult not read etot"

    return etot, nuclear_nuclear_repulsion, etot-nuclear_nuclear_repulsion



if __name__ == "__main__":

    df = []

    for file in Path(".").glob("sample*.out"):
        print("reading", file, "...")
        

        elements, coords = get_atomic_coordinates(file)
        etot, electronic, nuclear = get_energies(file)
        grad = get_gradient(file)
  
        grad_file = file.stem + "_grad.np"
        np.savetxt(grad_file, grad)

        coord_file = file.stem + "_coords.np"
        np.savetxt(coord_file, coords)

        df.append((file, etot, electronic, nuclear, grad_file, coord_file))

    df = pd.DataFrame(df, columns=["mopac_output_file", "etot", "electronic", "nuclear", "grad_file", "coord_file"])
    df.to_csv("training_data.csv", index=False)
    print(df)


