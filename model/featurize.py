import numpy as np
import torch
import pandas as pd
from pathlib import Path
from torch.utils.data import Dataset, DataLoader
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import KFold
import pickle


def angle(atom1, atom2, atom3):

    v1 = atom1 - atom2
    v2 = atom3 - atom2

    dp = np.dot(v1, v2)
    m1 = np.linalg.norm(v1)
    m2 = np.linalg.norm(v1)

    cos_theta = np.clip(dp / (m1 * m2), -1, 1)
    return np.arccos(cos_theta)


def dihedral(atom1, atom2, atom3, atom4):
    v_12 = atom1 - atom2
    v_23 = atom3 - atom2
    v_34 = atom4 - atom3

    n1 = np.cross(v_12, v_23)
    n2 = np.cross(v_23, v_34)

    n1 /= np.linalg.norm(n1)
    n2 /= np.linalg.norm(n2)

    m = np.cross(n1, v_23 / np.linalg.norm(v_23))

    x = np.dot(n1, n2)
    y = np.dot(m, n2)

    return np.arctan2(y, x)


def distance(atom1, atom2):
    return np.linalg.norm(atom1 - atom2)


def internal_coords(coords):

    F1 = coords[0]
    C1 = coords[1]
    C2 = coords[2]
    F2 = coords[3]
    H1 = coords[4]
    H2 = coords[5]
    H3 = coords[6]
    H4 = coords[7]

    feats = [

            1 / distance(C1, C2),
            1 / distance(C1, F1),
            angle(C2, C1, F1),

            1 / distance(C2, F2),
            angle(C1, C2, F2),

            dihedral(F1, C1, C2, F2),
            dihedral(H1, C1, C2, H3),
            dihedral(H2, C1, C2, H4),
            dihedral(H1, C1, C2, H4),

            1 / distance(F1, F2),

            1 / distance(C1, H1),
            1 / distance(C1, H2),
            angle(C2, C1, H1),
            angle(C2, C1, H2),
            angle(H1, C1, H2),

            1 / distance(C2, H3),
            1 / distance(C2, H4),
            angle(C1, C2, H3),
            angle(C1, C2, H4),
            angle(H3, C2, H4),
 
            ]

    return torch.tensor(feats).float()


class MolDataset(Dataset):

    def __init__(self, df):

        df = df.reset_index()

        self.n = len(df)
        self.data = {}

        for idx, row in df.iterrows():

            coords = row["coords"]
            energy = row["normalized_energy"]
            grad = row["grad"]

            feats = internal_coords(coords)
            Y = torch.tensor([energy])

            self.data[idx] = (feats, Y)

        self.feat_dim = len(feats)
        self.Y_dim = len(Y)

    def __len__(self):
        return self.n

    def __getitem__(self, idx):
        return self.data[idx]


##########################################################

def load_data(file, folds=2, batch_size=128):

    print("Reading df...", flush=True)
    df = pd.read_pickle(file)

    print(df, flush=True)

    # Normalize energies
    energy_scaler = StandardScaler()
    df["normalized_energy"] = energy_scaler.fit_transform(df["energy"].values.reshape((-1, 1))).reshape(-1)

    with open("energy_scaler.pkl", "wb") as pkl:
        pickle.dump(energy_scaler, pkl)

    #splitter = KFold(n_splits=folds, shuffle=True)
    #splits = splitter.split(df)
    df = df.sample(frac=1).reset_index(drop=True)

    train_indices = list(df.head(8500).index)
    test_indices = list(df.tail(1000).index)
    splits = [(train_indices, test_indices)]

    folds = []
    for fold, (train_indices, test_indices) in enumerate(splits):

        print("Loading fold", fold, "...", flush=True)

        fold_train_df = df.iloc[train_indices]
        fold_test_df = df.iloc[test_indices]

        train_ds = MolDataset(fold_train_df)
        test_ds = MolDataset(fold_test_df)

        train_dl = DataLoader(train_ds, batch_size=batch_size, shuffle=True)
        test_dl = DataLoader(test_ds, batch_size=batch_size, shuffle=True)

        folds.append((train_dl, test_dl))

    return folds, train_ds.feat_dim, train_ds.Y_dim, energy_scaler



