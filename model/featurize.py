from torch.utils.data import Dataset, random_split

from pathlib import Path
import pandas as pd
import numpy as np
import torch
import math
from sklearn.metrics import pairwise_distances
from sklearn.model_selection import KFold
import warnings


def inverse_squared_distance_matrix(coords):

    dist_mat = pairwise_distances(coords)
    sq_dist_mat = dist_mat ** 2

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        inverse = 1 / sq_dist_mat

    inverse = np.where(np.isnan(inverse), 0, inverse)

    return torch.tensor(inverse).flatten().float()

class MolDataset(Dataset):
    def __init__(self, df, dir, featurizer=inverse_squared_distance_matrix, shuffle=True, random_state=None):

        dir = Path(dir)
        self.featurizer = featurizer
        self.n = len(df)
        self.data = {}

        if shuffle:
            df.sample(frac=1, random_state=random_state)

        for idx, row in df.iterrows():
            
            coords = np.loadtxt(dir / row["coord_file"])
            grad = np.loadtxt(dir / row["grad_file"])
            etot = row["etot"]

            feats = self.featurizer(coords)
            Y = torch.tensor([etot] + list(grad.flatten())).float()

            self.data[idx] = (feats, Y)

        self.feat_dim = len(self.data[idx][0])
        self.Y_dim = len(self.data[idx][1])

    def __len__(self):
        return self.n - 1

    def __getitem__(self, idx):

        if not idx in self.data:
            print("cant find", idx)
            print(list(self.data.keys()), flush=True)
        return self.data[idx]


##################################################################################################################


def load_data(dir, K_folds=5, featurizer=inverse_squared_distance_matrix, shuffle=True, random_state=None):

    dir = Path(dir)
    df = pd.read_csv(dir / "training_data.csv")

    split_iter = KFold(n_splits=K_folds, shuffle=shuffle, random_state=random_state).split(df)
    folds = []
    print(f"Splitting into {K_folds} folds...", flush=True)

    for fold, (train_indices, test_indices) in enumerate(split_iter):

        print(f"    loading FOLD {fold}", flush=True)
        train_ds = MolDataset(df.iloc[train_indices].reset_index(drop=True), dir=dir, featurizer=featurizer, shuffle=shuffle, random_state=random_state)
        test_ds = MolDataset(df.iloc[test_indices].reset_index(drop=True), dir=dir, featurizer=featurizer, shuffle=shuffle, random_state=random_state)

        folds.append((train_ds, test_ds))

    print("Data loaded!")

    return folds



