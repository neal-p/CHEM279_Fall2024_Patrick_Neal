from torch.utils.data import Dataset, random_split, DataLoader

from pathlib import Path
import pickle
import pandas as pd
import numpy as np
import torch
import math
from sklearn.metrics import pairwise_distances
from sklearn.model_selection import KFold
from sklearn.covariance import EllipticEnvelope
from sklearn.preprocessing import MinMaxScaler, StandardScaler
import warnings


def inverse_squared_distance_matrix(coords):

    dist_mat = pairwise_distances(coords)
    sq_dist_mat = dist_mat ** 2

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")

    feats = []
    for i in range(len(coords)):
        for j in range(i+1, len(coords)):

            feats.append(sq_dist_mat[i,j])
    
    return torch.tensor(feats).float()

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
            etot = row["normalized_etot"]

            feats = self.featurizer(coords)
            #Y = torch.tensor([etot] + list(grad.flatten())).float()
            Y = torch.tensor([etot]).float()

            self.data[idx] = (feats, Y)

        self.feat_dim = len(self.data[idx][0])
        self.Y_dim = len(self.data[idx][1])

    def __len__(self):
        return self.n

    def __getitem__(self, idx):
        return self.data[idx]


##################################################################################################################


def load_data(dir, K_folds=2, batch_size=100, featurizer=inverse_squared_distance_matrix, shuffle=True, random_state=None):

    dir = Path(dir)
    df = pd.read_csv(dir / "training_data.csv")

    # use a stupid simple outlier detector to 
    # remove structures that are super high energy
    # and not reasonable
    #df["outlier"] = EllipticEnvelope(contamination=0.1).fit_predict(df["etot"].values.reshape((-1, 1)))

    #print("outliers")
    #print("mean", df["etot"].mean())
    #print("median", df["etot"].median())
    #print(df[df["outlier"] == -1][["etot"]])
    #print("removing", (df["outlier"] == -1).sum(), "outliers...")

    #df = df[df["outlier"] == 1].reset_index(drop=True)

    #Normalize 
    energy_scaler = StandardScaler()
    energy_scaler.fit(df["etot"].values.reshape((-1, 1)))
    df["normalized_etot"] = energy_scaler.transform(df["etot"].values.reshape((-1, 1))).reshape(-1)

    with open("energy_scaler.pkl", "wb") as pkl:
        pickle.dump(energy_scaler, pkl)

    split_iter = KFold(n_splits=K_folds, shuffle=shuffle, random_state=random_state).split(df)
    folds = []
    print(f"Splitting into {K_folds} folds...", flush=True)

    for fold, (train_indices, test_indices) in enumerate(split_iter):

        print(f"    loading FOLD {fold}", flush=True)
        train_ds = MolDataset(df.iloc[train_indices].reset_index(drop=True), dir=dir, featurizer=featurizer, shuffle=shuffle, random_state=random_state)
        test_ds = MolDataset(df.iloc[test_indices].reset_index(drop=True), dir=dir, featurizer=featurizer, shuffle=shuffle, random_state=random_state)

        train_dl = DataLoader(train_ds, batch_size=batch_size, shuffle=shuffle)
        test_dl = DataLoader(test_ds, batch_size=batch_size, shuffle=shuffle)

        train_dl.feat_dim = train_ds.feat_dim
        test_dl.feat_dim = test_ds.feat_dim

        train_dl.Y_dim = train_ds.Y_dim
        test_dl.Y_dim = test_ds.Y_dim

        folds.append((train_dl, test_dl))

    print("Data loaded!")

    return folds, df, energy_scaler



