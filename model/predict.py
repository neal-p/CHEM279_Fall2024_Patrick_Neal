import torch
import model
import featurize
import pickle
from sklearn.preprocessing import StandardScaler
from torch.utils.data import DataLoader
import pandas as pd
import numpy as np


if __name__ == "__main__":

    dir = "prepare_training_data"
    df_path = "prepare_training_data/training_data.csv"

    df = pd.read_csv(df_path)
    df["normalized_etot"] = 1
    ds = featurize.MolDataset(df, dir, shuffle=False)
    dl = DataLoader(ds, batch_size=100, shuffle=False)

    sizes = [
                (ds.feat_dim, 100),
                
                (100, 500),
                (500, 500),
                (500, 500),
                (500, 500),
                (500, 100),

                (100, ds.Y_dim)
            ]




    m = model.LinearModel(sizes=sizes)
    m.load_state_dict(torch.load("best_model.pt"))
    m.eval()

    with open("energy_scaler.pkl", "rb") as pkl:
        energy_scaler = pickle.load(pkl)

    predicted_energies = []

    with torch.no_grad():
        for feats, Y in dl:
            y_hat = m(feats)
            predicted_energies.extend(y_hat.detach().tolist())

    predicted_energies = np.array(predicted_energies).reshape((-1, 1))
    predicted_energies = energy_scaler.inverse_transform(predicted_energies)
    predicted_energies = predicted_energies.reshape(-1)

    df["predicted_etot"] = predicted_energies

    df.to_csv("predictions.csv", index=False)

