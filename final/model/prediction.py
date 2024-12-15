import sys
from pathlib import Path
sys.path.insert(0, str(Path(".").resolve()))

import torch
from model import model
from model import featurize
import pickle
from sklearn.preprocessing import StandardScaler

class Predictor:

    def __init__(self, results_path):

        self.root = Path(results_path)

        with open(results_path / "energy_scaler.pkl", "rb") as pkl:
            self.energy_scaler = pickle.load(pkl)

        with open(results_path / "grad_scaler.pkl", "rb") as pkl:
            self.grad_scaler = pickle.load(pkl)

        with open(results_path / "sizes.pkl", "rb") as pkl:
            sizes = pickle.load(pkl)

        model_path = self.root / "best_model.pt"
        self.m = model.LinearModel(sizes)
        self.m.load_state_dict(torch.load(model_path, weights_only=True))

    def predict(self, coords):

        feats = featurize.internal_coords_water(coords)

        with torch.no_grad():
            y_hat = self.m(feats).detach().numpy()
            
        energy = y_hat[0]
        grad = y_hat[1:]

        energy = self.energy_scaler.inverse_transform([[energy]])[0][0]
        grad = self.grad_scaler.inverse_transform(grad.reshape((-1, 1))).reshape((3, 3))

        return energy, grad


