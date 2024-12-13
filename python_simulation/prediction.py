import torch
import model
import featurize


class EnergyPredictor:

    def __init__(self, model_path, energy_scaler, sizes):
        self.model_path = model_path
        self.m = model.LinearModel(sizes)
        self.m.load_state_dict(torch.load(model_path, weights_only=True))
        self.m.eval()
        self.energy_scaler = energy_scaler

    def predict(self, coords):

        feats = featurize.internal_coords_water(coords)

        with torch.no_grad():
            y_hat = self.m(feats)
            raw_energy = y_hat[0].item()

        energy = self.energy_scaler.inverse_transform([[raw_energy]])[0][0]
        return energy


class GradPredictor:

    def __init__(self, model_path, sizes):
        self.model_path = model_path

        self.m = model.LinearModel(sizes)
        self.m.load_state_dict(torch.load(model_path, weights_only=True))
        self.m.eval()

    def predict(self, coords):

        feats = featurize.internal_coords_water(coords)

        with torch.no_grad():
            y_hat = self.m(feats).detach().numpy()

        grad = y_hat.reshape((3, 3))
        return grad


