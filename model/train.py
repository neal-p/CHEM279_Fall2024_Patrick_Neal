import torch
import featurize
import model





if __name__ == "__main__":

    folds, feat_dim, Y_dim, energy_scaler = featurize.load_data("data.pkl")

    sizes = [
                (feat_dim, 100),

                (100, 200),
                (200, 200),
                (200, 100),

                (100, Y_dim)
            ]

    for fold, (train_dl, test_dl) in enumerate(folds):

        m = model.LinearModel(sizes)
        loss_fn = torch.nn.MSELoss()


