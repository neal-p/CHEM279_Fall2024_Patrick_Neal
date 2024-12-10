import torch
import featurize
import model
from matplotlib import pyplot as plt
import numpy as np
from sklearn.metrics import root_mean_squared_error, r2_score
from tqdm import tqdm


def train_step(m, loss_fn, opt, dl):

    m.train(True)
    total_loss = 0.0
    count = 0

    for feats, Y in dl:
        opt.zero_grad()
        count += len(feats)

        y_hat = m(feats)

        loss = loss_fn(y_hat, Y)
        loss.backward()
        opt.step()
        total_loss += loss.item()

    return total_loss / count


def test_step(m, loss_fn, dl):

    m.eval()
    total_loss = 0.0
    count = 0

    with torch.no_grad():
        for feats, Y in dl:

            count += len(feats)

            y_hat = m(feats)

            loss = loss_fn(y_hat, Y)
            total_loss += loss.item()

    return total_loss / count

def energy_corelation(m, epoch, dl, energy_scaler, plot_file_stem):

    m.eval()

    true_energies = []
    pred_energies = []

    with torch.no_grad():
        for feats, Y in dl:

            for pred_row, Y_row in zip(m(feats), Y):
               
               pred_energies.append(pred_row[0].item())
               true_energies.append(Y_row[0].item())

    true_energies = np.array(true_energies).reshape((-1, 1))
    pred_energies = np.array(pred_energies).reshape((-1, 1))

    rescaled_true_energies = energy_scaler.inverse_transform(true_energies).reshape(-1)
    rescaled_pred_energies = energy_scaler.inverse_transform(pred_energies).reshape(-1)

    rmse = root_mean_squared_error(rescaled_true_energies, rescaled_pred_energies)
    r2 = r2_score(rescaled_true_energies, rescaled_pred_energies)

    plt.clf()
    plt.scatter(rescaled_true_energies, rescaled_pred_energies)
    plt.title(f"True vs Predicted Energy: RMSE={rmse:.3f}, R2={r2:.3f} (epoch={epoch})")
    plt.xlabel("True Energy")
    plt.ylabel("Predicted Energy")

    #plt.show()
    plt.savefig(f"{plot_file_stem}_energy_correlation.png")

    print(f"        RMSE={rmse:.3f}, R2={r2:.3f}", flush=True)


def training_loop(m, fold, loss_fn, train_dl, test_dl, epochs=1000, lr=1e-4, plot_file_stem=""):

    opt = torch.optim.Adam(m.parameters(), lr=lr)

    train_losses = []
    test_losses = []
    best_loss = 10000000

    train_loss = test_step(m, loss_fn, train_dl)
    test_loss = test_step(m, loss_fn, test_dl)

    train_losses.append(train_loss)
    test_losses.append(test_loss)

    for epoch in range(epochs):
        print(epoch, flush=True)

        # Training Step
        train_step(m, loss_fn, opt, train_dl)

        # Eval losses
        train_loss = test_step(m, loss_fn, train_dl)
        test_loss = test_step(m, loss_fn, test_dl)

        print(f"    EPOCH{epoch}: Train loss={train_loss:.5f}, Test loss={test_loss:.5f}", flush=True)

        train_losses.append(train_loss)
        test_losses.append(test_loss)

        if test_loss < best_loss:
            best_loss = test_loss
            print("        lowered test loss!")
            torch.save(m.state_dict(), f"best_model_{fold}.pt")

        plt.clf()
        plt.plot(list(range(len(train_losses))), train_losses, label="train", c="r")
        plt.plot(list(range(len(test_losses))), test_losses, label="test", c="g")
        plt.legend()
        plt.title(f"Loss Curve (epoch={epoch})")
        plt.xlabel("Epoch")
        plt.ylabel("Loss")

        #plt.show()
        plt.savefig(f"{plot_file_stem}_loss_curve.png")

        if epoch % 10 == 0:
            energy_corelation(m, epoch, test_dl, energy_scaler, plot_file_stem="rework_in_progress")


if __name__ == "__main__":

    print("START", flush=True)

    folds, feat_dim, Y_dim, energy_scaler = featurize.load_data("data.pkl")

    sizes = [
                (feat_dim, 100),

                (100, 200),
                (200, 200),
                (200, 100),

                (100, Y_dim)
            ]

    print("SIZES:")
    print(sizes, flush=True)

    for fold, (train_dl, test_dl) in enumerate(folds):

        print(f"FOLD {fold} training...", flush=True)

        m = model.LinearModel(sizes)
        loss_fn = torch.nn.MSELoss()

        training_loop(m, fold, loss_fn, train_dl, test_dl, plot_file_stem=f"rework_{fold}")

        m.load_state_dict(torch.load(f"best_model_{fold}.pt"))
        energy_corelation(m, test_dl, energy_scaler, plot_file_stem=f"rework_final_{fold}")


