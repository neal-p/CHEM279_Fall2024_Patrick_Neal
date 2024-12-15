import sys
from pathlib import Path
sys.path.insert(0, str(Path(".").resolve()))

import torch
import featurize
import model
from matplotlib import pyplot as plt
import numpy as np
from sklearn.metrics import root_mean_squared_error, r2_score
from tqdm import tqdm
import pickle
import shutil


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

    return rmse, r2


def grad_corelation(m, epoch, dl, grad_scaler, plot_file_stem):

    m.eval()

    true_grads = []
    pred_grads = []

    with torch.no_grad():
        for feats, Y in dl:

            for pred_row, Y_row in zip(m(feats), Y):
               
                pred_grads.extend(pred_row[1:].detach().tolist())
                true_grads.extend(Y_row[1:].detach().tolist())

    true_grads = np.array(true_grads).reshape((-1, 1))
    pred_grads = np.array(pred_grads).reshape((-1, 1))

    true_grads = grad_scaler.inverse_transform(true_grads).reshape(-1)
    pred_grads = grad_scaler.inverse_transform(pred_grads).reshape(-1)

    rmse = root_mean_squared_error(true_grads, pred_grads)
    r2 = r2_score(true_grads, pred_grads)

    plt.clf()
    plt.scatter(true_grads, pred_grads)
    plt.title(f"True vs Predicted Gradient Components: RMSE={rmse:.3f}, R2={r2:.3f} (epoch={epoch})")
    plt.xlabel("True Grad Component")
    plt.ylabel("Predicted Grad Component")

    #plt.show()
    plt.savefig(f"{plot_file_stem}_grad_correlation.png")

    print(f"        RMSE={rmse:.3f}, R2={r2:.3f}", flush=True)

    return rmse, r2

def training_loop(m, fold, loss_fn, train_dl, test_dl, plot_file_stem, model_file, epochs=100, lr=1e-4):

    opt = torch.optim.Adam(m.parameters(), lr=lr)

    train_losses = []
    test_losses = []
    best_loss = 10000000

    train_loss = test_step(m, loss_fn, train_dl)
    test_loss = test_step(m, loss_fn, test_dl)

    train_losses.append(train_loss)
    test_losses.append(test_loss)

    for epoch in range(epochs):
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
            torch.save(m.state_dict(), model_file)

        plt.clf()
        plt.plot(list(range(len(train_losses))), train_losses, label="train", c="r")
        plt.plot(list(range(len(test_losses))), test_losses, label="test", c="g")
        plt.legend()
        plt.title(f"Loss Curve (epoch={epoch})")
        plt.xlabel("Epoch")
        plt.ylabel("Loss")

        plt.savefig(f"{plot_file_stem}_loss_curve.png")


if __name__ == "__main__":

    print("START", flush=True)
    mode = "combined"

    print(str(Path(".").resolve()))
    print(list(Path(".").glob("*")))

    folds, feat_dim, Y_dim, energy_scaler, grad_scaler = featurize.load_data("data_dir/data.pkl", mode)

    with open("results/energy_scaler.pkl", "wb") as pkl:
        pickle.dump(energy_scaler, pkl)

    with open("results/grad_scaler.pkl", "wb") as pkl:
        pickle.dump(grad_scaler, pkl)

    sizes = [
                (feat_dim, 256),

                (256, 256),

                (256, Y_dim)
            ]

    print("SIZES:")
    print(sizes, flush=True)
    fold_grad_rmses = []
    fold_grad_r2s= []
    fold_energy_rmses = []
    fold_energy_r2s= []

    with open("results/sizes.pkl", "wb") as pkl:
        pickle.dump(sizes, pkl)

    for fold, (train_dl, test_dl) in enumerate(folds):

        print(f"FOLD {fold} training...", flush=True)

        m = model.LinearModel(sizes)
        loss_fn = torch.nn.MSELoss()

        training_loop(m, fold, loss_fn, train_dl, test_dl, plot_file_stem=f"results/{mode}_{fold}", model_file=f"results/{mode}_fold_{fold}_best_model.pt")
        m.load_state_dict(torch.load(f"results/{mode}_fold_{fold}_best_model.pt", weights_only=True))

        grad_rmse, grad_r2 = grad_corelation(m, "final", test_dl, grad_scaler, plot_file_stem=f"results/grad_{mode}_final_{fold}")

        energy_rmse, energy_r2 = energy_corelation(m, "final", test_dl, energy_scaler, plot_file_stem=f"results/energy_{mode}_final_{fold}")

        fold_grad_rmses.append(grad_rmse)
        fold_energy_rmses.append(energy_rmse)
        fold_grad_r2s.append(grad_r2)
        fold_energy_r2s.append(energy_r2)

    fold_grad_rmses = np.array(fold_grad_rmses)
    fold_energy_rmses = np.array(fold_energy_rmses)
    fold_grad_r2s = np.array(fold_grad_r2s)
    fold_energy_r2s = np.array(fold_energy_r2s)

    print("Energy")
    print(f"5 Fold rmse: mean=", np.mean(fold_energy_rmses), "std=", np.std(fold_energy_rmses))
    print(f"5 Fold r2: mean=", np.mean(fold_energy_r2s), "std=", np.std(fold_energy_r2s))

    best = np.where(fold_energy_rmses == fold_energy_rmses.max())[0][0]
    print("Best fold was", best)

    best_model = f"results/{mode}_fold_{best}_best_model.pt"
    shutil.copy(best_model, "results/best_model.pt")

    print("Grad")
    print(f"5 Fold rmse: mean=", np.mean(fold_grad_rmses), "std=", np.std(fold_grad_rmses))
    print(f"5 Fold r2: mean=", np.mean(fold_grad_r2s), "std=", np.std(fold_grad_r2s))

    print("Best fold was", np.where(fold_grad_rmses == fold_grad_rmses.max()))

    print("Done")


