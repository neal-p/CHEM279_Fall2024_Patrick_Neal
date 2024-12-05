import model
import torch
import featurize
from tqdm import tqdm
from matplotlib import pyplot as plt
import numpy as np
from sklearn.metrics import root_mean_squared_error


TQDM_PARAMS = {"leave":False}


def train_step(m, loss_fn, opt, dl):

    m.train(True)
    total_loss = 0.0
    count = 0
    
    for i, data in tqdm(enumerate(dl), 
                        **TQDM_PARAMS):

        feats, Y = data
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
        for i, data in tqdm(enumerate(dl), 
                            **TQDM_PARAMS):

            feats, Y = data
            y_hat = m(feats)
            count += len(feats)

            loss = loss_fn(y_hat, Y)

           # print("feats")
           # print(feats)
           # print()
           # print("yhat")
           # print(y_hat)
           # print()
           # print("Y")
           # print(Y)
           # print()
           # print(loss)
           # print(flush=True)

            total_loss += loss.item()

    return total_loss / count


def train(m, loss_fn, train_dl, test_dl, epochs=200, plot_file="train_progress", lr=1e-4):

    opt = torch.optim.Adam(m.parameters(), lr=lr)

    train_losses = []
    test_losses = []

    best_loss = 10000

    for i in range(epochs+1):

        # Get pure initialized weights' test performance
        # without doing any gradient updates
        if i == 0:
            avg_train_loss = test_step(m, loss_fn, train_dl)
            avg_test_loss = test_step(m, loss_fn, test_dl)

        else:
            avg_train_loss = train_step(m, loss_fn, opt, train_dl)
            avg_test_loss = test_step(m, loss_fn, test_dl)

        train_losses.append(avg_train_loss)
        test_losses.append(avg_test_loss)


        if avg_test_loss < best_loss:
            best_loss = avg_test_loss
            print("lowered loss!")
            torch.save(m.state_dict(), "best_model.pt")

        print(f"    EPOCH: {i}, train_loss={avg_train_loss:.5f}, test_loss={avg_test_loss:.5f}", flush=True)

        # with un-initialized
        plt.plot(list(range(i+1)), train_losses, label="train", c="r")
        plt.plot(list(range(i+1)), test_losses, label="test", c="b")
        plt.legend()
        plt.title(f"MSE Loss")
        plt.xlabel("Epoch")
        plt.ylabel("Avg loss")

        plt.savefig(plot_file + "_untrained.png")
        plt.close()
        plt.clf()

       # without un-initialized
        plt.plot(list(range(i)), train_losses[1:], label="train", c="r")
        plt.plot(list(range(i)), test_losses[1:], label="test", c="b")
        plt.legend()
        plt.title(f"MSE Loss")
        plt.xlabel("Epoch")
        plt.ylabel("Avg loss")

        plt.savefig(plot_file + ".png")
        plt.close()
        plt.clf()


    return m


def eval(m, dl, energy_scaler, plot_file):

    print("Loaded best model")
    m.load_state_dict(torch.load("best_model.pt"))

    true_energies = []
    pred_energies = []

    true_grad_norms = []
    pred_grad_norms = []

    with torch.no_grad():
        for feats, Y in dl:

            for row in m(feats):
                pred_etot = row[0].item()
               # pred_grad = row[1:].detach().numpy().reshape((int((len(row)-1)/3), 3))
               # pred_grad_norm = np.linalg.norm(pred_grad)

                pred_energies.append(pred_etot)
               # pred_grad_norms.append(pred_grad_norm)


            for row in Y:
                true_etot = row[0].item()
               # true_grad = row[1:].detach().numpy().reshape((int((len(row)-1)/3), 3))
               # true_grad_norm = np.linalg.norm(true_grad)

                true_energies.append(true_etot)
               # true_grad_norms.append(true_grad_norm)

    true_energies = np.array(true_energies).reshape((-1, 1))
    pred_energies = np.array(pred_energies).reshape((-1, 1))

    true_energies = energy_scaler.inverse_transform(true_energies).reshape(-1)
    pred_energies = energy_scaler.inverse_transform(pred_energies).reshape(-1)
     
    etot_rmse = root_mean_squared_error(true_energies, pred_energies)
   # grad_norm_rmse = root_mean_squared_error(true_grad_norms, pred_grad_norms)

    print(f"        etot rmse      ={etot_rmse:.5f}")
   # print(f"        grad_norm rmse ={grad_norm_rmse:.5f}", flush=True)

    plt.scatter(true_energies, pred_energies)
    plt.title(f"True Energy vs Predicted Energy: RMSE={etot_rmse:.5f}")
    plt.xlabel("True Energy")
    plt.ylabel("Predicted Energy")

    plt.savefig(plot_file)
    plt.close()
    plt.clf()

    print("Differences")
    print(true_energies)
    print(true_energies - pred_energies)
    print(flush=True)



if __name__ == "__main__":

    folds, df, energy_scaler = featurize.load_data("prepare_training_data")
    feat_dim = folds[0][0].feat_dim
    Y_dim = folds[0][0].Y_dim

    sizes = [
                (feat_dim, 100),
                
                (100, 500),
                (500, 500),
                (500, 500),
                (500, 500),
                (500, 100),

                (100, Y_dim)
            ]

    print("Layer sizes:")

    layer = 0
    for input, output in sizes:
        print(f"Layer{layer}: ({input}, {output})")
        layer += 1

    m = model.LinearModel(sizes=sizes)
    loss_fn = torch.nn.L1Loss()

    for fold, (train_dl, test_dl) in enumerate(folds):

        print(f"Starting FOLD {fold}...")
        print(f"    train_dl size={len(train_dl)}, test_dl size={len(test_dl)}", flush=True)

        m = train(m, loss_fn, train_dl, test_dl, plot_file=f"fold{fold}_training_graph.png")

        print(f"FOLD {fold} training complete...")
        print("    Train eval")
        eval(m, train_dl, energy_scaler, f"fold{fold}_train_energy_correlation.png")

        print("    Test eval")
        eval(m, test_dl, energy_scaler, f"fold{fold}_test_energy_correlation.png")

        break




