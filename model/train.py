import model
import torch
import featurize
from tqdm import tqdm


TQDM_PARAMS = {}


def train_step(m, loss_fn, opt, dl):

    m.train(True)
    total_loss = 0.0
    
    for i, data in tqdm(enumerate(dl), 
                        **TQDM_PARAMS):

        feats, Y = data
        opt.zero_grad()

        y_hat = m(feats)

        loss = loss_fn(y_hat, Y)
        loss.backward()

        opt.step()

        total_loss += loss.item()

    return total_loss / i


def test_step(m, loss_fn, dl):

    m.eval()
    total_loss = 0.0
   
    with torch.no_grad():
        for i, data in tqdm(enumerate(dl), 
                            **TQDM_PARAMS):

            feats, Y = data
            y_hat = m(feats)

            loss = loss_fn(y_hat, Y)
            total_loss += loss.item()

    return total_loss / i


def train(m, loss_fn, train_dl, test_dl, epochs=100, plot_file="train_progress.png"):

    opt = torch.optim.Adam(m.parameters(), lr=1e-3)

    train_losses = []
    test_losses = []

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

        print(f"    EPOCH: {i}, train_loss={avg_train_loss:.5f}, test_loss={avg_test_loss:.5f}", flush=True)

        plt.plot(list(range(i)), train_losses, label="train", c="r")
        plt.plot(list(range(i)), test_losses, label="test", c="b")
        plt.legend()
        plt.title(f"MSE Loss")
        plt.xlabel("Epoch")
        plt.ylabel("Avg loss")

        plt.savefig(plot_file)
        plt.close()
        plt.clf()


def eval(m, dl):

    true_energies = []
    pred_energies = []

    true_grad_norms = []
    pred_grad_norms = []

    for feats, Y in dl:

        for row in m(feats):
            pred_etot, *pred_grad = row
            pred_grad = np.array(pred_grad).reshape((None, 3))
            pred_grad_norm = np.linalg.norm(pred_grad)

            pred_energies.append(pred_etot)
            pred_grad_norms.append(pred_grad_norm)
            

        for row in Y:
            true_etot, *true_grad = row
            true_grad = np.array(true_grad).reshape((None, 3))
            true_grad_norm = np.linalg.norm(true_grad)

            true_energies.append(true_etot)
            true_grad_norms.append(true_grad_norm)
 
    etot_rmse = root_mean_squared_error(true_energies, pred_energies)
    grad_norm_rmse = root_mean_squared_error(true_grad_norms, pred_grad_norms)

    print(f"        etot rmse      ={etot_rmse:.5f}")
    print(f"        grad_norm rmse ={etot_rmse:.5f}", flush=True)


if __name__ == "__main__":

    folds = featurize.load_data("prepare_training_data")
    feat_dim = folds[0][0].feat_dim
    Y_dim = folds[0][0].Y_dim

    sizes = [
                (feat_dim, 100),
                
                (100, 200),
                (200, 100),

                (100, Y_dim)
            ]

    m = model.LinearModel(sizes=sizes)
    loss_fn = torch.nn.MSELoss()

    for fold, (train_dl, test_dl) in enumerate(folds):

        print(f"Starting FOLD {fold}...")
        print(f"    train_dl size={len(train_dl)}, test_dl size={len(test_dl)}", flush=True)

        m = train(m, loss_fn, train_dl, test_dl, plot_file=f"fold{fold}_training_graph.png")

        print(f"FOLD {fold} training complete...")
        print("    Train eval")
        eval(m, train_dl)

        print("    Test eval")
        eval(m, test_dl)



 


