import sys
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

if __name__ == "__main__":

   
    if len(sys.argv) < 2:
        raise ValueError("Must provide input csv file")

    df = pd.read_csv(sys.argv[1])

    df["fd_err_x"] = df["fd_x"] - df["analytical_x"]
    df["cd_err_x"] = df["cd_x"] - df["analytical_x"]

    df["fd_err_y"] = df["fd_y"] - df["analytical_y"]
    df["cd_err_y"] = df["cd_y"] - df["analytical_y"]

    df["fd_err_z"] = df["fd_z"] - df["analytical_z"]
    df["cd_err_z"] = df["cd_z"] - df["analytical_z"]

    df["fd_force_err"] = df.apply(lambda r: np.linalg.norm((r["fd_err_x"], r["fd_err_y"], r["fd_err_z"])), axis=1)
    df["cd_force_err"] = df.apply(lambda r: np.linalg.norm((r["cd_err_x"], r["cd_err_y"], r["cd_err_z"])), axis=1)


    fd_averages = {}
    cd_averages = {}


    for atom_idx, atom_df in df.groupby("atom_idx"):

        xs = atom_df["h"].tolist()
        y_fd = atom_df["fd_force_err"].tolist()
        y_cd = atom_df["cd_force_err"].tolist()

        plt.plot(xs, y_fd, color="blue")
        plt.plot(xs, y_cd, color="orange")

        for h, err in zip(xs, y_fd):
            if not h in fd_averages:
                fd_averages[h] = []

            fd_averages[h].append(err)

        for h, err in zip(xs, y_cd):
            if not h in cd_averages:
                cd_averages[h] = []

            cd_averages[h].append(err)


    plt.xscale("log")
    plt.xlabel("Finite Difference Step Size (log)")
    plt.yscale("log")
    plt.ylabel("Error compared to Analytical Force (log)")


    cd_averages = {h : np.mean(errs) for h,errs in cd_averages.items()}
    cd_hs, cd_errs = zip(*cd_averages.items())
    cd_hs = np.array(cd_hs)
    cd_errs = np.array(cd_errs)
    cd_hs_log = np.log10(cd_hs)
    cd_errs_log = np.log10(cd_errs)

    cd_valid = np.isfinite(cd_errs_log)
    cd_errs_log = cd_errs_log[cd_valid]
    cd_hs_log = cd_hs_log[cd_valid]

    fd_averages = {h : np.mean(errs) for h,errs in fd_averages.items()}
    fd_hs, fd_errs = zip(*fd_averages.items())
    fd_hs = np.array(fd_hs)
    fd_errs = np.array(fd_errs)
    fd_hs_log = np.log10(fd_hs)
    fd_errs_log = np.log10(fd_errs)

    fd_valid = np.isfinite(fd_errs_log)
    fd_errs_log = fd_errs_log[fd_valid]
    fd_hs_log = fd_hs_log[fd_valid]


    fd_m, fd_b = np.polyfit(fd_hs_log, fd_errs_log, 1)
    cd_m, cd_b = np.polyfit(cd_hs_log, cd_errs_log, 1)

    plt.plot(10**np.log10(fd_hs), 10**(fd_m * np.log10(fd_hs) + fd_b), color="blue", linestyle="--")
    plt.plot(10**np.log10(cd_hs), 10**(cd_m * np.log10(cd_hs) + cd_b), color="orange", linestyle="--")

    plt.legend(handles=[
                         Line2D([0], [0], color="blue", label="Forward Difference"),
                         Line2D([0], [0], color="blue", label=f"Forward Difference - Regression Slope={fd_m:.2f}", linestyle="--"),
                         Line2D([0], [0], color="orange", label="Central Difference"),
                         Line2D([0], [0], color="orange", label=f"Central Difference - Regression Slope={cd_m:.2f}", linestyle="--"),
               ])

    plt.savefig(sys.argv[1].split(".")[0] + ".png")


