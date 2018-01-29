import argparse
import numpy as np

import code_dir
from nbody import simple_plot
from read_sim_data import read_sim_data

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-p3d", "--plot_3d", action="store_true",
        help="plot simulation in 3 dimensions"
    )
    parser.add_argument(
        "-np", "--n_points", type=int, default=-1,
        help="maximum amount of points plotted"
    )
    args = parser.parse_args()
    _, _, sim_data, _ = read_sim_data()
    simple_plot([np.array(p).T for p in sim_data], args.plot_3d, args.n_points)

