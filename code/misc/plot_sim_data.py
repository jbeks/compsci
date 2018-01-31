import argparse
import numpy as np

import code_dir
from nbody import simple_plot
from read_sim_data import read_sim_data
from read_short_sim_data import read_short_sim_data

if __name__ == "__main__":
    # create command line argument parser
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-p3d", "--plot_3d", action="store_true",
        help="plot simulation in 3 dimensions"
    )
    parser.add_argument(
        "-np", "--n_points", type=int, default=-1,
        help="maximum amount of points plotted"
    )
    parser.add_argument(
        "-ld", "--long_data", action="store_true",
        help="read data from long format"
    )
    args = parser.parse_args()
    # read data to be plotted
    if args.long_data:
        _, _, sim_data, _ = read_sim_data()
    else:
        _, _, sim_data = read_short_sim_data()
    # plot data
    simple_plot([np.array(p).T for p in sim_data[-2:]], args.plot_3d, args.n_points)

