import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt

import code_dir
from nbody import *

# find longest and shortest distance of data_body to data_center
def find_dists(data_body, data_center):
    # initialize variables
    min_dist = np.linalg.norm(data_body[0] - data_center[0])
    max_dist = 0
    # for each pair of positions
    for d, d0 in zip(data_body, data_center):
        # determine distance
        dist = np.linalg.norm(d - d0)
        # and look for new longest and shortest distance
        if dist < min_dist: min_dist = dist
        if dist > max_dist: max_dist = dist
    return min_dist, max_dist

# find period of fist orbit of data_body around data_center
def find_period(data_body, data_center):
    # determine relative positions
    data = data_body - data_center
    # determine distance between datapoint and fist datapoint
    prev_dist = np.linalg.norm(data[1] - data[0])
    back = False
    for i in range(2, len(data)):
        # determine distance between datapoint and first datapoint
        cur_dist = np.linalg.norm(data[i] - data[0])
        # if not having found the half way point
        if not back:
            # check for half way point
            if cur_dist < prev_dist:
                back = True
        # else check for completed orbit
        elif cur_dist > prev_dist:
            return i
        # update previous distance
        prev_dist = cur_dist
    return -1

# find and print minimum and maximum distances to center and orbit period
# for all bodies in all_other
def evaluation(center, all_other, dt):
    for i in range(len(all_other)):
        # find minimum and maximum distance to center
        mi, ma = find_dists(all_other[i], center)
        # find period
        p = find_period(all_other[i], center) * dt / 86400.
        # print results
        print("Body " + str(i) + ":")
        print("  min dist to center = " + str(mi))
        print("  max dist to center = " + str(ma))
        print("  period of orbit    = " + str(p) + " (days)")

if __name__ == "__main__":
    # create parser
    parser = argparse.ArgumentParser()
    set_parser(parser)
    args = parser.parse_args()

    # create system
    G, sys = get_system_data()
    system = System(G, sys, args.itype.lower())

    # run simulation and evaluation
    sim_data = simulate(system, args.t_end, args.dt, args.t_dia, args.t_out)
    evaluation(sim_data[0], sim_data[1:], args.dt)
    # plot data if asked to
    if args.plot_2d or args.plot_3d:
        simple_plot([p.T for p in sim_data], args.plot_3d, args.n_points)

