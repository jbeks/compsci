import sys
import numpy as np
import matplotlib.pyplot as plt
from nbody import *

def find_dists(data_body, data_center):
    min_dist = np.linalg.norm(data_body[0] - data_center[0])
    max_dist = 0
    for d, d0 in zip(data_body, data_center):
        dist = np.linalg.norm(d - d0)
        if dist < min_dist: min_dist = dist
        if dist > max_dist: max_dist = dist
    return min_dist, max_dist

def find_period(data_body, data_center):
    data = data_body - data_center
    prev_dist = np.linalg.norm(data[1] - data[0])
    back = False
    for i in range(2, len(data)):
        cur_dist = np.linalg.norm(data[i] - data[0])
        if not back:
            if cur_dist < prev_dist:
                back = True
        elif cur_dist > prev_dist:
            return i
        prev_dist = cur_dist
    return -1

def evaluation(center, all_other, dt):
    for i in range(len(all_other)):
        mi, ma = find_dists(all_other[i], center)
        p = find_period(all_other[i], center) * dt / 86400.
        print("Body " + str(i) + ":")
        print("  min dist to center = " + str(mi))
        print("  max dist to center = " + str(ma))
        print("  period of orbit    = " + str(p) + " (days)")

if __name__ == "__main__":
    itype, t_end, dt, t_dia, t_out, plot = parse_arguments()
    G, sys = get_system_data()
    system = System(G, sys, itype)
    sim_data = simulate(system, t_end, dt, t_dia, t_out)
    evaluation(sim_data[0], sim_data[1:], dt)
    if plot:
        simple_plot([p.T for p in sim_data])

