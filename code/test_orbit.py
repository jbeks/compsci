import sys
import numpy as np
import matplotlib.pyplot as plt
from nbody import *

def v_error(data, index):
    return np.linalg.norm(data[0]) - np.linalg.norm(data[index])

def pos_data(data, orbit_i):
    return

def find_dists(data_body, data_center):
    min_dist = np.linalg.norm(data_body[0] - data_center[0])
    max_dist = 0
    for d, d0 in zip(data_body, data_center):
        dist = np.linalg.norm(d - d0)
        if dist < min_dist: min_dist = dist
        if dist > max_dist: max_dist = dist
    return min_dist, max_dist

def find_period(data_body, data_center, start):
    data = data_body - data_center
    prev_dist = np.linalg.norm(data[start + 1] - data[start])
    back = False
    for i in range(start + 2, len(data)):
        cur_dist = np.linalg.norm(data[i] - data[start])
        if not back:
            if cur_dist < prev_dist:
                back = True
        elif cur_dist > prev_dist:
            return i
        prev_dist = cur_dist
    return -1

def orbit_data(data, center):
    orbit_i = 0
    orbit_err = []
    while orbit_i < len(data) and orbit_i != -1:
        orbit_i = find_period(data, center, orbit_i)
        orbit_err.append(orbit_i)
    avg_err = abs(orbit_err[0] - np.mean(orbit_err))
    max_err = abs(orbit_err[0] - max(orbit_err))
    return orbit_err[0], avg_err, max_err, orbit_err

def evaluation(data, v_data, dt, system):
    find = lambda x: next((i for i in range(len(system)) if system[i].name == x.parent), None)
    for i in range(1, len(data)):
        body = system[i]
        p_index = find(body)
        mi, ma = find_dists(data[i], data[p_index])
        orbit_i, avg_err, max_err, points = orbit_data(data[i], data[p_index])
        p = orbit_i * dt / 86400.
        pos_diff = pos_data(data, points)
        p_err = np.linalg.norm(data[i][orbit_i] - data[i][0]) if orbit_i != 0 else -1
        v_err = v_error(v_data[i], orbit_i)
        print("Body " + str(i) + " (" + str(system[i]) + "):")
        print("  min dist to " + str(system[p_index]) + " = " + str(mi) + " (km)")
        print("  max dist to " + str(system[p_index]) + " = " + str(ma) + " (km)")
        if p > 0:
            print("  period of orbit: first " + str(p) + " and avg " + str(np.mean(points) * dt / 84600.) + " (days)")
            print("  std: " + str(np.std(points) * dt / 84600.) + " (s) and max: " + str(max(points) * dt / 84600.) + " (s) and min: " + str(min(points) * dt / 84600.) + " (s)")
            print("  orbit avg error: " + str(avg_err) + " (s) and max error " + str(max_err) + " (s) in " + str(len(points)) + " data points")
        else:
            print("  period of orbit not found")
        print("  positional error = " + str(p_err) + " (km)")
        print("  velocity error = " + str(v_err) + " (km/s)")

if __name__ == "__main__":
    orbit = {}
    itype, t_end, dt, t_dia, t_out, plot, planet, years = parse_arguments()
    if t_end == 1:
        t_end = years * orbit[planet] * 24 * 60 * 60
    G, sys = get_system_data()
    system = System(G, sys, itype)
    sim_data, v_data = simulate(system, t_end, dt, t_dia, t_out)
    evaluation(sim_data, v_data, dt, system.sys)
    if plot:
        simple_plot([p.T for p in sim_data])
        # full_plot([p.T for p in sim_data])
