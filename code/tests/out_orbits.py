import numpy as np
import argparse
import math
import os
import sys

import code_dir
from read_short_sim_data import read_short_sim_data

<<<<<<< HEAD
def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)

def angle(current, last):
    v1_u = unit_vector(current)
    v2_u = unit_vector(last)
    return math.degrees(np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0)))

=======
>>>>>>> 9b2a56a97e3968c61642bbe5b39b260d2065541b
def find_period(data_body, data_center, start):
    data = data_body - data_center
    prev_dist = np.linalg.norm(data[start + 1] - data[start])
    orig = prev_dist
    back = False
    for i in range(start + 2, len(data)):
        cur_dist = np.linalg.norm(data[i] - data[start])
        if not back:
            if cur_dist < prev_dist:
                back = True
        elif cur_dist > prev_dist: # and angle(cur_dist, orig) < 30:
            return i
        prev_dist = cur_dist
    return -1

def orbit_data(data, center):
    orbit_i = 0
    orbit_err = []
    while orbit_i < len(data) - 1 and orbit_i != -1:
        orbit_err.append(orbit_i)
        orbit_i = find_period(data, center, orbit_i)
    return orbit_err

def phelions(points, center_points):
    min_d = np.linalg.norm(points[0] - center_points[0])
    max_d = 0
    for i in range(1, len(points)):
        dist = np.linalg.norm(points[i] - center_points[i])
        if dist < min_d:
            min_d = dist
        if dist > max_d:
            max_d = dist
    return max_d, min_d

<<<<<<< HEAD
def no_orbit(planet_data, center_data, i, ori_pheli):
    aphelion_check = (np.linalg.norm(planet_data[i] - center_data[i]) > 2 * ori_pheli)
    angle_check = True
    if i > 300:
        old_position = planet_data[i - 300] - planet_data[i - 301]
        new_position = planet_data[i] - planet_data[i-1]
        angle_check = (angle(old_position, new_position) < 0.5)
    return aphelion_check and angle_check

def out_orbit_time(orbit_num, planet_data, center_data, ori_pheli, time):
    orbit_indices = orbit_data(planet_data, center_data)
    for i in range(orbit_indices[orbit_num], len(planet_data)):
        if no_orbit(planet_data, center_data, i, ori_pheli[orbit_num][0]):
            return time[i]
=======
def out_orbit_time(orbit_num, planet_data, center_data, ori_pheli, dt):
    orbit_indices = orbit_data(planet_data, center_data)
    for i in range(orbit_indices[orbit_num], len(planet_data)):
        if (
            np.linalg.norm(planet_data[i]-center_data[i]) \
            > 2 * (ori_pheli[orbit_num][0])
        ):
            return i * dt
>>>>>>> 9b2a56a97e3968c61642bbe5b39b260d2065541b
    return -1

def orbit_sta(new_phelis, ori_pheli, planet_data, center_data, dt):
    stabil = []
    if len(ori_pheli) > len(new_phelis):
        orbit_num = len(new_phelis)
        return out_orbit_time(
            orbit_num, planet_data, center_data, ori_pheli, dt
        )
    return -2

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "file", type=str,
        help="input file name"
    )
    parser.add_argument(
        "-init", "--initialization", action="store_true",
        help="run initialization"
    )
    parser.add_argument(
        "-si", "--sun_index", type=int, default=0,
        help="index of sun in the data"
    )
    parser.add_argument(
        "-bhi", "--black_hole_index", type=int, default=-1,
        help="index of black hole in the data"
    )
    parser.add_argument(
        "-rd", "--resdir", type=str, default="",
        help="path to directory to store results"
    )
    args = parser.parse_args()

    time, _, sim_data = read_short_sim_data(args.file)
    sim_data = np.array(sim_data)
    phel_list = []
    if args.initialization:
        sim_data_subset = [
            sim_data[i] for i in range(len(sim_data))
            if i != (args.sun_index % len(sim_data))
        ]
    else:
        sim_data_subset = [
            sim_data[i] for i in range(len(sim_data))
            if i != (args.sun_index % len(sim_data))
            and i != (args.black_hole_index % len(sim_data))
        ]
    for planet in sim_data_subset:
        cur_phel = []
        orbit_list = orbit_data(planet, sim_data[args.sun_index])
        print("orbits: " + str(len(orbit_list) - 1))
        for i in range(len(orbit_list) - 1):
            cur_phel.append(
                (
                    phelions(
                        planet[orbit_list[i]:orbit_list[i + 1]],
                        sim_data[0][orbit_list[i]:orbit_list[i + 1]]
                    ),
                    time[orbit_list[i]]
                )
            )
        phel_list.append(cur_phel)
    sep = "\\" if os.name == "nt" else "/"
    outdir = os.path.abspath(
        os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
    ) + sep + "output" + sep + "out"
    reff = outdir + sep + "ori_list_data"
    if args.initialization:
        np.save(reff, phel_list)
    else:
        ori_list = np.load(reff + ".npy")
        orbit_stabi = []
        for i in range(len(phel_list)):
            sys.stderr.write("planet no. " + str(i) + "\n")
            ori_list_subset = np.array([
                ori_list[i][j][0] for j in range(len(ori_list[i]))
                if ori_list[i][j][1] <= time[-1]
            ])
            phel_list_i = np.array([lst[0] for lst in phel_list[i]])
            orbit_stabi.append(
                orbit_sta(
                    phel_list_i, ori_list_subset,
                    sim_data_subset[i], sim_data[0],
                    time[1] - time[0]
                )
            )
<<<<<<< HEAD
            orbit_stabi.append('next_planet')
        print(orbit_stabi)
=======
        fname = os.path.splitext(os.path.basename(args.file))[0]
        if args.resdir == "":
            resdir = outdir
        else:
            resdir = args.resdir
        with open(resdir + sep + fname + ".res", "w") as f:
            for o in orbit_stabi:
                f.write(str(o) + "\n")

>>>>>>> 9b2a56a97e3968c61642bbe5b39b260d2065541b
