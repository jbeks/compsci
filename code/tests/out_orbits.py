import numpy as np
import argparse
import math
import os
import sys

import code_dir
from read_short_sim_data import read_short_sim_data

# out_orbits.py:
# compares a simulation with black hole
# to one without and checks the differences
# between the orbits of the planets
# and whether they have detached from orbit

def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)

def angle(current, last):
    """ Calculates the angle between two vectors. """
    u1 = unit_vector(current)
    u2 = unit_vector(last)
    return math.degrees(np.arccos(np.clip(np.dot(u1, u2), -1.0, 1.0)))

def find_period(data_body, data_center, start):
    """
    Find the timestep where the body 'data_body' orbiting
    around 'data_center', starting on timestep 'start,
    has made a full orbit.
    """
    data = data_body - data_center
    prev_dist = np.linalg.norm(data[start + 1] - data[start])
    orig = prev_dist
    back = False
    # loop through the body datapoints per timestep
    for i in range(start + 2, len(data)):
        cur_dist = np.linalg.norm(data[i] - data[start])
        # check when orbit is halfway done
        if not back:
            if cur_dist < prev_dist:
                back = True
        # return timestep when orbit is done
        elif cur_dist > prev_dist:
            return i
        prev_dist = cur_dist
    return -1

def orbit_data(data, center):
    """
    Find the orbit locations in all the datapoints
    of body 'data' orbiting around body 'center'.
    """
    orbit_i = 0
    orbit_err = []
    while orbit_i < len(data) - 1 and orbit_i != -1:
        orbit_err.append(orbit_i)
        orbit_i = find_period(data, center, orbit_i)
    return orbit_err

def phelions(points, center_points):
    """
    Find the aphelion and perihelion (smallest and biggest
    distance to the orbit center) of all bodies.
    """
    min_d = np.linalg.norm(points[0] - center_points[0])
    max_d = 0
    for i in range(1, len(points)):
        dist = np.linalg.norm(points[i] - center_points[i])
        if dist < min_d:
            min_d = dist
        if dist > max_d:
            max_d = dist
    return max_d, min_d

def no_orbit(planet_data, center_data, i, ori_pheli):
    """ Check whether a planet has detached from orbit. """
    aphelion_check = (np.linalg.norm(planet_data[i] - center_data[i]) > 2 * ori_pheli)
    angle_check = True
    if i > 300:
        old_position = planet_data[i - 300] - planet_data[i - 301]
        new_position = planet_data[i] - planet_data[i-1]
        angle_check = (angle(old_position, new_position) < 0.5)
    return aphelion_check and angle_check

def out_orbit_time(orbit_num, planet_data, center_data, ori_pheli, dt):
    """ Find the timestep of detachment of orbit for body 'planet_data'. """
    orbit_indices = orbit_data(planet_data, center_data)
    for i in range(orbit_indices[orbit_num], len(planet_data)):
        if no_orbit(planet_data, center_data, i, ori_pheli[orbit_num][0]):
            return i * dt
    return -1

def orbit_sta(new_phelis, ori_pheli, planet_data, center_data, dt):
    """
    Check the amount of orbits of body 'planet_data' with the reference simulation
    and if this amount differs, find out whether the planet has detached.
    """
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

    # read in the simulation data
    time, _, sim_data = read_short_sim_data(args.file)
    sim_data = np.array(sim_data)
    phel_list = []
    # if there is no reference simulation file, create it
    if args.initialization:
        sim_data_subset = [
            sim_data[i] for i in range(len(sim_data))
            if i != (args.sun_index % len(sim_data))
        ]
    # run the tests and find out the differences in orbit
    else:
        sim_data_subset = [
            sim_data[i] for i in range(len(sim_data))
            if i != (args.sun_index % len(sim_data))
            and i != (args.black_hole_index % len(sim_data))
        ]
    # go through all orbiting bodies
    for planet in sim_data_subset:
        cur_phel = []
        # find the timesteps of completed orbits
        orbit_list = orbit_data(planet, sim_data[args.sun_index])
        print("orbits: " + str(len(orbit_list) - 1))
        # find the aphelion and perihelion for every orbit
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

    # save the found datapoints if it is initialisation
    if args.initialization:
        np.save(reff, phel_list)
    # else compare the datapoints with the reference
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
        fname = os.path.splitext(os.path.basename(args.file))[0]
        if args.resdir == "":
            resdir = outdir
        else:
            resdir = args.resdir
        with open(resdir + sep + fname + ".res", "w") as f:
            for o in orbit_stabi:
                f.write(str(o) + "\n")
