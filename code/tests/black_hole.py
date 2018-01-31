import argparse
import numpy as np
from warnings import warn

import code_dir
from nbody import *

# adds black hole arguments (speed and distance) to parser
def set_parser_bh(parser):
    parser.add_argument(
        "dist", type=float,
        help="distance of black hole from solar system"
    )
    parser.add_argument(
        "speed", type=float,
        help="speed of black hole (km/s)"
    )

# runs a solar system simulation
# with a black hole at the given distance and with the given speed
def simulate_bh(dist, speed, args, G, sys):
    # longest distance of black hole from solar system
    start_dist = 1.5e11

    # axis on which the black hole is placed
    e1 = np.array([0,1,0], dtype=float)
    e2 = np.array([0,0,1], dtype=float)
    e1 /= np.linalg.norm(e1)
    e2 /= np.linalg.norm(e2)

    # minimum heigt for a simulation (takes 1.2 * orbit neptune in time)
    min_height = 1.2 * 5201280000. * speed / 2
    # check whether given distance is smaller than maximum distance
    if start_dist < dist:
        warn("Given distance is larger than assumed largest distance")
        height = min_height
    else:
        # calculate height for simulation
        # (where dist from solar system is start_dist)
        height = np.sqrt(start_dist ** 2 - dist ** 2)
        # if height is less than min_height, set height to min_height
        if height < min_height:
            height = min_height

    # calculate position and velocity of black hole
    vec_dist = dist * e1
    vec_height = height * e2
    bh_p = vec_dist + vec_height
    bh_v = -e2 * speed

    # create black hole
    sun_m = 1.989e+30
    bh = Body(
        6.5 * sun_m,                                    # stellar black hole
        bh_p,
        bh_v,
        ("Black_Hole", "None")
    )

    # create system with black hole
    system = System(G, sys+[bh], args.itype.lower())
    # if no time is given, run for the time it takes
    # for the black hole to move 2 * height
    if args.t_end == 0:
        t_end = 2 * np.linalg.norm(vec_height) / speed
    else:
        t_end = args.t_end

    # return output of simulation
    return simulate(system, t_end, args.dt, args.t_dia, args.t_out)


if __name__ == "__main__":
    # create parser
    parser = argparse.ArgumentParser()
    set_parser(parser)
    set_parser_bh(parser)
    args = parser.parse_args()

    # get system from standard input
    G, sys = get_system_data()

    # run black hole simulation
    sim_data = simulate_bh(args.dist, args.speed, args, G, sys)
    # plot data if asked to
    if args.plot_2d or args.plot_3d:
        simple_plot([p.T for p in sim_data], args.plot_3d, args.n_points)

