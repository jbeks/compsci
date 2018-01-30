import argparse
import numpy as np
from warnings import warn

import code_dir
from nbody import *


def set_parser_bh(parser):
    parser.add_argument(
        "dist", type=float,
        help="distance of black hole from solar system"
    )
    parser.add_argument(
        "speed", type=float,
        help="speed of black hole (km/s)"
    )


def simulate_bh(dist, speed, args, G, sys):
    start_dist = 1.5e11

    e1 = np.array([0,1,0], dtype=float)
    e2 = np.array([0,0,1], dtype=float)

    e1 /= np.linalg.norm(e1)
    e2 /= np.linalg.norm(e2)

    # for time of 2.5 * orbit neptune
    min_height = 1.2 * 5201280000. * speed / 2
    if start_dist < dist:
        warn("Given distance is larger than assumed largest distance")
        height = min_height
    else:
        height = np.sqrt(start_dist ** 2 - dist ** 2)
        if height < min_height:
            height = min_height

    vec_dist = dist * e1
    vec_height = height * e2
    bh_p = vec_dist + vec_height
    bh_v = -e2 * speed

    sun_m = 1.989e+30
    bh = Body(
        6.5 * sun_m,                                    # stellar black hole
        bh_p,
        bh_v,
        ("Black_Hole", "None")
    )

    system = System(G, sys+[bh], args.itype.lower())
    if args.t_end == 0:
        t_end = 2 * np.linalg.norm(vec_height) / speed
    else:
        t_end = args.t_end

    return simulate(system, t_end, args.dt, args.t_dia, args.t_out)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    set_parser(parser)
    set_parser_bh(parser)
    args = parser.parse_args()

    G, sys = get_system_data()
    system = System(G, sys, args.itype.lower())

#    sim_data_orig = simulate(system, args.t_end, args.dt, args.t_dia, args.t_out)

    sim_data = simulate_bh(args.dist, args.speed, args, G, sys)
    if args.plot_2d or args.plot_3d:
        simple_plot([p.T for p in sim_data], args.plot_3d, args.n_points)

