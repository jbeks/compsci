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
    sun_m = 1.989e+30
    b = Body(
        6.5 * sun_m,                                    # stellar black hole
        [dist, speed * args.t_end / 2, 0],
        [0, -speed, 0]
    )
    system = System(G, sys+[b], args.itype.lower())
    return simulate(system, args.t_end, args.dt, args.t_dia, args.t_out)


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

