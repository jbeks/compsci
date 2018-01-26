import code_dir
from nbody import *

def simulate_bh(dist, speed, args, G, sys):
    sun_m = 1.989e+30
    b = Body(
        6.5 * sun_m,                                    # stellar black hole
        [dist, speed * args[1] / 2, 0],
        [0, -speed, 0]
    )
    system = System(G, sys+[b], args[0])
    return simulate(system, *args[1:5])


if __name__ == "__main__":
    args = parse_arguments()
    G, sys = get_system_data()
    system = System(G, sys, args[0])
    sim_data_orig = simulate(system, *args[1:5])

    dist = 5e8
    speed = 2000

    sim_data = simulate_bh(dist, speed, args, G, sys)
    if args[5] or args[6]:
        simple_plot(
            [p.T for p in sim_data], args[6], args[7]
        )
#        simple_plot(
#            [p.T for p in list(np.array(sim_data[:-1]) \
#                - np.array(sim_data_orig)) + sim_data[:-1]],
#            args[6], args[7]
#        )

