import os
import sys
import argparse
import tempfile
import numpy as np
import shutil
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes

import code_dir
from nbody import *

def set_parser_m(parser):
    parser.add_argument(
        "-v", "--verbose", action="store_true",
        help="print verbose energy calculations"
    )

def read_ref(name):
    data = []
    stepsize = 0
    with open("ref/" + name + ".ref", 'r') as f:
        for line in f:
            if not stepsize and line.startswith("stepsize"):
                stepsize = int(line.split()[1])
                continue
            part = line.split()
            data.append([float(part[0]), float(part[1]), float(part[2])])
    return data, stepsize

def test_ref(data, body, dt):
    ref_data, stepsize = read_ref(body)
    plot = []
    for index in range(len(data)):
        if (index * dt) % stepsize == 0:
            # print("DAY: " + str(int((index * dt) / (60 * 60 * 24))))
            point = np.linalg.norm(np.array(ref_data)[int((index * dt) / stepsize)] - np.array(data)[index])
            plot.append(point)
    return plot

def evaluation(data, dt, system):
    diff_data = []
    find = lambda x: next((i for i in range(len(system)) if system[i].name == x.parent), None)
    for i in range(len(data)):
        body = system[i]
        print("body " + str(i))
        print("name " + str(system[i].name))
        p_index = find(body)
        if p_index != None:
            print("parent " + str(system[p_index].name))
        ref_check = test_ref(data[i], body.name, dt)
        diff_data.append(ref_check)
        print()
    return diff_data

def plot_diff(data, methods, t_end):
    # http://akuederle.com/matplotlib-zoomed-up-inset

    t_energy = np.linspace(0, t_end, len(data[3][0]))
    t_refcheck = np.linspace(0, t_end, len(data[3][1][0]))

    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()

    ax1.set_xlabel('time (s)')
    ax1.set_ylabel('ref error')
    ax2.set_ylabel('energy error')

    lines1 = []
    for i in range(len(methods)):
        line, = ax1.plot(t_refcheck, data[i][1][0], c='C'+str(i % len(methods)), linestyle="-")
        lines1.append(line)
    ax1.legend(lines1, methods, loc=1)

    lines2 = []
    for i in range(1, len(methods)):
        line, = ax2.plot(t_energy, data[i][0], c='C'+str(i % len(methods)), linestyle="-.")
        lines2.append(line)
    ax2.legend(lines2, methods[1:], loc=3)

    plt.grid(True)

    axins = zoomed_inset_axes(ax1, 5000, loc=2)
    for i in range(len(methods)):
        line, = axins.plot(t_refcheck, data[i][1][0], c='C'+str(i % len(methods)), linestyle="-")
    x1, x2, y1, y2 = (1.576e9 + 760000), (1.576e9 + 840000), 1.99585e7, 1.9961e7
    axins.set_xlim(x1, x2)
    axins.set_ylim(y1, y2)
    plt.yticks(visible=False)
    plt.xticks(visible=False)

    mark_inset(ax1, axins, loc1=1, loc2=4, fc="none", ec="0.5")

    plt.grid(True)
    plt.show()

if __name__ == "__main__":
    # read arguments
    parser = argparse.ArgumentParser()
    set_parser(parser)
    set_parser_m(parser)
    args = parser.parse_args()

    # copy input file for multiple uses
    with open(".tmp.input", 'w+') as f:
        shutil.copyfileobj(sys.stdin, f)

    # check if we need to test just one or all of the methods
    methods = ["euler", "verlet", "rk4", "hermite"]
    if args.itype.lower() != "all":
        methods = [args.itype.lower()]

    # run simulation with the preferred method(s) & save the data
    method_data = []
    for method in methods:
        with open(".tmp.input", 'r+') as f:
            G, sys = get_system_data(f, True)
            system = System(G, sys, method)
            print("itype " + str(system.itype) + "\n")
            sim_data = simulate(system, args.t_end, args.dt, args.t_dia,
                                args.t_out, args.verbose)
            ref_diff = evaluation(sim_data, args.dt, system.sys)
            method_data.append([system.e_diff, ref_diff])

    # delete the extra input file
    os.remove(".tmp.input")

    plot_diff(method_data, methods, args.t_end)

    if args.plot_2d or args.plot_3d:
        simple_plot([p.T for p in sim_data], args.plot_3d, args.n_points)
