from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
import matplotlib.pyplot as plt
import numpy as np
import argparse
import tempfile
import shutil
import os
import sys

import code_dir
from nbody import *

def set_parser_m(parser):
    parser.add_argument(
        "-v", "--verbose", action="store_true",
        help="print verbose energy calculations"
    )

def read_ref(name):
    """
    Read in the reference datafile
    of planet with name 'name'.
    """
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

def read_all_refs(sys):
    """
    Read in all reference datafiles.
    """
    data = []
    stepsizes = []
    for body in sys:
        ref_data, stepsize = read_ref(body.name)
        data.append(ref_data)
        stepsizes.append(stepsize)
    return data, stepsizes

def test_ref(data, ref_data, stepsize, dt):
    """
    Compare the simulation data with the reference files.
    """
    if len(ref_data) < len(data) * dt // stepsize:
        raise IndexError("Insufficient reference data")
    plot = []
    # check all bodies
    for index in range(len(data)):
        # take a datapoint every year
        if (index * dt) % (stepsize * 365) == 0:
            point = np.linalg.norm(
                np.array(ref_data)[int((index * dt) / stepsize)] \
                - np.array(data)[index]
            )
            # get difference in percentage
            point = point / np.linalg.norm(np.array(ref_data)[int((index * dt) / stepsize)]) * 100
            plot.append(point)
    return plot

def evaluation(data, ref_data, stepsizes, dt, system):
    """
    Evaluate the simulation data.
    """
    diff_data = []
    find = lambda x: next((i for i in range(len(system)) if system[i].name == x.parent), None)
    for i in range(len(data)):
        body = system[i]
        print("body " + str(i))
        print("name " + str(system[i].name))
        p_index = find(body)
        if p_index != None:
            print("parent " + str(system[p_index].name))
        ref_check = test_ref(data[i], ref_data[i], stepsizes[i], dt)
        diff_data.append(ref_check)
        print()
    return diff_data

def plot_diff(data, methods, t_end):
    # convert t_end from seconds to (Earth) years
    t_end = t_end / (60 * 60 * 24 * 365)

    methods = np.array(methods)

    # get the x values for both datasets
    t_energy = np.linspace(0, t_end, len(data[0][0]))
    t_refcheck = np.linspace(0, t_end, len(data[0][1][0]))

    # copy y axis
    fig, ax1 = plt.subplots()
    if len(methods) > 1: ax1.set_xlim(0, 1.2 * t_end)
    ax2 = ax1.twinx()

    # set ax titles
    ax1.set_xlabel('time (years)')
    ax1.set_ylabel('accuracy error (%)\n for Mercury')
    ax2.set_ylabel('energy error (%)')

    # set plot environment
    plt.title("Accuracy and energy preservation of different integration methods")
    plt.gcf().subplots_adjust(bottom=0.15)
    figManager = plt.get_current_fig_manager()
    figManager.window.showMaximized()
    plt.grid(True)

    lines1 = []
    # plot 1st dataset based on left y axis
    for i in range(len(methods)):
        line, = ax1.plot(
            t_refcheck, data[i][1][1],
            c='C'+str(i % len(methods)), linestyle="-"
        )
        lines1.append(line)

    ax1.legend(
        lines1, methods, title="accuracy", loc=2, bbox_to_anchor=(0.01, 0.993)
    )

    skip = []
    lines2 = []
    # ignore euler and rk2 method because they skew the other results
    if "euler" in methods:
        skip.append(np.where(np.array(methods) == "euler")[0][0])
    if "rk2" in methods:
        skip.append(np.where(np.array(methods) == "rk2")[0][0])
    # plot 2nd dataset based on right y axis
    for i in range(len(methods)):
        if i in skip:
            continue
        line, = ax2.plot(
            t_energy, data[i][0],
            c='C'+str(i % len(methods)), linestyle="-."
        )
        lines2.append(line)

    ax2.legend(
        lines2,
        methods[[i for i in range(len(methods)) if i not in skip]],
        title="energy preservation",
        loc=2,
        bbox_to_anchor=(0.13, 0.993)
    )

    # if there is just one method tested,
    # don't create an inset (not needed)
    if len(methods) == 1:
        plt.show()
        return

    # find the zoom value
    zoom = []
    for i in range(len(methods)):
        if methods[i] == "hermite" or methods[i] == "rk4":
            zoom.append(data[i][1][1][-1])
    ydist = max(zoom) - min(zoom)

    # create an inset which is 0.9 times as long as the y axis
    axins = zoomed_inset_axes(ax1, (1 / ((1.5 * ydist) / max(data[0][1][1]))), loc=1)
    for i in range(len(methods)):
        if methods[i] == "hermite" or methods[i] == "rk4":
            line, = axins.plot(t_refcheck, data[i][1][1], c='C'+str(i % len(methods)), linestyle="-")

    # calculate the x and y limit of the inset
    xdist = 0.0005
    x1, x2 = t_refcheck[-1] - xdist, t_refcheck[-1] + xdist
    y1, y2 = min(zoom) - 0.25 * ydist, max(zoom) + 0.25 * ydist
    axins.set_xlim(x1, x2)
    axins.set_ylim(y1, y2)
    plt.xticks(visible=False)

    # place the inset
    mark_inset(ax1, axins, loc1=2, loc2=3, fc="none", ec="0.5")
    # show plot
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
    methods = ["euler", "verlet", "rk2", "rk4", "hermite"]
    if args.itype.lower() != "all":
        methods = [args.itype.lower()]

    method_data = []
    ref_data = []
    stepsizes = []
    # run simulation with the preferred method(s) & save the data
    for method in methods:
        with open(".tmp.input", 'r+') as f:
            # set up system
            G, sys = get_system_data(f, True)
            system = System(G, sys, method)
            # import all reference data
            if method == methods[0]:
                ref_data, stepsizes = read_all_refs(system.sys)
            # run the simulation
            print("itype " + str(system.itype) + "\n")
            sim_data = simulate(system, args.t_end, args.dt, args.t_dia,
                                args.t_out, args.verbose)
            # evaluate: compare with reference data
            ref_diff = evaluation(sim_data, ref_data, stepsizes, args.dt, system.sys)
            method_data.append([system.e_diff, ref_diff])

    # delete the extra input file
    os.remove(".tmp.input")

    # generate plot
    plot_diff(method_data, methods, args.t_end)

    if args.plot_2d or args.plot_3d:
        simple_plot([p.T for p in sim_data], args.plot_3d, args.n_points)
