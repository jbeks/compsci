import sys
import numpy as np
import matplotlib.pyplot as plt
from nbody import *

def v_error(data, index):
    return np.linalg.norm(data[0]) - np.linalg.norm(data[index])

def p_error(data, center, orbit_i):
    change = data[orbit_i] - data[0]
    center_change = center[orbit_i] - center[0]
    return np.linalg.norm(change) - np.linalg.norm(center_change)

def v_change(data):
    norm = []
    for point in data:
        norm.append(np.linalg.norm(point))
    return (norm - np.roll(norm, -1))[:-1]

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
    while orbit_i < len(data) - 1 and orbit_i != -1:
        orbit_err.append(orbit_i)
        orbit_i = find_period(data, center, orbit_i)
    orbit_diff = (np.roll(orbit_err, -1) - orbit_err)[:-1]
    return (orbit_err[1] if len(orbit_err) > 1 else -1), orbit_err, orbit_diff

def orbit_type(a):
    same = True
    for i in range(a.size - 1):
        if a[i + 1] < a[i]:
                return "cyclic"
        if same and a[i + 1] != a[i]:
            same = False
    return "equal" if same else "divergent"

def make_plot(data, dt_i):
    # print(ap)
    for i in range(len(data)):
        x = np.linspace(0, 176, len(data[i]))
        plt.plot(x, data[i], label=str("line: " + str(dt_i[i])), c='C'+str(i%9))
        # t = np.linspace(0, 176, len(ap[i]))
        # plt.scatter([(x[0] * dt_i[i]) / (60 * 60 * 24) for x in ap[i]], t, c='C'+str(i%9))
        # plt.scatter([(x[1] * dt_i[i]) / (60 * 60 * 24) for x in ap[i]], t, c='C'+str(i%9))
    plt.legend(['dt = 100', 'dt = 300', 'dt = 500', 'dt = 700', 'dt = 1000', 'dt = 2000', 'dt = 3000', 'dt = 4000'])
    plt.grid(True)
    plt.show()

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

def test_ref(data, v_data, body, dt):
    if body != "Phobos":
        return []
    ref_data, stepsize = read_ref(body.lower())
    plot = []
    if body == "Phobos":
        for index in range(len(data)):
            if (index * dt) % stepsize == 0:
                print("DAY: " + str(int((index * dt) / (60 * 60 * 24))))
                print(ref_data[int((index * dt) / stepsize)])
                print(data[index])
                point = np.linalg.norm(np.array(ref_data)[int((index * dt) / stepsize)] - np.array(data)[index])
                print(point)
                plot.append(point)
                print()
    return plot

def apoperi(data_body, data_center, orbit_i):
    print(orbit_i)
    dists = []
    for start, end in zip(orbit_i[:-1], orbit_i[1:]):
        min_dist = np.linalg.norm(data_body[start] - data_center[start])
        max_dist = 0
        min_i = max_i = -1
        for i, (d, d0) in enumerate(zip(data_body[start:end], data_center[start:end])):
            dist = np.linalg.norm(d - d0)
            if dist < min_dist:
                min_dist = dist
                min_i = i
            if dist > max_dist:
                max_dist = dist
                max_i = i
        dists.append([min_i, max_i])
    return dists


def evaluation(data, v_data, dt, system):
    find = lambda x: next((i for i in range(len(system)) if system[i].name == x.parent), None)
    ref_check_m = []
    point_m = []
    for i in range(len(data)):
        body = system[i]
        print("b" + str(i))
        print("name " + str(system[i]))
        p_index = find(body)
        if p_index != None:
            print("parent " + str(system[p_index]))
            mi, ma = find_dists(data[i], data[p_index])
            print("periapsis " + str(mi))
            print("apoapsis " + str(ma))
            orbit_i, points, diff = orbit_data(data[i], data[p_index])
            print("orbit " + str(len(points) - 1))
        #     if len(points) - 1 > 0:
        #         print("avg " + str(np.mean(diff) * dt / 86400.))
        #         print("max " + str(max(diff) * dt / 86400.))
        #         print("min " + str(min(diff) * dt / 86400.))
        #         print("type " + orbit_type(diff))
        #         p_err = p_error(data[i], data[p_index], orbit_i)
        #         print("p_err " + str(p_err))
        #         v_err = v_error(v_data[i], orbit_i)
        #         print("v_err " + str(v_err))
        # else:
        #     v_ch = v_change(v_data[i])
        #     print("v_change ")
        #     print(v_ch)
        ref_check = test_ref(data[i], v_data[i], body.name, dt)
        # print("ref-check " + str(ref_check))
        # print()
        if ref_check != []:
            ref_check_m = ref_check
        if i == 5:
            point_m = points
    return ref_check_m, point_m

if __name__ == "__main__":
    orbit = {'day': 1, 'Mercury': 87.97, 'Venus': 224.70, 'Earth': 365.25, 'Moon': 27.32, 'Mars': 686.98, 'Phobos': 0.3189}
    time = [100, 300, 500, 700, 1000, 2000, 3000, 4000] #, 6000, 8000, 10000]
    itype, t_end, dt, t_dia, t_out, plot, save, planet, years = parse_arguments()
    if save or plot:
        print("Sorry, we don't do saving and/or plotting...")
    if t_end == 1:
        t_end = int(years) * orbit[planet] * 24 * 60 * 60
    dt_i = time[:int(dt)]
    plots = []
    ap = []
    G, sys = get_system_data()
    for dt in dt_i:
        system = System(G, sys, itype)
        print("itype " + itype)
        print("b " + str(len(system.sys)))
        print("dt " + str(dt))
        print("t_end " + str(t_end))
        print("E_init " + str(system.e0))
        print()
        sim_data, v_data = simulate(system, t_end, dt, t_dia, t_out)
        dt_data, orbit_i = evaluation(sim_data, v_data, dt, system.sys)
        # ap.append(apoperi(sim_data[1], sim_data[0], orbit_i))
        plots.append(dt_data)
        G, sys = get_system_data(True)
    make_plot(plots, dt_i)
    if save:
        f = open(save + ".npy", 'wb')
        np.save(f, sim_data)
    if plot:
        #simple_plot([p.T for p in sim_data])
        full_plot([p.T for p in sim_data])
