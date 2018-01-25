import sys as syspy
import argparse
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import system_interpolation as cmod

def simple_plot(data, plt_3d=False):
    fig = plt.figure()
    if plt_3d:
        ax = fig.add_subplot(111, projection="3d")
        for i in range(len(data)):
            if syspy.version_info[0] < 3:
                ax.scatter(
                    data[i][0], data[i][1], data[i][2], s=1
                )
            else:
                ax.scatter(
                    data[i][0], data[i][1], data[i][2], c='C'+str(i%9), s=1
                )
    else:
        ax = fig.add_subplot(111)
        for i in range(len(data)):
            ax.scatter(data[i][0], data[i][1], c='C'+str(i%9), s=1)
    plt.grid(True)
    plt.axis("equal")
    plt.show()

class Body:
    def __init__(self, m, p, v):
        self.dim = len(p)
        self.m = float(m)
        self.p = np.array(p)
        self.v = np.array(v)
    def __repr__(self):
        s = ""
        s += "mass = " + str(self.m) + "\n"
        s += " pos = " + " ".join([str(x) for x in self.p]) + "\n"
        s += " vel = " + " ".join([str(x) for x in self.v])
        return s
    def get_ek(self):
        return .5 * self.m * np.linalg.norm(self.v) ** 2
    def get_ep(self, sys, G):
        self.p = np.array(self.p)
        ep = 0
        for b in sys:
            if b == self:
                continue
            vector = b.p - self.p
            ep -= b.m / np.linalg.norm(vector)
        return ep * self.m * G

class System:
    def __init__(self, G, sys, itype):
        self.G = G
        self.sys = sys
        self.e0 = self.get_ek() + self.get_ep()
        self.itype = itype
    def __repr__(self):
        s = ""
        for b in self.sys:
            s += str(b) + "\n"
        return s[:-1]
    def step(self, dt):
        sys_list = [[b.m, list(b.p), list(b.v)] for b in self.sys]
        pos, vel = cmod.interpolate(self.itype, dt, self.G, sys_list)
        assert (
            len(pos) == len(self.sys) and len(vel) == len(self.sys)
        ), "Insufficient interpolation output"
        for i in range(len(self.sys)):
            self.sys[i].p = pos[i]
            self.sys[i].v = vel[i]
    def get_ek(self):
        ek = 0
        for b in self.sys:
            ek += b.get_ek()
        return ek
    def get_ep(self):
        ep = 0
        for b in self.sys:
            ep += b.get_ep(self.sys, self.G)
        return 0.5 * ep
    def print_energy(self, t, n):
        ek = self.get_ek()
        ep = self.get_ep()
        etot = ek + ep
        print('at time t =',t,', after',n,'steps :')
        print('  E_kin =',ek,', E_pot =',ep,', E_tot =',etot)
        print('             E_tot - E_init = ',etot-self.e0)
        print('  (E_tot - E_init) / E_init =', (etot - self.e0) / self.e0)

def simulate(system, t_end, dt, dt_dia=-1, dt_out=-1):
    n = 0
    t = 0
    t_out = dt_out - 0.5 * dt
    t_dia = dt_dia - 0.5 * dt
    if dt_out > 0:
        print(system)
    if dt_dia > 0:
        system.print_energy(t, n)
    lst = [[x.p for x in system.sys]]
    while t < t_end:
        system.step(dt)
        t += dt
        n += 1
        if dt_dia > 0 and t >= t_dia:
            system.print_energy(t, n)
            t_dia += dt_dia
        if dt_out > 0 and t >= t_out:
            print(system)
            t_out += dt_out
        lst.append([b.p for b in system.sys])
    return [
        p.squeeze() for p in np.split(np.array(lst), len(system.sys), axis=1)
    ]

def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "itype", type=str,
        help="integration type"
    )
    parser.add_argument(
        "t_end", type=float,
        help="total runtime simulation"
    )
    parser.add_argument(
        "dt", type=float,
        help="timestep for calculations in simulation"
    )
    parser.add_argument(
        "-d", "--t_dia", type=float, default=-1,
        help="timestep for energy evaluation output"
    )
    parser.add_argument(
        "-o", "--t_out", type=float, default=-1,
        help="timestep for body-data output"
    )
    parser.add_argument(
        "-p2d", "--plot_2d", action="store_true",
        help="plot simulation"
    )
    parser.add_argument(
        "-p3d", "--plot_3d", action="store_true",
        help="plot simulation in 3 dimensions"
    )
    args = parser.parse_args()
    return (
        args.itype.lower(), args.t_end, args.dt,
        args.t_dia, args.t_out, args.plot_2d, args.plot_3d
    )

def get_system_data():
    input_dim = -1
    b_data = []
    G = float(syspy.stdin.readline().strip())
    for line in syspy.stdin:
        if line[0] == "#":
            continue
        line_strip = line.strip()
        if line_strip == '':
            break
        b_params = [float(x) for x in line_strip.split()]
        b_dim, r = divmod(len(b_params[1:]), 2)
        assert r == 0, "Inconsistent input dimensions"
        if input_dim < 0: input_dim = b_dim
        else: assert input_dim == b_dim, "Inconsistent input dimensions"
        b_data.append(
            Body(
                b_params[0],
                b_params[1:b_dim+1],
                b_params[b_dim+1:2*b_dim+1]
            )
        )
    return G, b_data

if __name__ == "__main__":
    itype, t_end, dt, t_dia, t_out, p2d, p3d = parse_arguments()
    G, sys = get_system_data()
    system = System(G, sys, itype)
    sim_data = simulate(system, t_end, dt, t_dia, t_out)
    if p2d or p3d:
        simple_plot([p.T for p in sim_data], p3d)

