import sys as syspy
import argparse
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import system_interpolation as cmod


def simple_plot(data, plt_3d=False, n_points=-1):
    if n_points >= 0 and n_points < len(data[0][0]):
        data = [[dim[::(len(dim)//n_points+1)] for dim in b] for b in data]
    fig = plt.figure()
    if plt_3d:
        ax = fig.add_subplot(111, projection="3d")
        for i in range(len(data)):
            ax.scatter(data[i][0], data[i][1], data[i][2], c='C'+str(i%9), s=1)
            ax.scatter(
                data[i][0][-1], data[i][1][-1], data[i][2][-1],
                c='C'+str(i%9), s=17.5
            )
    else:
        ax = fig.add_subplot(111)
        for i in range(len(data)):
            ax.scatter(data[i][0], data[i][1], c='C'+str(i%9), s=1)
            ax.scatter(data[i][0][-1], data[i][1][-1], c='C'+str(i%9), s=17.5)
    plt.grid(True)
    plt.axis("equal")
    plt.show()

class Body:
    def __init__(self, m, p, v, names):
        self.dim = len(p)
        self.m = float(m)
        self.p = np.array(p, dtype=float)
        self.v = np.array(v, dtype=float)
        if len(names) == 2:
            self.name = names[0].lower()
            self.parent = names[1].lower()
    def to_string(self):
        s = ""
        s += str(self.m) + "\n"
        s += " ".join([str(x) for x in self.p]) + "\n"
        s += " ".join([str(x) for x in self.v])
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
        self.e_diff = []
    def to_string(self):
        s = ""
        for b in self.sys:
            s += b.to_string() + "\n"
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
    def print_energy(self, t, n, verbose=True):
        ek = self.get_ek()
        ep = self.get_ep()
        etot = ek + ep
        if verbose:
            syspy.stdout.write(
                "at time t = " + str(t) + ", after " + str(n) + " steps :\n"
            )
            syspy.stdout.write(
                "  E_kin = " + str(ek) + ", E_pot = " + str(ep) + \
                ", E_tot = " + str(etot) + "\n"
            )
            syspy.stdout.write(
                "             E_tot - E_init = " + str(etot-self.e0) + "\n"
            )
            syspy.stdout.write(
                "  (E_tot - E_init) / E_init =  " + \
                str((etot - self.e0) / self.e0) + "\n"
            )
            self.e_diff.append((etot - self.e0) / self.e0)
        else:
            self.e_diff.append((etot - self.e0) / self.e0)

def simulate(system, t_end, dt, dt_dia=-1, dt_out=-1, verbose=True):
    n = 0
    t = 0
    t_out = dt_out - 0.5 * dt
    t_dia = dt_dia - 0.5 * dt
    if dt_out > 0:
        syspy.stdout.write(str(t) + "\n")
        for l in system.to_string().split("\n"):
            syspy.stdout.write(l + "\n")
    if dt_dia > 0:
        system.print_energy(t, n, verbose)
    lst = [[x.p for x in system.sys]]
    while t < t_end:
        system.step(dt)
        t += dt
        n += 1
        if dt_dia > 0 and t >= t_dia:
            system.print_energy(t, n, verbose)
            t_dia += dt_dia
        if dt_out > 0 and t >= t_out:
            syspy.stdout.write(str(t) + "\n")
            for l in system.to_string().split("\n"):
                syspy.stdout.write(l + "\n")
            t_out += dt_out
        lst.append([b.p for b in system.sys])
    return [
        p.squeeze() for p in np.split(np.array(lst), len(system.sys), axis=1)
    ]

def set_parser(parser):
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
    parser.add_argument(
        "-np", "--n_points", type=int, default=-1,
        help="maximum amount of points plotted"
    )

def get_system_data(f=None, names=False):
    input_dim = -1
    b_data = []
    if not f:
        f = syspy.stdin
    G = float(f.readline().strip())
    for line in f:
        if line[0] == "#":
            continue
        line_strip = line.strip()
        if line_strip == '':
            break
        b_params = [float(x) for x in line_strip.split()[:7]]
        b_dim, r = divmod(len(b_params[1:]), 2)
        assert r == 0, "Inconsistent input dimensions"
        if input_dim < 0: input_dim = b_dim
        else: assert input_dim == b_dim, "Inconsistent input dimensions"
        if names:
            assert len(line_strip.split()) == 9, "Inconsistent input dimensions"
            b_names = [str(x) for x in line_strip.split()[7:]]
        else: b_names = []
        b_data.append(
            Body(
                b_params[0],
                b_params[1:b_dim+1],
                b_params[b_dim+1:2*b_dim+1],
                b_names
            )
        )
    return G, b_data

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    set_parser(parser)
    args = parser.parse_args()
    G, sys = get_system_data()
    system = System(G, sys, args.itype.lower())
    sim_data = simulate(system, args.t_end, args.dt, args.t_dia, args.t_out)
    if args.plot_2d or args.plot_3d:
        simple_plot(
            [np.array(p).T for p in sim_data], args.plot_3d, args.n_points
        )

