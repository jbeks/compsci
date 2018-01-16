import sys
import numpy as np
import matplotlib.pyplot as plt

def simple_plot(xs, ys, type=0):
    for i in range(len(xs)):
        if type == 0:
            plt.scatter(xs[i], ys[i], c='C'+str(i%9), s=10)
        elif type ==1:
            plt.plot(xs[i], ys[i], c='C'+str(i%9), lw=2)
    plt.grid(True)
#    plt.axvline(0, c='k')
#    plt.axhline(0, c='k')
#    plt.axis("equal")
#    plt.ylim([-2,12])
#    plt.xlim([-2,12])
    plt.show()

class Body:
    def __init__(self, m, p, v):
        assert len(p) == 3, "Insufficient location parameters"
        assert len(v) == 3, "Insufficient velocity parameters"
        self.m = float(m)
        self.p = np.array(p)
        self.v = np.array(v)
        self.a = np.zeros(3)
        self.j = np.zeros(3)
#        self.p_tmp = np.zeros(3)
#        self.v_tmp = np.zeros(3)
    def __repr__(self):
        return "({}, {}, {})".format(self.m, self.p, self.v)
#    def set_update(self, p, v):
#        self.p_tmp = p
#        self.v_tmp = v
#    def confirm_update(self):
#        self.p = self.p_tmp
#        self.v = self.v_tmp
    def get_ek(self):
        # TODO removed * self.m
        return .5 * np.linalg.norm(self.v) ** 2
    def get_ep(self, sys, G):
        ep = 0
        for b in sys:
            if b == self:
                continue
            vector = b.p - self.p
            ep -= b.m / np.linalg.norm(vector)
        return ep * self.m * G
    def set_p(self, p):
        self.p = p
    def set_v(self, v):
        self.v = v
    def calc_a(self, sys, G):
        self.a = np.zeros(3)
        for b in sys:
            if b == self:
                continue
            vector = b.p - self.p
            self.a += vector * b.m / float(np.linalg.norm(vector) ** 3)
        self.a *= G
    def calc_j(self, sys, G):
        self.j = np.zeros(3)
        for b in sys:
            if b == self:
                continue
            p_vector = b.p - self.p
            v_vector = b.v - self.v
            len_p_vector = np.linalg.norm(p_vector)
            self.j += b.m * (v_vector / float(len_p_vector ** 3)
                - 3 * np.dot(p_vector, v_vector) * p_vector / len_p_vector ** 5)
        self.j *= G

class System:
    def __init__(self, sys, G, itype):
        self.G = G
        self.sys = sys
        self.e0 = self.get_ek() + self.get_ep()
        self.ifunc = None
        assert itype in [
            "euler",
            "verlet",
            "leapfrog",
            "rk2",
            "rk4",
            "hermite",
        ], "Integration type not supported"
        if itype == "euler":
            self.ifunc = self.euler
        elif itype == "verlet" or itype == "leapfrog":
            self.ifunc = self.verlet
        elif itype == "rk2":
            self.ifunc = self.rk2
        elif itype == "rk4":
            self.ifunc = self.rk4
        elif itype == "hermite":
            self.ifunc = self.hermite
    def __repr__(self):
        s = ''
        for b in self.sys:
            s += str(b) + '\n'
        return s[:-1]
    # make sure sys does not contain the body with position p
    def step(self, dt):
        self.ifunc(dt)
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
    def euler(self, dt):
        for b in self.sys: b.calc_a(self.sys, self.G); break
        for b in self.sys: b.set_p(b.p + b.v * dt); break
        for b in self.sys: b.set_v(b.v + b.a * dt); break
    def verlet(self, dt):
        for b in self.sys: b.calc_a(self.sys, self.G); break
        for b in self.sys: b.set_v(b.v + .5 * b.a * dt); break
        for b in self.sys: b.set_p(b.p + b.v * dt); break
        for b in self.sys: b.calc_a(self.sys, self.G); break
        for b in self.sys: b.set_v(b.v + .5 * b.a * dt); break
    def rk2(self, dt):
        p0 = []
        v0_5 = []
        for b in self.sys: b.calc_a(self.sys, self.G); break
        for b in self.sys:
            p0.append(b.p)
            v0_5.append(b.v + .5 * b.a * dt)
            break
        for b in self.sys: b.set_p(b.p + .5 * b.v * dt); break
        for b in self.sys: b.calc_a(self.sys, self.G); break
        for b in self.sys: b.set_v(b.v + b.a * dt); break
        for i in range(len(self.sys)):
            self.sys[i].set_p(p0[i] + v0_5[i] * dt)
            break
    def rk4(self, dt):
        p0 = []
        v0 = []
        k1 = []
        k2 = []
        k3 = []
        for b in self.sys:
            p0.append(b.p)
            v0.append(b.v)
            break
        for b in self.sys:
            b.calc_a(self.sys, self.G)
            k1.append(b.a * dt)
            break
        for i in range(len(self.sys)):
            self.sys[i].set_p(p0[i] + .5 * v0[i] * dt + .125 * k1[i] * dt)
            self.sys[i].calc_a(self.sys, self.G)
            k2.append(b.a * dt)
            break
        for i in range(len(self.sys)):
            self.sys[i].set_p(p0[i] + v0[i] * dt + .5 * k2[i] * dt)
            self.sys[i].calc_a(self.sys, self.G)
            k3.append(b.a * dt)
            break
        for i in range(len(self.sys)):
            self.sys[i].set_p(p0[i] + v0[i] * dt + (k1[i] + 2*k2[i]) * dt / 6.)
            self.sys[i].set_v(v0[i] + (k1[i] + 4*k2[i] + k3[i]) / 6.)
            break
    def hermite(self, dt):
        p0 = []
        v0 = []
        a0 = []
        j0 = []
        for b in self.sys:
            p0.append(b.p)
            v0.append(b.v)
            break
        for b in self.sys:
            b.calc_a(self.sys, self.G)
            b.calc_j(self.sys, self.G)
            a0.append(b.a)
            j0.append(b.j)
            break
        for b in self.sys:
            b.set_p(b.p + b.v * dt + .5 * b.a * dt*dt + b.j * dt*dt*dt / 6.)
            b.set_v(b.v + b.a * dt + .5 * b.j * dt*dt)
            break
        for b in self.sys:
            b.calc_a(self.sys, self.G)
            b.calc_j(self.sys, self.G)
            break
        for i in range(len(self.sys)):
            self.sys[i].set_v(v0[i] + .5 * (a0[i] + self.sys[i].a) * dt \
                + (j0[i] - self.sys[i].j) * dt*dt / 12.)
            break
        for i in range(len(self.sys)):
            self.sys[i].set_p(p0[i] + .5 * (v0[i] + self.sys[i].v) * dt \
                + (a0[i] - self.sys[i].a) * dt*dt / 12.)
            break

def simulate(system, dt, t_end, dt_out=-1, dt_dia=-1, show_plot=False):
    n = 0
    t =  .5 * dt
    t_out = dt_out
    t_dia = dt_dia

    if dt_dia > 0:
        system.print_energy(0, n)

    lst = [[x.p for x in system.sys]]
    while t < t_end:

        system.step(dt)
        t += dt
        n += 1

        if dt_dia > 0 and t >= t_dia:
            system.print_energy(t_dia, n)
            t_dia += dt_dia

        if dt_out > 0 and t >= t_out:
            print(system)
            t_out += dt_out

        lst.append([x.p for x in system.sys])

    lst = np.array(lst).T

    if show_plot:
        simple_plot(lst[0], lst[1], 0)

    return None

def get_params():
    line_index = 0
    params = [[]]
    for line in sys.stdin:
        line_strip = line.strip()
        if line_strip == '':
            break
        if line_index < 5:
            params.append(float(line_strip))
        else:
            b_params = [float(x) for x in line_strip.split(' ')]
            params[0].append(
                Body(b_params[0], b_params[1:4], b_params[4:])
            )
        line_index += 1
    return params

if __name__ == "__main__":
    assert len(sys.argv) == 2, "Insufficient amount of arguments"
    itype = sys.argv[1].lower()
#    assert itype in [
#        "euler",
#        "verlet",
#        "leapfrog",
#        "rk2",
#        "rk4",
#        "hermite",
#    ], "Integration type not supported"
#    if itype == "euler":
#        ifunc = Integration.euler
#    elif itype == "verlet" or itype == "leapfrog":
#        ifunc = Integration.verlet
#    elif itype == "rk2":
#        ifunc = Integration.rk2
#    elif itype == "rk4":
#        ifunc = Integration.rk4
#    elif itype == "hermite":
#        ifunc = Integration.hermite
    params = get_params()
    system = System(*params[:2], itype)
    simulate(system, *params[2:], False)

