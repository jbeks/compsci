import sys
import numpy as np
import matplotlib.pyplot as plt

def simple_plot(xs, ys, type=0):
    for i in range(len(xs)):
        if type == 0:
            plt.scatter(xs[i], ys[i], c='C'+str(i%9), s=10)
        elif type ==1:
            plt.plot(xs[i], ys[i], c='C'+str(i%9), lw=2)
    plt.axvline(0, c='k')
    plt.axhline(0, c='k')
    plt.grid(True)
    plt.axis("equal")
#    plt.ylim([-2,12])
#    plt.xlim([-2,12])
    plt.show()

class Integration:
    def euler(dt, p0, v0, afunc):
        a0 = afunc(p0)
        p1 = p0 + v0 * dt
        v1 = v0 + a0 * dt
        return p1, v1
    def verlet(dt, p0, v0, afunc):
        a0 = afunc(p0)
        p1 = p0 + v0 * dt + 0.5 * a0 * dt * dt
        a1 = afunc(p1)
        v1 = v0 + 0.5 * (a0 + a1) * dt
        return p1, v1
    def rk2(dt, p0, v0, afunc):
        a0 = afunc(p0)
        v0_5 = v0 + 0.5 * a0 * dt
        p0_5 = p0 + 0.5 * v0 * dt
        a0_5 = afunc(p0_5)
        v1 = v0 + a0_5 * dt
        p1 = p0 + v0_5 * dt
        return p1, v1
    def rk4(dt, p0, v0, afunc):
        k1 = afunc(p0) * dt
        k2 = afunc(p0 + 0.5 * v0 * dt + 0.125 * k1 * dt) * dt
        k3 = afunc(p0 + v0 * dt + 0.5 * k2 * dt) * dt
        p1 = p0 + v0 * dt + (1 / 6.) * (k1 + 2 * k2) * dt
        v1 = v0 + (1 / 6.) * (k1 + 4 * k2 + k3)
        return p1, v1

class Body:
    def __init__(self, p, m, v):
        assert len(p) == 3, "Insufficient location parameters"
        assert len(v) == 3, "Insufficient velocity parameters"
        self.m = float(m)
        self.p = np.array(p)
        self.v = np.array(v)
        self.p_tmp = np.zeros(3)
        self.v_tmp = np.zeros(3)
        self.e0 = 0
    def __repr__(self):
        return "({}, {}, {})".format(self.m, self.p, self.v)
    def set_update(self, p, v):
        self.p_tmp = p
        self.v_tmp = v
    def confirm_update(self):
        self.p = self.p_tmp
        self.v = self.v_tmp
    def set_e0(self, sys, G):
        self.e0 = self.get_ek(sys) + self.get_ep(sys, G)
    def get_ek(self, sys):
        return .5 * self.m * np.linalg.norm(self.v) ** 2
    def get_ep(self, sys, G):
        ep = 0
        for b in sys:
            if b == self:
                continue
            vector = b.p - self.p
            ep -= b.m / np.linalg.norm(vector)
        return ep * self.m * G

class System:
    def __init__(self, sys, G, ifunc):
        self.G = G
        self.sys = sys
        self.ifunc = ifunc
        for b in sys:
            b.set_e0(sys, G)
    def __repr__(self):
        s = ''
        for b in self.sys:
            s += str(b) + '\n'
        return s[:-1]
    # make sure sys does not contain the body with position p
    def get_acceleration(self, p, sys, G):
        a = 0
        for b in sys:
            vector = b.p - p
            a += vector * G * b.m / float(np.linalg.norm(vector) ** 3)
        return a
    def step(self, dt):
        for b in self.sys:
            b.set_update(
                *self.ifunc(
                    dt,
                    b.p,
                    b.v,
                    lambda p: self.get_acceleration(
                        p, [x for x in self.sys if x != b], self.G
                    )
                )
            )
            break                                                       #TODO
        for b in self.sys:
            b.confirm_update()
            break                                                       #TODO
    def print_energy(self, t, n):
        for b in self.sys:
            ek = b.get_ek(self.sys)
            ep = b.get_ep(self.sys, self.G)
            etot = ek + ep
            print('at time t =',t,', after',n,'steps :')
            print('  E_kin =',ek,', E_pot =',ep,', E_tot =',etot)
            print('             E_tot - E_init = ',etot-b.e0)
            print('  (E_tot - E_init) / E_init =', (etot - b.e0) / b.e0)
#            print('t          ', t)
#            print('ek         ', ek)
#            print('ep         ', ep)
#            print('etot       ', etot)
#            print('e0         ', b.e0)
#            print('etot - e0  ', etot - b.e0)
#            print('etot-e0/e0 ', (etot - b.e0) / b.e0)
            break                                                       #TODO

def simulate(sys, dt, t_end, dt_out=-1, dt_dia=-1, show_plot=False):
    n = 0
    t =  .5 * dt
    t_out = dt_out
    t_dia = dt_dia

    if dt_dia > 0:
        sys.print_energy(0, n)

    lst = [[x.p for x in sys.sys]]
    while t < t_end:

        sys.step(dt)
        t += dt
        n += 1

        if dt_dia > 0 and t >= t_dia:
            sys.print_energy(t_dia, n)
            t_dia += dt_dia

        if dt_out > 0 and t >= t_out:
            print(sys)
            t_out += dt_out

        lst.append([x.p for x in sys.sys])

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
                Body(b_params[1:4], b_params[0], b_params[4:])
            )
        line_index += 1
    return params

if __name__ == "__main__":
    assert len(sys.argv) == 2, "Insufficient amount of arguments"
    itype = sys.argv[1].lower()
    assert itype in [
        "euler",
        "verlet",
        "leapfrog",
        "rk2",
        "rk4",
    ], "Integration type not supported"
    if itype == "euler":
        ifunc = Integration.euler
    elif itype == "verlet" or itype == "leapfrog":
        ifunc = Integration.verlet
    elif itype == "rk2":
        ifunc = Integration.rk2
    elif itype == "rk4":
        ifunc = Integration.rk4
    params = get_params()
    system = System(*params[:2], ifunc)
    simulate(system, *params[2:], False)

