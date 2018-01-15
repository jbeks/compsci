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

class Body:
    def __init__(self, p, m, v):
        assert len(p) == 3, "Insufficient location parameters"
        assert len(v) == 3, "Insufficient velocity parameters"
        self.p = np.array(p)
        self.m = float(m)
        self.v = np.array(v)
        self.p_tmp = np.zeros(3)
        self.e0 = 0
    def __repr__(self):
        return "({}, {}, {})".format(self.m, self.p, self.v)
    def acceleration(self, sys, G):
        a = 0
        for b in sys:
            if b == self:
                continue
            vector = b.p - self.p
            a += vector * G * b.m / float(np.linalg.norm(vector) ** 3)
        return a
    def calculate_update(self, sys, dt, G):
        assert False, "Implement in child class"
    def confirm_update(self):
        self.p = self.p_tmp
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
    def __init__(self, sys, G):
        self.G = G
        self.sys = sys
        for b in sys:
            b.set_e0(sys, G)
    def __repr__(self):
        s = ''
        for b in self.sys:
            s += str(b) + '\n'
        return s[:-1]
    def step(self, dt):
        for b in self.sys:
            b.calculate_update(self.sys, dt, self.G)
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

class Body_Euler(Body):
    def calculate_update(self, sys, dt, G):
        a = self.acceleration(sys, G)
        self.p_tmp = self.p + self.v * dt
        self.v = self.v + a * dt

class Body_Verlet(Body):
    def calculate_update(self, sys, dt, G):
        p0 = self.p
        v0 = self.v
        a0 = self.acceleration(sys, G)
        self.v = self.v + 0.5 * a0 * dt
        self.p = self.p + self.v * dt
        a1 = self.acceleration(sys, G)
        self.v = v0 + 0.5 * (a0 + a1) * dt
        self.p_tmp = self.p
        self.p = p0

# TODO next Bodies copied from artcopsci, change to fit code more
#      might not work when updating multiple bodies

class Body_Leapfrog(Body):
    def calculate_update(self, sys, dt, G):
        p0 = self.p
        a0 = self.acceleration(sys, G)
        self.p = self.p + self.v * dt + .5 * a0 * dt * dt
        a1 = self.acceleration(sys, G)
        self.v = self.v + .5 * (a0 + a1) * dt
        self.p_tmp = self.p
        self.p = p0

class Body_RK2(Body):
    def calculate_update(self, sys, dt, G):
        p0 = self.p
        v0 = self.v + .5 * a0 * dt
        a0 = self.acceleration(sys, G)
        self.p = self.p + .5 * self.v * dt
        a1 = self.acceleration(sys, G)
        self.v = self.v + a1 * dt
        self.p = p0 + v0 * dt
        self.p_tmp = self.p
        self.p = p0

class Body_RK4(Body):
    def calculate_update(self, sys, dt, G):
        p0 = self.p
        a0 = self.acceleration(sys, G)
        self.p = p0 + .5 * self.v * dt + .125 * a0 * dt * dt
        a1 = self.acceleration(sys, G)
        self.p = p0 + self.v * dt + 0.5 * a1 * dt * dt
        a2 = self.acceleration(sys, G)
        self.p = p0 + self.v * dt + (1 / 6.) * (a0 + 2 * a1) * dt * dt
        self.v = self.v + (1 / 6.) * (a0 + 4 * a1 + a2) * dt
        self.p_tmp = self.p
        self.p = p0

#class System_Verlet(System):
#    def step(self, dt):
#        for _ in range(2):
#            super().step(dt)

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

def get_params(body_type):
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
                body_type(b_params[1:4], b_params[0], b_params[4:])
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
#    if itype == "verlet":
#        params = get_params(Body_Verlet)
#        system = System_Verlet(*params[:2])
#    else:
    if itype == "euler":
        btype = Body_Euler
    elif itype == "verlet":
        btype = Body_Verlet
    elif itype == "leapfrog":
        btype = Body_Leapfrog
    elif itype == "rk2":
        btype = Body_RK2
    elif itype == "rk4":
        btype = Body_RK4
    params = get_params(btype)
    system = System(*params[:2])
    simulate(system, *params[2:], False)

