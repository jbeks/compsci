import sys
import numpy as np
import matplotlib.pyplot as plt

def simple_plot(xs, ys, type=0):
    for i in range(len(xs)):
        if type == 0:
            plt.scatter(xs[i], ys[i], c='C'+str(i%9), s=5)
        elif type ==1:
            plt.plot(xs[i], ys[i], c='C'+str(i%9), lw=2)
    plt.axvline(0, c='k')
    plt.axhline(0, c='k')
    plt.grid(True)
    plt.axis("equal")
#    plt.ylim([-2,12])
#    plt.xlim([-2,12])
    plt.show()

class Body_Verlet:
    def __init__(self, p, m, v):
        assert len(p) == 3, "Insufficient location parameters"
        assert len(v) == 3, "Insufficient velocity parameters"
        self.p = np.array(p)
        self.m = float(m)
        self.v = np.array(v)
        self.p_tmp = np.zeros(3)
    def __repr__(self):
        return "({}, {}, {})".format(self.m, self.p, self.v)
    def acceleration(self, sys, G):
        a = 0
        for b in sys:
            if b == self:
                continue
            vector = (b.p - self.p)
            a += vector * G * b.m / float(np.linalg.norm(vector) ** 3)
        return a
    def calculate_update(self, sys, dt, G):
        a = self.acceleration(sys, G)
        self.v = self.v + 0.5 * a * dt
        self.p_tmp = self.p + self.v * dt
    def confirm_update(self):
        self.p = self.p_tmp
# TODO old energy code (might still be useful)
#    def potential_E(self, system, G):
#        result = 0
#        for i in range(len(system)):
#            for j in range(i):
#                result += G* system[i].m * system[j].m / \
#                    np.linalg.norm(system[j].p - system[i].p)
#        return result
#    def kinetic_E(self, system):
#        result = 0
#        for i in range(len(system)):
#            result += np.linalg.norm(system[i].m * system[i].v) ** 2 / \
#                float(2 * system[i].m)
#        return result

class System_Verlet:
    def __init__(self, sys, G):
        self.G = G
        self.sys = sys
    def __repr__(self):
        s = ''
        for b in self.sys:
            s += str(b) + '\n'
        return s[:-1]
    def step(self, dt):
        for _ in range(2):
            for b in self.sys:
                b.calculate_update(self.sys, dt, self.G)
            for b in self.sys:
                b.confirm_update()

def verlet(system, dt, t, pr=False):
# TODO old energy code (might still be useful)
#    Eko = .5*system[0].m*np.linalg.norm(system[0].v)*np.linalg.norm(system[0].v)
#    Eka = .5*system[1].m*np.linalg.norm(system[1].v)*np.linalg.norm(system[1].v)
#    Epo = system[0].m*G*np.linalg.norm(system[1].p - system[0].p)
#    Epa = system[1].m*G*np.linalg.norm(system[0].p - system[1].p)

    lst = []
    for i in range(int(t/float(dt))+1):
# TODO old energy code (might still be useful)
#        for b in system:
#            p = b.potential_E(system, G)
#            k = b.kinetic_E(system)
#            print(p, k, p+k)

        if pr:
            print(i*dt)
            print(system)

        lst.append([x.p for x in system.sys])
        system.step(dt)

# TODO old energy code (might still be useful)
#        Eko_tmp = .5*system[0].m*np.linalg.norm(system[0].v)*np.linalg.norm(system[0].v)
#        Eka_tmp = .5*system[1].m*np.linalg.norm(system[1].v)*np.linalg.norm(system[1].v)
#        Epo_tmp = system[0].m*G*np.linalg.norm(system[1].p - system[0].p)
#        Epa_tmp = system[1].m*G*np.linalg.norm(system[0].p - system[1].p)
#        print(-Eko+Eko_tmp-Epo+Epo_tmp)
#        print(-Eka+Eka_tmp-Epa+Epa_tmp)
#        Eko = Eko_tmp
#        Eka = Eka_tmp
#        Epo = Epo_tmp
#        Epa = Epa_tmp

    lst = np.array(lst).T
    simple_plot(lst[0], lst[1], 0)
    return None

def get_params(body_type):
    line_index = 0
    params = [[]]
    for line in sys.stdin:
        line_strip = line.strip()
        if line_strip == '':
            break
        if line_index < 3:
            params.append(float(line_strip))
        else:
            b_params = [float(x) for x in line_strip.split(' ')]
            params[0].append(
                body_type(b_params[1:4], b_params[0], b_params[4:])
            )
        line_index += 1
    return params

if __name__ == "__main__":
    params = get_params(Body_Verlet)
    system = System_Verlet(*params[:2])
    verlet(system, *params[2:], False)


