import sys
import numpy as np

def read_short_sim_data():
    t = []
    m = []
    p = []
    mass_done = False
    i = 0
    for line in sys.stdin:
        if not mass_done:
            lst = [float(x) for x in line.strip().split()]
            if len(lst) != 1:
                t.append(m[-1])
                m = m[:-1]
                for _ in range(len(m)):
                    p.append([])
                p[i].append(lst)
                i += 1
                mass_done = True
            else:
                m.append(lst[0])
        else:
            if i < len(m):
                p[i].append([float(x) for x in line.strip().split()])
                i += 1
            else:
                t.append(float(line.strip()))
                i = 0
    return t, m, p

if __name__ == "__main__":
    for x in read_short_sim_data():
        print(np.array(x).shape)

