import sys
import numpy as np

def read_sim_data():
    t = []
    m = []
    p = []
    v = []
    i = 0
    j = 0
    for line in sys.stdin:
        if i == 0:
            m.append(float(line.strip()))
        elif i == 1:
            lst = [float(x) for x in line.strip().split()]
            if len(lst) == 1:
                t.append(m[-1])
                m[-1] = lst[0]
                i-= 1
                j = 0
            else:
                if len(t) == 1:
                    p.append([])
                p[j].append(lst)
        elif i == 2:
            if len(t) == 1:
                v.append([])
            v[j].append([float(x) for x in line.strip().split()])
            i = -1
            j += 1
        i += 1
    return t, m[:len(m)//len(t)], p, v

if __name__ == "__main__":
    for x in read_sim_data():
        print(np.array(x).shape)

