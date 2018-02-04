import sys
import numpy as np

def read_sim_data():
    """
    Reads simulation data in normal format from standard input
    Returns read data (time, mass, position, velocity)
    """
    # initialize variables
    t = []
    m = []
    p = []
    v = []
    i = 0
    j = 0
    for line in sys.stdin:
        # record mass
        if i == 0:
            m.append(float(line.strip()))
        # record position
        elif i == 1:
            lst = [float(x) for x in line.strip().split()]
            # if a single value was found,
            # previous value was a time value
            if len(lst) == 1:
                # last recorded value in "m" was a time
                t.append(m[-1])
                m[-1] = lst[0]
                i-= 1
                j = 0
            else:
                if len(t) == 1:
                    p.append([])
                p[j].append(lst)
        # record velocity
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

