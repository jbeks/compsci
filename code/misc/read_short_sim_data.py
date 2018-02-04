import sys
import numpy as np

# reads simulation data in short format from standard input
# returns read data (time, mass, position)
def read_short_sim_data(fname=None):
    # initialize variables
    t = []
    m = []
    p = []
    mass_done = False
    i = 0
    # use file if filename is given, else use standard input
    if fname:
        f = open(fname)
    else:
        f = sys.stdin
    for line in f:
        # first record all mass values
        if not mass_done:
            lst = [float(x) for x in line.strip().split()]
            # if more than one value is found,
            # this is assumed to be the first position of the data
            if len(lst) != 1:
                # last recorded value was first time value
                t.append(m[-1])
                m = m[:-1]
                # record position
                for _ in range(len(m)):
                    p.append([])
                p[i].append(lst)
                i += 1
                # signal all mass has been recorded
                mass_done = True
            else:
                m.append(lst[0])
        else:
            # record position for all bodies
            if i < len(m):
                p[i].append([float(x) for x in line.strip().split()])
                i += 1
            # record  time value
            else:
                t.append(float(line.strip()))
                i = 0
    # close file if file was used
    if fname:
        f.close()
    return t, m, p

if __name__ == "__main__":
    for x in read_short_sim_data():
        print(np.array(x).shape)

