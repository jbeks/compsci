import sys
from read_sim_data import read_sim_data

def shorten_sim_data():
    """
    Shorten normal format simulation data
    to short format simulation data
    """
    # read data from file
    t, m, p, _ = read_sim_data()
    # start file by listing all masses
    s = "\n".join([str(mass) for mass in m])
    # add time and positions at that time
    for i in range(len(t)):
        s += "\n" + str(t[i])
        for b in p:
            s += "\n" + " ".join([str(dim) for dim in b[i]])
    return s

if __name__ == "__main__":
    short_data = shorten_sim_data()
    sys.stdout.write(short_data)

