import sys
from read_sim_data import read_sim_data

def shorten_sim_data():
    t, m, p, _ = read_sim_data()
    s = "\n".join([str(mass) for mass in m])
    for i in range(len(t)):
        s += "\n" + str(t[i])
        for b in p:
            s += "\n" + " ".join([str(dim) for dim in b[i]])
    return s

if __name__ == "__main__":
    short_data = shorten_sim_data()
    sys.stdout.write(short_data)

