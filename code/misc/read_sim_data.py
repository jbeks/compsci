import sys

def read_sim_data():
    t = []
    m = []
    p = []
    v = []
    i = 0
    for line in sys.stdin:
        if i == 0:
            m.append(float(line.strip()))
        elif i == 1:
            lst = [float(x) for x in line.strip().split()]
            if len(lst) == 1:
                t.append(m[-1])
                m[-1] = lst[0]
                p.append([])
                v.append([])
                i-= 1
            else:
                p[-1].append(lst)
        elif i == 2:
            v[-1].append([float(x) for x in line.strip().split()])
            i = -1
        i += 1
    print(t)
    return t, m[:len(m)//len(t)], p, v

if __name__ == "__main__":
    for x in read_sim_data():
        print(x)

