import os
import numpy as np
import matplotlib.pyplot as plt

def read_res(resdir):
    sep = "\\" if os.name == "nt" else "/"
    res = []
    for fname in [n for n in os.listdir(resdir) if n.endswith(".res")]:
        params = [float(x) for x in os.path.splitext(fname)[0].split("_")[2:]]
        with open(resdir + sep + fname, "r") as f:
            i = 0
            for line in f:
                if i >= len(res):
                    res.append({})
                if params[0] not in res[i]:
                    res[i][params[0]] = {}
                res[i][params[0]][params[1]] = float(line.strip())
                i += 1
    return res

def plot_heatmap(mp, xticks=[], yticks=[], r=[], title="", xl="", yl=""):
    plt.figure()
    plt.imshow(mp)
    if xticks != []:
        plt.xticks(list(range(len(xticks))), xticks, rotation=30, ha='right')
    if yticks != []:
        plt.yticks(list(range(len(yticks))), yticks)
    if len(r) != 2:
        plt.imshow(mp, cmap="RdYlBu_r")
    else:
        plt.imshow(mp, vmin=r[0], vmax=r[1], cmap="RdYlBu_r")
    plt.title(title)
    plt.xlabel(xl)
    plt.ylabel(yl)
    plt.tight_layout()
    plt.show()
    return

if __name__ == "__main__":
    sep = "\\" if os.name == "nt" else "/"
    resdir = os.path.abspath(
        os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
    ) + sep + "output" + sep + "res"
    res = read_res(resdir)

    dists = []
    speeds = []
    for r in res:
        dists += list(r.keys())
        for v in r.values():
            speeds += list(v.keys())
    dists = list(set(dists))
    speeds = list(set(speeds))
    dists.sort()
    speeds.sort()
    speeds = speeds[::-1]

    mps = []
    for r in res:
        mp = []
        for s in speeds:
            tmp = []
            for d in dists:
                tmp.append(r[d][s])
            mp.append(tmp)
        mps.append(mp)
    mps = np.array(mps)

    nplanets = np.sum(mps > 0, axis=0)

    plot_heatmap(
        nplanets, ["{:.2e}".format(dist) for dist in dists], speeds,
        (0, len(mps)), "Amount of detached planets",
        "distance (km)", "speed (km/s)"
    )


