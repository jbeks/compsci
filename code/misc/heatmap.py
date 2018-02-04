import os
import numpy as np
import matplotlib.pyplot as plt

# read res files into list of dictionaries
def read_res(resdir):
    # determine appropriate parameter
    sep = "\\" if os.name == "nt" else "/"
    res = []
    for fname in [n for n in os.listdir(resdir) if n.endswith(".res")]:
        # get parameters (dist, speed) from filename
        params = [float(x) for x in os.path.splitext(fname)[0].split("_")[2:]]
        with open(resdir + sep + fname, "r") as f:
            i = 0
            for line in f:
                # add dictionary to list for each planet
                if i >= len(res):
                    res.append({})
                # add dictionary for first parameter
                if params[0] not in res[i]:
                    res[i][params[0]] = {}
                # add data to dict
                res[i][params[0]][params[1]] = float(line.strip())
                i += 1
    return res

# create heatmap of data in mp
def plot_heatmap(
    mp, xticks=[], yticks=[], r=[], title="", xl="", yl="", save=""
):
    # set range if no range is given
    if len(r) != 2:
        r = [np.min(mp), np.max(mp)]
    # create heatmap using imshow
    plt.figure()
    plt.imshow(mp, vmin=r[0], vmax=r[1], cmap="RdYlBu_r")
    # create empty heatmap with created/given scale to use for colorbar
    heatmap = plt.pcolor(
        [[0]], vmin=r[0], vmax=r[1], cmap="RdYlBu_r", visible=False
    )
    plt.title(title)
    plt.xlabel(xl)
    plt.ylabel(yl)
    # add tick-labels if given
    if xticks != []:
        plt.xticks(list(range(len(xticks))), xticks)
    if yticks != []:
        plt.yticks(list(range(len(yticks))), yticks)
    cb = plt.colorbar(heatmap)
    cb.set_label("N planets")
    plt.tight_layout()
    # save figure if asked to
    if save != "":
        plt.savefig(save)
    plt.show()
    return

if __name__ == "__main__":
    # determine appropriate parameter
    sep = "\\" if os.name == "nt" else "/"

    # create dictionaries with res-data
    rootdir = os.path.abspath(
        os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
    )
    resdir = rootdir + sep + "output" + sep + "res"
    res = read_res(resdir)

    figdir = rootdir + sep + "output" + sep + "figs"
    figloc = figdir + sep + "heatmap.png"

    # convert dictionary into spead and distance lists
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

    # create heatmap data per planet
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

    # create heatmap for number of planets
    nplanets = np.sum(mps > 0, axis=0)

    # plot heatmap
    plot_heatmap(
        nplanets,
        [dist / 1e9 for dist in dists],
        [int(speed) for speed in speeds],
        (0, len(mps)),
        "Number of detached planets",
        r"Distance from barycenter ($10^9$ km)",
        "Speed (km/s)",
        figloc
    )

