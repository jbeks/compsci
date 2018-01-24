import argparse

def read_ref(name):
    sep = 0
    data = []
    ax = ["X", " Y", " Z"]
    newfile = open("ref/" + name + ".ref", 'w+')
    with open("ref/" + name + "_data.ref", 'r') as f:
        for line in f:
            if line.startswith("Step-size"):
                split = line.split()
                stepsize = int(split[2]) * 60
                newfile.write("stepsize " + str(stepsize) + "\n")
            if "A.D." not in line:
                continue
            s = line.split(",")
            if len(s) < 8:
                continue
            newfile.write(s[2] + " " + s[3] + " " + s[4])
            newfile.write("\n")
    newfile.close()

def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "file", type=str,
        help="data file to parse"
    )
    args = parser.parse_args()
    return args.file

if __name__ == "__main__":
    filename = parse_arguments()
    read_ref(filename)
