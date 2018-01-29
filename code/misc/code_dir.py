import sys
from os import path, walk, name

code_dir = path.dirname(path.dirname(path.abspath(__file__)))
for directory in next(walk(code_dir))[1]:
    if name == 'nt':
        separator = "\\"
    else:
        separator = "/"
    sys.path.append(code_dir + separator + directory)

