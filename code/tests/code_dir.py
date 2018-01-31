import sys
from os import path, walk, name

# get path to "code" directory
code_dir = path.dirname(path.dirname(path.abspath(__file__)))
# determine system separator
if name == 'nt':
    separator = "\\"
else:
    separator = "/"
# add all directories in the "code" directory to the includepath
for directory in next(walk(code_dir))[1]:
    sys.path.append(code_dir + separator + directory)

