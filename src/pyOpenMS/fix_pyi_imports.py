import sys


if __name__ == '__main__':
    args = sys.argv[1:]
    with open(args[0], 'r') as contents:
        save = contents.read()
    with open(args[0]+".bup", 'w') as contents:
        contents.write('from pyopenms import *')
        contents.write(save)

