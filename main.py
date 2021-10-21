# This file is the command line for the QCEMU.
import sys
from Parser import readFile

if __name__ == '__main__':
    filename = sys.argv[1]
    readFile(filename)
    