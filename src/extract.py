"""
This file has routines that extract information from dump lammps files
"""
from numpy import *

def positions(fname, step):
    f = open(fname)
    n = 0
    for line in f.readlines():
        n += 1
        if n < 4: continue
        if n == 4:
            npart = int(line)
            s1 = 9 * (step + 1) +  step * npart + 1
            x = zeros((npart, 3))
            continue
        if n < s1: continue
        np = n - s1
        for j in range(3): x[np, j] = line.split()[j + 2]
        if n == s1 + npart -1: break
    f.close()
    return x

if __name__ == "__main__":
    x = positions("dump.lammpstrj", 1)
