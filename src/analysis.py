"""
This files has analysis routines. They have to be in Python, because we want to tamper with them
"""

import numpy as np
import itertools as it
import ctypes as C

librdf = C.CDLL('./librdf.so')

def _f1(l):
    return np.arctan(np.sqrt(4*l**2 - 2))

def _f2(l):
    ang = 2*l*(4*l**2-3)/(np.sqrt(4*l**2-2)*(4*l**2+1))
    return 8*l*np.arctan(ang)
    
def _vol(l):
    if l < 0.5:
        v = 4.0/3.0 * np.pi * l**3
    elif l < 1.0/np.sqrt(2):
        v = -np.pi/12*(3-36*l**2+32*l**3)
    elif l < np.sqrt(3)/2:
        v = -np.pi/4 + 3 * np.pi * l**2
        v += np.sqrt(4 * l**2 - 2)
        v += (1 - 12 * l**2) * _f1(l)
        v += 2.0/3.0 * l**2 * _f2(l)
    return v

def rdf(x, size):
    """
    rdf that gets a box and calculates with PBC in 3d.

    Parameters:
    
    - x: 2D numpy array
         First dimension is the number of particles, second one is x y and z
    """
    nbins = 200
    tmp = (C.c_double * (nbins * 2))()
    natoms = np.shape(x)[0]
    x_p = x.ctypes.data_as(C.POINTER(C.c_double))
    librdf.rdf.argtypes = [C.POINTER(C.c_double), C.c_int, C.c_int, C.c_double, C.POINTER(C.c_double)]
    librdf.rdf(x_p, natoms, nbins, size, tmp)
    gr = np.frombuffer(tmp, dtype = np.double, count = nbins * 2)
    return gr.reshape((nbins, 2))


def rdf_py(x, box):
    """
    rdf that gets a box and calculates with PBC in 3d.

    Parameters:
    
    - x: 2D numpy array
         First dimension is the number of particles, second one is x y and z
    """
   
    if type(box) == float:
        xlo = -box/2
        ylo = -box/2
        zlo = -box/2
        xhi = box/2
        yhi = box/2
        zhi = box/2
    else:
        raise ValueError("Only support cubic boxes")
    
    vol = (xhi - xlo) * (yhi - ylo) * (zhi - zlo)
    rmax = np.sqrt(3)*box/2
    npoints = 200
    delr = rmax/(npoints+1)
    rdf = np.zeros((npoints, 2))
    rdf[:, 0] = np.linspace(0, rmax, npoints + 1)[:-1]
    natoms = np.shape(x)[0]
    for i in range(natoms):
        for j in range(i+1, natoms):
            dx = x[j] - x[i]
            if dx[0] > box/2: dx[0] -= box
            elif dx[0] < -box/2: dx[0] += box         
            if dx[1] > box/2: dx[1] -= box
            elif dx[1] < -box/2: dx[1] += box
            if dx[2] > box/2: dx[2] -= box
            elif dx[2] < -box/2: dx[2] += box
            r = np.linalg.norm(dx)
            idx = int(r*npoints/rmax)
            rdf[idx, 1] += 1
    density = natoms / vol
    for i in rdf:
        l = i[0]/box
        p = _vol(l + delr/box) - _vol(l)
        i[1] = i[1] / (p * density * box**3 * (natoms/2))
    return rdf


def pbc_rdf(x, box):
    """
    rdf that gets a box and calculates with PBC in 3d, replicating as needed

    Parameters:
    
    - x: 2D numpy array
         First dimension is the number of particles, second one is x y and z
    """
    if type(box) == float:
        xlo = -box/2
        ylo = -box/2
        zlo = -box/2
        xhi = box/2
        yhi = box/2
        zhi = box/2
    else:
        raise ValueError("Only support cubic boxes")
    
    natoms = np.shape(x)[0]
    base = (box * np.linspace(0, 2, 3))**2
    rmax = 0.5*box
    vol = rmax * rmax * rmax
    npoints = 1000
    delr = rmax/npoints
    rdf = np.zeros((npoints, 2))
    rdf[:, 0] = np.linspace(0, rmax, npoints)
    hits = 0
    for i in range(natoms):
        for j in range(natoms):
            for ix, iy, iz in it.product(range(-1, 2), range(-1, 2), range(-1, 2)):
                if i == j and ix == 0 and iy == 0 and iz == 0: continue
                d2 = base[abs(ix)] + base[abs(iy)] + base[abs(iz)]
                dx = x[j, 0] - x[i, 0]
                dy = x[j, 1] - x[i, 1]
                dz = x[j, 2] - x[i, 2]
                dd = ix * dx + iy * dy + iz * dz
                dd *= 2*box
                dr = x[j] - x[i]
                r = np.sqrt(np.dot(dr, dr) + d2 + dd)
                idx = int(r*npoints/rmax)
                if r < rmax:
                    hits += 1
                    rdf[idx, 1] += 1
    
    c = 4*np.pi / (3 * vol)
    for i in rdf:
        dV = (i[0] + delr) ** 3 - i[0]**3
        i[1] = i[1]  / (dV * c * hits)

    return rdf

def ssf(x):
    """
    get the form_factor for atoms in all three dimensions
    """
    npoints = 1000
    natoms, dim = np.shape(x)
    s = np.zeros((npoints/2+1, 2*dim), dtype=np.complex)
    for i in range(dim):
        h = np.histogram(x[:, i], npoints)
        dens = h[0]/float(natoms)
        r = h[1]
        s[:, 2*i] = np.fft.rfftfreq(npoints, (r[1] - r[0]))
        s[:, 2*i+1] = np.fft.rfft(dens)
    return s
    

if __name__ == '__main__':
    import extract as E
    x = E.positions('dump.lammpstrj', 0)
    rdf = rdf(x, 23.9397*2)
