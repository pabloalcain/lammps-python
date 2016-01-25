"""
Module for post-process in neutronstars
"""
import numpy as np
from collections import OrderedDict
import os
import analysis as A

PREFIX = OrderedDict()
PREFIX['V'] = ''
PREFIX['l'] = 'l'
PREFIX['x'] = 'x'
PREFIX['N'] = 'N'
PREFIX['d'] = 'd'
PREFIX['T'] = 'T'

COLS = {'x': (2, 3, 4), 'v': (5, 6, 7), 'clus': (8,), 'type': (1,)}
DATA_TYPE = {'x': np.double, 'v': np.double, 'clus': np.int32, 'type': np.int32}

def _check(path):
  """Check whether path has only one possible value. If so, we don't
  need to set it from the outside and returns the only choice

  """

  if len(os.listdir(path)) > 1:
    helper = "Multiple options for directory {0}: {1}. Please specify a value"
    raise AttributeError(helper.format(path, os.listdir(path)))
  return os.listdir(path)[0]

class Extraction(object):
  """
  Extraction class. We can extract either scalars, vector or
  per-particle magnitudes

  """
  # TODO: make getter attributes and abstract the array/particle part
  def __init__(self, path='.'):
    self.path = path

  def _set_path(self, parameters):
    """Parameters is a dictionary with entries for the 6 parameters we
    use when we search in a root data directory; V, l, x, N, d, T

    Here we *believe* on the dir structure:
    {V}/l{l}/x{x}/N{N}/d{d}/T{T}
    """

    _path = self.path
    # for every possible parameter, we check whether it was passed or not
    for _param in PREFIX:
      if _param in parameters:
        direct = PREFIX[_param] + str(parameters[_param])
      else:
        direct = _check(_path)
      _path += '/' + direct
    return _path

  def scalar(self, mag, parameters, full=False):
    """
    Extract a single collective value with its stdev.
    It doesn't support anymore height, that came from the
    fit method. Now lambda makes the conversion from k to
    lambda. Possible magnitudes:
    - breadth
    - S_absorption
    - pressure
    - size_avg
    - k_absorption
    - surface
    - volume
    - size_std
    - potential
    - kinetic
    - energy
    - euler
    - temperature

    - lambda (conversion from k_absorption)
    
    if full, return all the values instead of avg +/- std
    """
    _path = self._set_path(parameters)
    try:
      _file = open(_path + "/thermo.dat")
    except IOError:
      return 0, 0
    header = _file.readline()[2:-1]
    _file.close()

    if mag == 'lambda':
      idx = header.split(', ').index('k_absorption')
    else:
      idx = header.split(', ').index(mag)
    therm = np.loadtxt(_path +'/thermo.dat')
    if full:
      if mag == 'lambda':
        return 2*np.pi/therm[:, idx]
      return therm[:, idx]
    avg = np.mean(therm, axis=0)[idx]
    std = np.std(therm, axis=0)[idx]
    if mag == 'lambda':
      std = 2*np.pi/(avg*avg) * std
      avg = 2*np.pi/avg
    return avg, std

  def array(self, mag, parameters):
    """
    Extract an array magnitude, that is in a textfile inside the
    directory. Possible magnitudes are:

    - cluster
    - rdf
    - ssf

    Obviously each one has a different structure, so we only pass the
    parsed array as a matrix
    """
    _path = self._set_path(parameters)
    try:
      out = np.loadtxt('{0}/{1}.dat'.format(_path, mag))
    except ValueError:
      out = np.loadtxt('{0}/{1}.dat'.format(_path, mag), delimiter=',')
    return out

  def mste(self, parameters, idx=0):
    """
    Compute MSTE from the dump file
    """
    x = self.particle('x', parameters, idx)
    v = self.particle('v', parameters, idx)
    t = self.particle('type', parameters, idx)
    box = self.box(parameters, idx)
    size = box[0, 1] - box[0, 0]
    index = A.cluster(x, v, t, size, energy=False)
    return index

  def rdf(self, parameters, idx=0):
    """
    Compute RDF from the dump file
    """
    x = self.particle('x', parameters, idx=idx)
    t = self.particle('type', parameters, idx=idx)
    box = self.box(parameters, idx=idx)
    size = box[0, 1] - box[0, 0]
    return A.rdf(x, t, size)

  def particle(self, mag, parameters, idx=0):
    """
    Extract a 'per particle' magnitude that is in a lammps dump file
    inside the directory. Get only the first timestep. Possible
    magnitudes are:

    - x
    - v
    - type
    - clus
    """
    _path = self._set_path(parameters)
    _size = len(COLS[mag])
    _file = open('{0}/dump.lammpstrj'.format(_path))
    npart = _file.readlines
    _nline = 0
    for line in _file.readlines():
      _nline += 1
      if _nline < 4:
        continue
      if _nline == 4:
        npart = int(line)
        out = np.zeros((npart, _size), dtype=DATA_TYPE[mag])
        offset = (npart + 9) * idx
        continue
      if _nline < offset + 10:
        continue
      if _nline == offset + npart + 10:
        break
      idxp = _nline - offset - 10
      for i, j in enumerate(COLS[mag]):
        out[idxp, i] = line.split()[j]
    _file.close()
    return out

  def box(self, parameters, idx=0):
    """
    Extract box at the first timestep from a lammps file
    """
    _path = self._set_path(parameters)
    _file = open('{0}/dump.lammpstrj'.format(_path))
    size = np.zeros((3, 2))
    _nline = 0
    for line in _file.readlines():
      _nline += 1
      if _nline < 4:
        continue
      if _nline == 4:
        npart = int(line)
        offset = (npart + 9) * idx
        continue
      if _nline > offset + 5 and _nline < offset + 9:
        size[_nline - offset - 6, :] = line.split()
        continue
      if _nline == offset + 9:
        break
    return size

def main():
  """Main function with a use case, we extract the energy and compare
  between horowitz and medium eos.
  """
  ext = Extraction('/home/pablo/artemis/eos/data')
  mag = 'energy'
  parameters = {'x': 0.5, 'd': 0.08, 'T': 0.5, 'V': 'horowitz'}
  horowitz = ext.scalar(mag, parameters)
  parameters['V'] = 'medium'
  medium = ext.scalar(mag, parameters)
  print 'Energy for medium: {0} +/- {1}'.format(*medium)
  print 'Energy for horowitz: {0} +/- {1}'.format(*horowitz)
  typ = ext.particle('type', parameters)
  print typ


if __name__ == "__main__":
  main()
