"""
Module for post-process in lammps
"""
import numpy as np

class Extraction(object):
  """
  Extraction class. We can extract either scalars, vector or
  per-particle magnitudes
  """
  def __init__(self, path='.'):
    self.path = path
    self.db = []
    try:
      f = open('{0}/db.dat'.format(path))
    except IOError:
      raise IOError("File db.dat not found in path {0}".format(path))
    for line in f.readlines():
      pa = line[:-1].split(', ')
      item = {}
      for it in pa:
        key, val = it.split(' ')
        try:
          val = float(val)
        except ValueError:
          pass
        item[key] = val
      self.db.append(item)

  def entries(self, parameters):
    """
    Return all the entries from the database that match the parameters.

    Parameters
    ----------

    parameters : dict
        A dictionary of parameters to look for in the database

    Returns
    -------

    ret : list
        A list of the entries (directories) in the database.
    """
    ret = []
    for entry in self.db:
      append = True
      for key, val in zip(parameters.keys(), parameters.values()):
        append = append and (val == entry[key])
      if append:
        ret.append(entry)
    return ret

  def particle(self, cols, parameters, dtype=np.float64, idx=0):
    """Extract a 'per particle' magnitude that is in a lammps dump file
    that matches the parameters.

    Parameters
    ----------

    cols : tuple
        Columns to extract

    parameters : dict
        A dictionary of parameters from which to extract

    dtype : datatype, optional
        Datatype of the array

    idx : datatype, optional
        Index of the timestep to look for [not the timestep itself]

    Returns
    -------

    output : list
        A list of numpy arrays with the requested data that match all
        parameters.
    """
    _size = len(cols)
    output = []
    for item in self.entries(parameters):
      path = '{0}/{1}'.format(self.path, item['id'])
      _file = open('{0}/dump.lammpstrj'.format(path))
      npart = _file.readlines
      _nline = 0
      for line in _file.readlines():
        _nline += 1
        if _nline < 4:
          continue
        if _nline == 4:
          npart = int(line)
          out = np.zeros((npart, _size), dtype=dtype)
          offset = (npart + 9) * idx
          continue
        if _nline < offset + 10:
          continue
        if _nline == offset + npart + 10:
          break
        idxp = _nline - offset - 10
        for i, j in enumerate(cols):
          out[idxp, i] = line.split()[j]
      _file.close()
      output.append(out)
    return output

  def x(self, parameters, idx=0):
    """Extract the positions that are in a lammps dump file
    that matches the parameters.

    Parameters
    ----------

    parameters : dict
        A dictionary of parameters from which to extract

    idx : datatype, optional
        Index of the timestep to look for [not the timestep itself]

    Returns
    -------

    x : list
        A list of numpy arrays with the positions that match all
        parameters.
    """
    return self.particle((2, 3, 4), parameters, np.float64, idx)

  def v(self, parameters, idx=0):
    """Extract the velocities that are in a lammps dump file
    that matches the parameters.

    Parameters
    ----------

    parameters : dict
        A dictionary of parameters from which to extract

    idx : datatype, optional
        Index of the timestep to look for [not the timestep itself]

    Returns
    -------

    v : list
        A list of numpy arrays with the velocities that match all
        parameters.
    """
    return self.particle((5, 6, 7), parameters, np.float64, idx)

  def t(self, parameters, idx=0):
    """Extract the types that are in a lammps dump file
    that matches the parameters.

    Parameters
    ----------

    parameters : dict
        A dictionary of parameters from which to extract

    idx : datatype, optional
        Index of the timestep to look for [not the timestep itself]

    Returns
    -------

    t : list
        A list of numpy arrays with the types that match all
        parameters.
    """
    return self.particle((1,), parameters, np.int32, idx)

  def box(self, parameters, idx=0):
    """Extract the boxes that are in a lammps dump file
    that matches the parameters.

    Parameters
    ----------

    parameters : dict
        A dictionary of parameters from which to extract

    idx : datatype, optional
        Index of the timestep to look for [not the timestep itself]

    Returns
    -------

    output : list
        A list of numpy arrays with the boxes that match all
        parameters.
    """
    output = []
    for item in self.entries(parameters):
      path = '{0}/{1}'.format(self.path, item['id'])
      _file = open('{0}/dump.lammpstrj'.format(path))
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
      output.append(size)
    return output
