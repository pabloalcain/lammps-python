"""
Logger class
"""

import fcntl
import os
import random
from collections import OrderedDict as od

def _create_path(system, style='folder'):
  """
  Create a path from the system.

  Parameters
  ----------

  system : neutronstar System
      The system to create path from

  style : {'folder', 'hash'}
      Style to use. Folder creates a path that mimics a folder
      structure similar to that of the system parameters. Hash creates
      a hash integer from the system dictionary.

  Returns
  -------

  path : string
      Path relative to root_path in which we will insert the data
  """
  prefix = od([('lambda', 'l'), ('N', 'N'), ('expansion', 'exp'),
               ('potential', ''), ('x', 'x'), ('d', 'd'), ('T', 'T')])

  if style == 'folder':
    path = ''
    for key in prefix.keys():
      if key in system.keys():
        path = path + ''.join((prefix[key], str(system[key]))) + '/'
  elif style == 'hash':
    path = str(hash(frozenset(system)))
  else:
    raise ValueError("The style {0} wasn't found".format(style))

  return path

class Logger(object):
  """
  Main logger class. Can either plot or write the data on a text file
  """
  def __init__(self, system, root_path='./data'):
    """
    Constructor.

    Parameters
    ----------

    system : System
        The molecular dynamics system. From it we will extract the
        global parameters that define the simulation

    root_path : str, optional
        Root path in which we will store all the information
    """
    identifier = [' '.join(map(str, (i, system[i]))) for i in system]
    path = _create_path(system, 'folder')
    self.path = '{0}/{1}'.format(root_path, path)
    os.makedirs(self.path)
    identifier.append('id {0}'.format(path))
    success = False
    while not success:
      try:
        fdb = open('{0}/db.dat'.format(root_path), 'a')
        fcntl.flock(fdb, fcntl.LOCK_EX)
        print>>fdb, ', '.join(identifier)
        fcntl.flock(fdb, fcntl.LOCK_UN)
        fdb.close()
        success = True
      except IOError:
        pass
    fkey = open('{0}/key.dat'.format(self.path), 'w')
    print>>fkey, ', '.join(identifier)
    fkey.close()
    system.lmp.command('log {0}/log.lammps'.format(self.path))

  def dump(self, system, style='text'):
    """
    Dump the system information.

    Parameters
    ----------

    system : System
        The molecular dynamics system. It needs to have a dump method
        in it. It can be eventually changed so it can be ported
        outside lammps.

    style : {'text', 'image'}
        Style to use when dumping
    """
    system.dump('{0}'.format(self.path), style)

  def log(self, analyzer):
    """
    Write the logging of the analyzer.

    Parameters
    ----------

    analyzer : Analyzer
        The analyzer to be logged
    """
    analyzer.log('{0}'.format(self.path))

  def plot(self, analyzer):
    """
    Plot the result of the analyzer.

    Parameters
    ----------

    analyzer : Analyzer
        The analyzer to be logged
    """
    analyzer.plot('{0}'.format(self.path))
