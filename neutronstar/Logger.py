"""
Logger class
"""

import fcntl
import random

class Logger(object):
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
    key = hash(frozenset(system)) + random.randint(0, 1e10)
    identifier.append('id {0}'.format(key))
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
    self.path = '{0}/{1}'.format(root_path, key)

    fkey = open('{0}/key.dat'.format(self.path))
    print>>fkey, ', '.join(identifier)
    fkey.close()

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
