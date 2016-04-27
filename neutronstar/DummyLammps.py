"""
Dummy lammps class, good for debugging purposes. It only logs
everything on a file
"""

class Lammps(object):
  def __init__(self, args, vals):
    """
    Constructor. Can be overriden with lammps real class

    Parameters
    ----------

    args : iterable
        The iterable of arguments

    vals : iterable
        The iterable of values
    """
    self.fout = open('log_dummy', 'w')
    print>>self.fout, args, vals

  def command(self, command):
    """
    Run a command. The dummy class only logs it.

    Parameters
    ----------

    command :
        Command to run
    """
    print>>self.fout, command
