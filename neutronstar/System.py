"""
System class: the Lammps Molecular Dynamics main object
"""

class System(object):
  def __init__(self, gpu=False, silent=True):
    """
    Constructor: Instantiate LAMMPS Class, that belongs to this
    object. Therefore, LAMMPS system will always belong to it.

    Parameters
    ----------

    gpu : boolean
        Use gpu acceleration package

    silent : boolean
        Don't print information on the screen
    """
    pass

  @property
  def parameters(self):
    """
    Parameters of the simulation
    """
    pass

  @parameters.setter
  def parameters(self, value):
    """
    Set parameters of the simulation
    """
    pass

  def minimize(self, etol, ftol, maxiter, maxeval):
    """
    Minimize to remove some of the potential energy, probably due to
    initial random configuration, and set temperature to Gaussian
    according to the system temperature.

    Parameters
    ----------

    etol : float
        Stopping tolerance for energy (unitless)
    ftol : float
        Stopping tolerance for force (force units)
    maxiter : int
        Maximum number of iterations of minimizer
    maxeval : int
        Maximum number of force evaulations
    """
    pass

  def expand(self, rate):
    """
    Set initial condition for an expansion as seen in [Dorso]_

    Parameters
    ----------

    rate : float
        Expansion rate, in box per seconds units

    .. [Dorso] Dorso and Strachan, Phys. Rev. B 54, 236
    """
    pass

  def unexpand(self):
    """
    Stop expansion
    """
    pass

  def read_dump(self, fname):
    """
    Read information from the last snapshot of a dump file. It purges
    all previous particles, so the particles are exactly the type and
    ids of the new particles are those of the dump file.

    Parameters
    ----------

    fname : str
        LAMMPS-style dump filename to read

    Raises
    ------

    IOError
        In case there is no snapshot in the filename to read
    """
    pass

  def thermalize(self, freq, wind):
    """Runs the system until it has thermalized. The criterion for
    stability is:

    1.- The average temperature of the last freq*steps is close to
    the set temperature by a standard deviation and

    2.- The energy stops decreasing (slope > 0)

    Parameters
    ----------

    freq: int
        Number of timesteps between each temperature to be considered
        in the average and slope

    wind: int
       Number of timesteps in the slope and average calculation of the
       temperature.

    Notes
    -----

    Due to criterion 2, this works only when setting temperature from
    high to low. Nevertheless, when going from low to high
    temperatures, one would expect that just setting the temperatures
    and complying with 1 should be enough.
    """
    pass
