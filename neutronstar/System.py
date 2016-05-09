"""
System class: the Lammps Molecular Dynamics main object
"""

from lammps import lammps#neutronstar.DummyLammps import Lammps as lammps
import neutronstar
from neutronstar import Potential

from random import randint
import numpy as np

class System(dict):
  """
  Any system will be a dictionary with the parameters of the
  simulation. It inherits from a dictionary to keep the syntax as
  clean as possible
  """
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

    _args = []
    if silent:
      _args += ['screen', 'none', '-nocite']

    if gpu:
      _args += ['-pk', 'gpu 1 neigh no', '-sf', 'gpu']

    self.lmp = lammps("", _args)
    super(System, self).__init__()

  def __setitem__(self, key, value):
    """
    The overriding of the setting of the item. We set it as a
    dictionary item and make all the changes needed in the
    simulation. We can set here any information of the system and, if
    it is identified as something that must change the simulation
    characteristic.
    """
    super(System, self).__setitem__(key, value)

  def minimize(self, etol, ftol, maxiter, maxeval):
    """
    Minimize to remove some of the potential energy, probably due to
    initial random configuration, and set temperature to Gaussian
    according to the system temperature if it exists.

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
    self.lmp.command("min_style hftn")
    _mi = "minimize {etol} {ftol} {maxiter} {maxeval}"
    self.lmp.command(_mi.format(etol=etol, ftol=ftol,
                                maxiter=maxiter, maxeval=maxeval))
    try:
      _temp = self['T']
      _seed = randint(0, 10000)
      _vel = "velocity all create {T} {s}"
      self.lmp.command(_vel.format(T=_temp, s=_seed))
    except KeyError:
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
    _npart = self['N']
    _density = self['d']
    _vol = float(_npart)/_density
    _size = _vol ** (1.0 / 3.0) / 2
    _vel = rate * _size
    _vc = ('velocity all ramp {vi} -{v} {v} {ri} -{s} {s} '
           'sum yes units box')
    self.lmp.command('fix expansion all deform 1')
    self.lmp.command(_vc.format(vi='vx', ri='x', v=_vel, s=_size))
    self.lmp.command(_vc.format(vi='vy', ri='y', v=_vel, s=_size))
    self.lmp.command(_vc.format(vi='vz', ri='z', v=_vel, s=_size))

  def unexpand(self):
    """
    Stop expansion
    """
    self.lmp.command("unfix expansion")

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

    _ts = None
    with open(fname, 'r') as filep:
      for line in filep:
        if line == "ITEM: TIMESTEP\n":
          _ts = filep.next()[:-1]

    if not _ts:
      raise IOError("File {0} does not look correct.".format(fname))

    self.lmp.command('read_dump {0} {1} x y z vx vy vz purge yes '
                     'add yes replace no'.format(fname, _ts))

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
    thermo = neutronstar.computes.Thermo.Thermo.Thermo()
    energy = np.zeros(wind)
    temperature = np.zeros(wind)
    step = np.zeros(wind)
    i = 0
    while True:
      i = i + 1
      self.run(freq)
      # Extract thermo values
      [temp, _, _, etot, _] = analysis.thermo(self)
      # Add to the window
      energy[i % wind] = etot
      temperature[i % wind] = temp
      step[i % wind] = i
      # Only update values if the window is complete
      if i < wind:
        continue
      [slope, _] = np.polyfit(step, energy, 1)
      diff = abs(self['T'] - np.mean(temperature))
      std = np.std(temperature)
      if slope > 0 and diff < std:
        break
    self.lmp.command("reset_timestep 0")
    # To tally everything since we reset timestep
    self.lmp.command("run 0 pre yes post no")

  @property
  def x(self):
    """
    Get x as a numpy array from the simulation
    """
    _data = np.array(self.lmp.gather_atoms("x", 1, 3))
    natoms = np.shape(_data)[0]
    return _data.reshape(natoms, 3)

  @property
  def v(self):
    """
    Get v as a numpy array from the simulation
    """
    _data = np.array(self.lmp.gather_atoms("v", 1, 3))
    natoms = np.shape(_data)[0]
    return _data.reshape(natoms, 3)

  @property
  def t(self):
    """
    Get types as a numpy array from the simulation
    """
    _data = np.array(self.lmp.gather_atoms("type", 0, 1), dtype=np.int)
    return _data

  def run(self, steps):
    """
    Wrapper for the "run" command in lammps.

    Parameters
    ----------

    steps: number of steps to run
    """
    self.lmp.command("run {ns} pre no post yes".format(ns=steps))

class NeutronStarSystem(System):
  """
  The Neutron Star System. It sets up the typical usage for Neutron
  Star Simulations
  """
  def __init__(self, gpu=False, silent=True):
    super(NeutronStarSystem, self).__init__(gpu, silent)
    _sc = ("#Nuclear model",
           "units             lj",
           "atom_style        atomic",
           "timestep          0.10",
           "region            box block 0 1.0 0 1.0 0 1.0",
           "create_box        2 box",
           "mass              1 938.0",
           "mass              2 938.0",
           "pair_style        table linear {ninter}".format(ninter=5000),
           "neighbor          5.2 multi",
           "neigh_modify      every 1 delay 0 check yes one 40000 page 400000",
           "comm_modify       vel yes",
           "thermo_style      custom step temp ke epair etotal press",
           "thermo            1000",
           "compute           mste all mste/atom 5.4",)
    for cmd in _sc:
      self.lmp.command(cmd)

  def __setitem__(self, key, value):
    """
    We define how to work when we override certain parameters of the
    system.
    """
    super(NeutronStarSystem, self).__setitem__(key, value)
    if key in ['potential', 'lambda']:
      # This is because we don't know if potential and lambda have
      # been set
      try:
        _pot = self['potential']
        _lambda = self['lambda']
        _name = Potential.build_table(_pot, _lambda)
        base = 'pair_coeff {t1} {t2} {fname} {pair} {cutoff}'
        _pp_cutoff = max(5.4, float(_lambda))
        _nn = base.format(t1=1, t2=1, fname=_name, pair='NN', cutoff=5.4)
        _np = base.format(t1=1, t2=2, fname=_name, pair='NP', cutoff=5.4)
        _pp = base.format(t1=2, t2=2, fname=_name, pair='PP', cutoff=_pp_cutoff)
        self.lmp.command(_nn)
        self.lmp.command(_np)
        self.lmp.command(_pp)
      except KeyError:
        pass

    elif key in ['x', 'N']:
      # This is because we don't know if potential and lambda have
      # been set
      try:
        _fraction = self['x']
        _npart = self['N']
        _np = int(_fraction * _npart)
        _nn = int(_npart) - _np
        _cr = 'create_atoms {t} random {n1} {s} box'
        self.lmp.command('delete_atoms group all')
        self.lmp.command(_cr.format(t=1, n1=_nn, s=randint(0, 10000)))
        self.lmp.command(_cr.format(t=2, n2=_np, s=randint(0, 10000)))
      except KeyError:
        pass

    elif key == "T":
      _temp = 'fix 1 all nvt temp {T} {T} 1000.0'
      self.lmp.command(_temp.format(T=value))

    elif key == "d":
      try:
        _npart = self['N']
        _vol = float(_npart)/value
        _size = _vol ** (1.0 / 3.0) / 2
        _cb = ('change_box all x final -{s} {s} y final -{s} {s} '
               'z final  {s} {s} remap')
        self.lmp.command(_cb.format(s=_size))
      except KeyError:
        pass
