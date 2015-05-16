"""
MDSys: LAMMPS python wrapper for neutronstars
"""
import numpy as np
from lammps import lammps
from random import randint
from os import makedirs, listdir
import neutronstar.analysis as analysis
import neutronstar.graphics as graphics
import neutronstar.results as results
import neutronstar.potential as potential

_param_keys = ('potential', 'lambda', 'x', 'N', 'd', 'T')
_comp_keys = ('rdf', 'ssf', 'fit', 'mste', 'mink', 'thermo')
_var = {'rdf': (),
        'ssf': ('k_absorption', 'S_absorption'),
        'fit': ('height', 'del_height', 
                'lambda', 'del_lambda'),
        'mste': ('size_avg', 'size_std'),
        'mink': ('volume', 'surface', 'breadth', 'euler'),
        'thermo': ('temperature', 'kinetic',
                   'potential', 'energy', 
                   'pressure'),}

_comp = {'rdf': 0,
         'ssf': 0,
         'fit': (),
         'mste': [0, 0],
         'mink': (),
         'thermo': (),}
_prefix = {'potential': '',
           'lambda': 'l',
           'x': 'x',
           'N': 'N',
           'd': 'd',
           'T': 'T',}


class MDSys(object):
  def __init__(self, gpu=False, silent=True, root='./data', log='log.lammps'):
    """
    Constructor: Instantiate the lammps class, so the system is always
    aware of the object it has.

    Pass dummy info to lammps and set root directory

    gpu = Use LAMMPS gpu package
    silent = Run in silent mode (no output on screen)
    root = Root directory for file hierarchy
    """

    _args = ['-log', log]

    if silent:
      _args += ['-screen', 'none', '-nocite']

    if gpu:
      _args += ['-pk', 'gpu 1', '-sf', 'gpu']

    self.lmp = lammps("", _args)

    script = (
      "#Nuclear model",
      "units             lj",
      "atom_style        atomic",
      "timestep          0.10",
      "region            box block 0 1.0 0 1.0 0 1.0",
      "create_box        2 box",
      "mass              1 938.0",
      "mass              2 938.0",
      "pair_style        table linear {ninter}".format(ninter=5000),
      "neighbor          1.2 bin",
      "neigh_modify      every 1 delay 0 check yes one 8000 page 80000",
      "comm_modify       vel yes",
      "thermo_style      custom step temp ke epair etotal press",
      "thermo            1000",
      "compute           mste all mste/atom 5.4",
    )

    for cmd in script:
      self.lmp.command(cmd)

    self.root = root
    self.collective = {}
    self.npairs = 4
    self.parameters = {}
    for i in _param_keys:
      self.parameters[i] = None


  def setup(self, parameters, list_computes):
    """
    Setup: set the parameters and computes.

    * parameters: dictionary with values of lambda, N, potential, x, density
    and temperature
    * list_computes: list of computes
    """
    self.set_path(parameters)
    self.set_parameters(parameters)
    self.set_computes(list_computes)
    # In the beginning, tally everything
    self.lmp.command("run 0 pre yes post no")

  def set_path(self, parameters):
    """
    Update this_path according to value of parameters
    """

    self.path = self.root
    for i in _param_keys:
      self.path += "/{0}{1}".format(_prefix[i], parameters[i])
    self.path += "/"
    
    try:
      makedirs(self.path)
    except OSError:
      if len(listdir(self.path)) != 0:
        msg = ("Directory {0} already exists:"
               "rename base path or delete old files")
        raise OSError(msg.format(self.path))

  def set_parameters(self, parameters):
    """
    We only updates those that change
    """

    for i in _param_keys:
      if parameters[i] != self.parameters[i]:
        self.update(i, parameters[i])

  def update(self, key, value):
    """
    Different commands for the update of parameters
    """

    if key not in _param_keys:
      raise KeyError("Key {0} does not exist.")

    self.parameters[key] = value
    if key in ['potential', 'lambda']:
      _pot = self.parameters['potential']
      _l = self.parameters['lambda']
      if _l and _pot: 
        potential.build_table(_pot, _l)
        base = 'pair_coeff {t1} {t2} potential.table {pair} {cutoff}'
        command = (base.format(t1=1, t2=1, pair='NN', cutoff=5.4),
                   base.format(t1=1, t2=2, pair='NP', cutoff=5.4),
                   base.format(t1=2, t2=2, pair='PP', cutoff=max(5.4, float(_l))))
      else:
        command = ()
    elif key in ['x', 'N']:
      _x = self.parameters['x']
      _N = self.parameters['N']
      if _x and _N:
        nprot = int(_x * _N)
        nneut = int(_N) - nprot
        command = ('delete_atoms group all',
                   'create_atoms 1 random {n1} {s} box'.format(n1=nneut, s=randint(0, 10000)),
                   'create_atoms 2 random {n2} {s} box'.format(n2=nprot, s=randint(0, 10000)))
      else:
        command = ()
        
    elif key == "T":
      command = ('fix 1 all nvt temp {T} {T} 100.0'.format(T=value),)

    elif key == "d":
      _N = self.parameters['N']
      _vol = float(_N)/value
      _size = _vol ** (1.0 / 3.0)
      _size = _size / 2
      command = (('change_box all x final -{s} {s} '
                  'y final -{s} {s} '
                  'z final -{s} {s} remap').format(s=_size),)

    for cmd in command:
      self.lmp.command(cmd)

  def set_computes(self, list_computes):
    """
    Set computes and re-initialize them if needed
    """

    _r = 'rdf' in list_computes
    _s = 'ssf' in list_computes
    _f = 'fit' in list_computes
    if not _r and (_s or _f):
      raise AttributeError("Cannot calculate structure factor nor fit without rdf")

    self.n_tally = 0
    self.computes = {}
    
    for c in list_computes:
      self.computes[c] = _comp[c]
      for v in _var[c]:
        self.collective[v] = []

  def minimize(self):
    """
    Minimize to remove some of the potential energy, probably due to
    initial random configuration, and set temperature from Gaussian
    """
    _l = self.parameters['lambda']
    _coff = max(_l, 5.4)
    _T = self.parameters['T']
    _s = randint(0, 10000)
    command = ('pair_coeff 2 2 potential.table PP 5.4',
               'min_style hftn',
               'minimize 0 1.0 1000 100000',
               'pair_coeff 2 2 potential.table PP {c}'.format(c=_coff),
               'velocity all create {T} {s}'.format(T=_T, s=_s),
               'reset_timestep 0')

    for cmd in command:
      self.lmp.command(cmd)

  def expand(self, vel):
    """
    Set initial conditions for an expansion as seen in 

    Dorso and Strachan, Phys. Rev. B 54, 236
    """
    _N = self.parameters['N']
    _d = self.parameters['d']
    _vol = float(_N)/_d
    _size = _vol ** (1.0 / 3.0)
    _size = _size / 2
    command = (('fix expansion all deform 1 x vel {0} '
                'y vel {0} z vel {0} remap none units box').format(vel),
               'velocity all ramp vx 0 {0} x -{1} {1} sum yes units box'.format(vel, _size),
               'velocity all ramp vy 0 {0} y -{1} {1} sum yes units box'.format(vel, _size),
               'velocity all ramp vz 0 {0} z -{1} {1} sum yes units box'.format(vel, _size),
               'velocity all zero linear')
               
    for cmd in command:
      self.lmp.command(cmd)

  def unexpand(self):
    """
    Stop expansion
    """
    
    self.lmp.command("unfix expansion")

  def read_dump(self, fname, mste=False):

    """Read information from the last snapshot of a dump file. In case
    the dump file has the mste column, set mste to True

    """

    _t = None
    with open(fname, 'r') as fp:
      for line in fp:
        if line == "ITEM: TIMESTEP\n":
          _t = fp.next()[:-1]

    if not _t:
      raise IOError("File {0} does not look correct.".format(fname))

    if mste: mste_string = "c_mste"
    else: mste_string = ""
    
    cmd = ('read_dump {0} {1} x y z vx vy vz {2} purge yes '
           'add yes replace no'.format(fname, _t, mste_string))

    self.lmp.command(cmd)

  def run(self, Nsteps):
    """
    Wrapper for the "run" command in lammps
    """
    self.lmp.command("run {ns} pre no post yes".format(ns=Nsteps))


  def equilibrate(self, nfreq=300, wind=20):
    """
    This method takes care of the thermalization, with a Langevin
    thermostat.

    nfreq is how many timesteps to take between runs, and wind is the
    size of the window. The criterion for stability is that the
    average temperature of the last 100 steps is close to the set
    temperature by a standard deviation, while the energy stops
    decreasing (slope > 0). This works only when setting temperature
    from high to low (while for going from low to high, one would
    expect that just setting the temperatures should be enough).

    Thermalization doesn't write any log or dump file.
    """
    energy = np.zeros(wind)
    temperature = np.zeros(wind)
    step = np.zeros(wind)
    i = 0
    _N = self.parameters["N"]
    while True:
      i = i + 1
      self.run(nfreq)
      # Extract thermo values
      [temp, ke, epair, etot, press] = analysis.thermo(self.lmp, _N)
      # Add to the window
      energy[i % wind] = etot
      temperature[i % wind] = temp
      step[i % wind] = i
      # Only update values if the window is complete
      if i < wind:
        continue
      # Slow, we are calculating things again. Could be updated
      # and deleted each time. Check profiling!
      [slope, aux] = np.polyfit(step, energy, 1)
      diff = abs(self.parameters['T'] - np.mean(temperature))
      std = np.std(temperature)
      if slope > 0 and diff < std:
        break
    self.lmp.command("reset_timestep 0")
    # To tally everything since we reset timestep
    self.lmp.command("run 0 pre yes post no")

  def results(self,
	      r_mink=1.8, r_cell=0.5,
	      nbins=200, rmax=None):
    """
    Method to take all the results that have been set in the setup()
    method
    """
    _N = self.parameters['N']
    _d = self.parameters['d']
    if not rmax: rmax = (float(_N)/_d)**(1.0/3)*0.5

    self.n_tally += 1
    n = float(self.n_tally)

    # Go through the analysis
    if "rdf" in self.computes: rdf = analysis.rdf(self.lmp, nbins, rmax, self.npairs)
    if "ssf" in self.computes: ssf = analysis.structure(rdf, _d, self.npairs)
    if "fit" in self.computes: fit = analysis.fit(rdf)
    if "mste" in self.computes: mste = analysis.mste(self.lmp, _N)
    if "mink" in self.computes: mink = analysis.minkowski(self.lmp, r_mink, r_cell)
    if "thermo" in self.computes: thermo = analysis.thermo(self.lmp, _N)
    if "lind" in self.computes: lind = analysis.lind(self.lmp)

    # Tally everything
    if "rdf" in self.computes: results.rdf(self.computes, self.collective, rdf, n)
    if "ssf" in self.computes: results.ssf(self.computes, self.collective, ssf, n)
    if "fit" in self.computes: results.fit(self.computes, self.collective, fit, n)
    if "mste" in self.computes: results.mste(self.computes, self.collective, mste, n)
    if "mink" in self.computes: results.mink(self.computes, self.collective, mink, n)
    if "thermo" in self.computes: results.thermo(self.computes, self.collective, thermo, n)
    if "lind" in self.computes: results.lind(self.computes, self.collective, lind, n)

  def dump(self):
    """
    Wrapper to dump positions
    """
    path = self.path
    dump_fname = path + 'dump.lammpstrj'
    tmp = "dump myDUMP all custom 1 {0} id type x y z vx vy vz c_mste"
    self.lmp.command(tmp.format(dump_fname))
    self.lmp.command("dump_modify myDUMP sort id append yes")
    self.lmp.command("run 0 pre yes post no")
    self.lmp.command("undump myDUMP")

  def flush(self):
    """
    Write log in the data file and plot
    """
    path = self.path
    if "mste" in self.computes: graphics.mste(self.computes['mste'], path=path)
    if "rdf" in self.computes: graphics.rdf(self.computes['rdf'], path=path)
    if "ssf" in self.computes: graphics.ssf(self.computes['ssf'], path=path)
    if "thermo" in self.computes: graphics.thermo(self.collective, path=path)
    self.set_computes(self.computes)
    self.lmp.command("reset_timestep 0")
