"""
MDSys: LAMMPS python wrapper for neutronstars
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import pylab as pl
from lammps import lammps
from random import randint
from os import makedirs, listdir
import neutronstar.analysis as A

_keys = ["potential", "lambda", "x", "N", "d", "T"]

def build_table(pot, l):
  """
  This method builds the actual potential table that will be read
  in lammps. So far, three different potentials are going to be
  supported:

  - Pandha medium
  - Pandha stiff
  - Horowitz
  """

  rc_nuc = 5.4
  rc_cou = max(l, rc_nuc)
  N = 5000
  pairs = ['NN', 'NP', 'PP']
  r = {}
  V = {}
  F = {}
  descr = {}
  Ncou = N
  r['NN'] = np.linspace(0, rc_nuc, N+1)[1:]
  r['NP'] = np.linspace(0, rc_nuc, N+1)[1:]
  r['PP'] = np.linspace(0, rc_cou, Ncou+1)[1:]

  Vc = 1.44
  uc = 1.0/l

  if pot in ["medium", "stiff", "newmed"]:
    if pot == "medium":
      (Vr, Va, V0) = (3088.118, 2666.647, 373.118)
      (ur, ua, u0) = (1.7468, 1.6, 1.5)

    if pot == "newmed": ## PARECE QUE ESTA MAL NEWMED
      (Vr, Va, V0) = (3097.0, 2696.0, 379.5)
      (ur, ua, u0) = (1.648, 1.528, 1.628)

    if pot == "stiff":
      (Vr, Va, V0) = (3601.482, 2834.338, 17630.256)
      (ur, ua, u0) = (2.2395, 2.0, 3.25)

    def Vnn(r):
      return V0 * np.exp(-u0 * r) / r
    
    def Vnp(r):
      return Vr * np.exp(-ur * r) / r - Va * np.exp(-ua * r) / r
      
    def Vpp(r):
      return V0 * np.exp(-u0 * r) / r + Vc * np.exp(-uc * r) / r


  elif pot == "horowitz":
    a = 110.0
    b = -26.0
    c = 24.0
    L = 1.25
    def Vnn(r):
      return a * np.exp(-r**2 / L) + (b + c) * np.exp(-r**2 / (2*L))

    def Vnp(r):
      return a * np.exp(-r**2 / L) + (b - c) * np.exp(-r**2 / (2*L))

    def Vpp(r):
      return a * np.exp(-r**2 / L) + (b + c) * np.exp(-r**2 / (2*L)) + Vc * np.exp(-uc * r) / r

  else:
    raise AttributeError("Option {0} for potential not found".format(pot))

  V['NN'] = Vnn(r['NN'])
  V['NP'] = Vnp(r['NP'])
  V['PP'] = Vpp(r['PP'])

  descr['NN'] = "# {0} potential for same species".format(pot)
  descr['NP'] = "# {0} potential for different species".format(pot)
  descr['PP'] = ("# {0} potential for same species "
                 "with Coulomb interaction lambda = {1}").format(pot, l)


  with open("potential.table", 'w') as fp:
    for p in pairs:
      print>>fp, descr[p]+"\n"
      print>>fp, p
      print>>fp, "N {0}\n".format(N)
      _r = r[p]
      _V = V[p]
      _F = - np.diff(_V) / np.diff(_r)
      _F.resize(N)
      for i in xrange(N): # improve with enumerate
        print>>fp, i+1, _r[i], _V[i], _F[i]



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
    for i in _keys:
      self.parameters[i] = None


  def setup(self, parameters, computes):
    """
    Setup: set the parameters and computes.

    * parameters: dictionary with values of lambda, N, potential, x, density
    and temperature
    * computes: list of computes
    * root: root directory for all information
    """
    self.set_path(parameters)
    self.set_parameters(parameters)
    self.set_computes(computes)
    # In the beginning, tally everything
    self.lmp.command("run 0 pre yes post no")

  def set_path(self, parameters):
    """
    Update this_path according to value of parameters
    """

    _prefix = {'potential': '',
               'lambda': 'l',
               'x': 'x',
               'N': 'N',
               'd': 'd',
               'T': 'T'}
    self.path = self.root
    for i in _keys:
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

    for i in _keys:
      if parameters[i] != self.parameters[i]:
        self.update(i, parameters[i])

  def update(self, key, value):
    """
    Different commands for the update of parameters
    """

    if key not in _keys:
      raise KeyError("Key {0} does not exist.")

    self.parameters[key] = value
    if key in ['potential', 'lambda']:
      _pot = self.parameters['potential']
      _l = self.parameters['lambda']
      if _l and _pot: 
        build_table(_pot, _l)
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
                   'create_atoms 1 random {n1} {s} box'.format(n1=nprot, s=randint(0, 10000)),
                   'create_atoms 2 random {n2} {s} box'.format(n2=nneut, s=randint(0, 10000)))
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

  def set_computes(self, computes):
    """
    Set computes and re-initialize them if needed
    """
    self.n_tally = 0
    if "rdf" in computes:
      self.c_rdf = 0
    if "mste" in computes:
      self.c_mste = 0
    if "ssf" in computes:
      self.c_ssf = 0
    self.variables = []

    self.computes = computes
    if not 'rdf' in computes:
      if 'ssf' in computes:
        raise AttributeError("Cannot calculate structure factor without rdf")
      if 'fit' in computes:
        raise AttributeError("Cannot fit without rdf")

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

  def read_dump(self, fname):

    """
    Read information from the last snapshot of a dump file
    """

    _t = None
    with open(fname, 'r') as fp:
      for line in fp:
        if line == "ITEM: TIMESTEP\n":
          _t = fp.next()[:-1]

    if not _t:
      raise IOError("File {0} does not look correct.".format(fname))

    cmd = ('read_dump {0} {1} x y z vx vy vz purge yes '
           'add yes replace no'.format(fname, _t))

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
      [temp, ke, epair, etot, press] = A.thermo(self.lmp, _N)
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
    if rmax == None: 
      rmax = (float(_N)/_d)**(1.0/3)*0.5

    self.n_tally += 1
    n = float(self.n_tally)

    if "rdf" in self.computes:
      t_r = A.rdf(self.lmp, nbins, rmax, self.npairs)
      self.c_rdf *= (n-1)/n
      self.c_rdf += t_r/n

    if "ssf" in self.computes:
      [t_s, a, b] = A.structure(t_r, _d, self.npairs)
      self.c_ssf *= (n-1)/n
      self.c_ssf += t_s/n
      self.collective['k_absorption'] = a
      self.collective['S_absorption'] = b

    if "fit" in self.computes:
      [a, b, c, d] = A.fit(t_r)
      self.collective['height'] = a
      self.collective['del_height'] = b
      self.collective['lambda'] = c
      self.collective['del_lambda'] = d

    if "mste" in self.computes:
      [t_c, a, b] = A.mste(self.lmp, _N)
      self.c_mste *= (n-1)/n
      self.c_mste += t_c/n
      self.collective['size_avg'] = a
      self.collective['size_std'] = b

    if "mink" in self.computes:
      [a, b, c, d] = A.minkowski(self.lmp, r_mink, r_cell)
      self.collective['volume'] = a
      self.collective['surface'] = b
      self.collective['breadth'] = c
      self.collective['euler'] = d


    if "thermo" in self.computes:
      [a, b, c, d, e] = A.thermo(self.lmp, _N)
      self.collective['temperature'] = a
      self.collective['kinetic'] = b
      self.collective['potential'] = c
      self.collective['energy'] = d
      self.collective['pressure'] = e

    if "lind" in self.computes:
      self.lind()

    self.log()

  def dump(self):
    """
    Wrapper to dump positions
    """
    path = self.path
    dump_fname = path + 'dump.lammpstrj'
    tmp = "dump myDUMP all custom 1 {0} id type x y z vx vy vz"
    self.lmp.command(tmp.format(dump_fname))
    self.lmp.command("dump_modify myDUMP sort id append yes")
    self.lmp.command("run 0 pre yes post no")
    self.lmp.command("undump myDUMP")

  def flush(self):
    """
    Write log in the data file and plot
    """
    path = self.path
    if "mste" in self.computes:
      # To add cluster size to file
      x = np.array(range(len(self.c_mste)))
      temp = np.vstack((x, self.c_mste))
      mste_fname = path + 'cluster.dat'
      np.savetxt(mste_fname, temp.T, header='size, number', fmt='%6i %1.4e')
      idx = self.c_mste.nonzero()[0]
      pl.figure()
      pl.loglog(x[idx], self.c_mste[idx], 'o-')
      pl.xlabel('Cluster size')
      pl.ylabel('Frequency')
      pl.tight_layout()
      pl.savefig(path + 'cluster.pdf')
      pl.close()


    if "rdf" in self.computes:
      rdf_fname = path + 'rdf.dat'
      h = 'r, a-a, ia-a, 1-1, i1-1, 1-2, i1-2, 2-2, i2-2'
      np.savetxt(rdf_fname, self.c_rdf, header=h, fmt='%1.4e')
      pairs = ['a-a', '1-1', '1-2', '2-2']
      for i in range(4):
        fig = pl.figure()
        pl.plot(self.c_rdf[:, 0], self.c_rdf[:, i*2+1], 'o-')
        pl.xlabel('Distance [fm]')
        pl.ylabel('RDF({0})'.format(pairs[i]))
        pl.tight_layout()
        pl.savefig(path + 'rdf_{0}'.format(pairs[i]))
        pl.close()

    if "ssf" in self.computes:
      ssf_fname = path + 'ssf.dat'
      h = 'r, a-a, 1-1, 1-2, 2-2'
      np.savetxt(ssf_fname, self.c_ssf, header=h)
      pairs = ['a-a', '1-1', '1-2', '2-2']
      for i in range(4):
        fig = pl.figure()
        pl.plot(self.c_ssf[:, 0], self.c_ssf[:, i+1], 'o-')
        pl.xlabel(r'Wave number [fm$^{-1}$]')
        pl.ylabel('SSF_{0}'.format(pairs[i]))
        pl.tight_layout()
        pl.savefig(path + 'ssf_{0}'.format(pairs[i]))
        pl.close()

    if "thermo" in self.computes:
      thermo_fname = path + 'thermo.dat'
      h = ', '.join(self.collective.keys())
      np.savetxt(thermo_fname, self.variables, header=h, fmt='%1.4e')

    self.set_computes(self.computes)
    self.lmp.command("reset_timestep 0")

  def log(self):
    """
    Append new data to variables. This has a problem with memory
    usage, and the fact that append may render slow for big
    arrays. Should be checked afterwards for long runs.

    Thoughts: log might upgrade the mean of each column + std dev on
    the fly.
    """
    self.variables.append(self.collective.values())
