"""
MDSys: LAMMPS python wrapper for neutronstars
"""
import numpy as np
from lammps import lammps
from random import randint
from os import makedirs, listdir
import neutronstar.analysis as analysis
import neutronstar.graphics as graphics

_param_keys = ('potential', 'lambda', 'x', 'N', 'd', 'T')
_comp_keys = ('rdf', 'ssf', 'fit', 'mste', 'mink', 'thermo')


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
      for i, (ri, vi, fi) in enumerate(zip(_r, _V, _F)):
        print>>fp, i+1, ri, vi, fi


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
    * root: root directory for all information
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

    _prefix = {'potential': '',
               'lambda': 'l',
               'x': 'x',
               'N': 'N',
               'd': 'd',
               'T': 'T'}
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

    _is_rdf = 'rdf' in list_computes
    _is_ssf = 'ssf' in list_computes
    _is_fit = 'fit' in list_computes
    if not _is_rdf and (_is_ssf or _is_fit):
      raise AttributeError("Cannot calculate structure factor nor fit without rdf")

    self.n_tally = 0
    self.computes = {}
    
    _var = {'rdf': (),
            'ssf': ('k_absorption', 'S_absorption'),
            'fit': ('height', 'del_height', 
                    'lambda', 'del_lambda'),
            'mste': ('size_avg', 'size_std'),
            'mink': ('volume', 'surface', 'breadth', 'euler'),
            'thermo': ('temperature', 'kinetic',
                       'potential', 'energy', 
                       'pressure'),
            }

    _comp = {'rdf': 0,
             'ssf': 0,
             'fit': (),
             'mste': [0, 0],
             'mink': (),
             'thermo': (),
             }
    
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
    if rmax == None: 
      rmax = (float(_N)/_d)**(1.0/3)*0.5

    self.n_tally += 1
    n = float(self.n_tally)

    if "rdf" in self.computes:
      t_r = analysis.rdf(self.lmp, nbins, rmax, self.npairs)
      self.computes['rdf'] *= (n-1)/n
      self.computes['rdf'] += t_r/n

    if "ssf" in self.computes:
      [t_s, k, s] = analysis.structure(t_r, _d, self.npairs)
      self.computes['ssf'] *= (n-1)/n
      self.computes['ssf'] += t_s/n
      self.collective['k_absorption'].append(k)
      self.collective['S_absorption'].append(s)

    if "fit" in self.computes:
      [hei, dh, lam, dl] = analysis.fit(t_r)
      self.collective['height'].append(hei)
      self.collective['del_height'].append(dh)
      self.collective['lambda'].append(lam)
      self.collective['del_lambda'].append(dl)

    if "mste" in self.computes:
      [t_c, avg, std] = analysis.mste(self.lmp, _N)

      # Normalization of the proton fraction
      self.computes['mste'][1] *= (n-1) * self.computes['mste'][0]
      self.computes['mste'][1] += t_c[:,1] * t_c[:, 0]
      self.computes['mste'][1] /= (n-1) * self.computes['mste'][0] + t_c[:, 0]
      self.computes['mste'][1] = np.nan_to_num(self.computes['mste'][1])
      
      # Normalization of the occupancy
      self.computes['mste'][0] *= (n-1)/n
      self.computes['mste'][0] += t_c[:, 0]/n

      self.collective['size_avg'].append(avg)
      self.collective['size_std'].append(std)

    if "mink" in self.computes:
      [vol, sur, bre, eul] = analysis.minkowski(self.lmp, r_mink, r_cell)
      self.collective['volume'].append(vol)
      self.collective['surface'].append(sur)
      self.collective['breadth'].append(bre)
      self.collective['euler'].append(eul)

    if "thermo" in self.computes:
      [tem, kin, pot, ene, pre] = analysis.thermo(self.lmp, _N)
      self.collective['temperature'].append(tem)
      self.collective['kinetic'].append(kin)
      self.collective['potential'].append(pot)
      self.collective['energy'].append(ene)
      self.collective['pressure'].append(pre)

    if "lind" in self.computes:
      self.lind()

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
    if "mste" in self.computes: 
      graphics.mste(self.computes['mste'][0], 
                    self.computes['mste'][1], path=path)

    if "rdf" in self.computes: graphics.rdf(self.computes['rdf'], path=path)
    if "ssf" in self.computes: graphics.ssf(self.computes['ssf'], path=path)
    if "thermo" in self.computes: graphics.thermo(self.collective, path=path)

    self.set_computes(self.computes)
    self.lmp.command("reset_timestep 0")
    
