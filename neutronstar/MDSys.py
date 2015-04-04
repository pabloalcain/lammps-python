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
from scipy.stats.distributions import t
from scipy.optimize import curve_fit
import neutronstar.analysis as A

_keys = ["potential", "lambda", "x", "N", "d", "T"]

def build_table(pot, l):
  """
  This method builds the actual potential table that will be read
  in lammps. So far, three different potentials are going to be
  supported:

  - Pandha medium
  - Pandha stiff
  - Horowitz (not yet implemented)

  The code here wasn't given much thought so maybe it can be improved
  by whoever wants to give it a try, and, at least, make it cleaner!
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
      (ur, ua, u0) = (1.5, 1.6, 1.7468)

    if pot == "newmed": ## PARECE QUE ESTA MAL NEWMED
      (Vr, Va, V0) = (2696.0, 3097.0, 379.5)
      (ur, ua, u0) = (1.528, 1.648, 1.628)

    if pot == "stiff":
      (Vr, Va, V0) = (2834.338, 3601.482, 17630.256)
      (ur, ua, u0) = (2.2398, 2.0, 3.25)


    descr['NN'] = "# Pandha {0} potential for same species".format(pot)
    V['NN'] = V0 * np.exp(-u0 * r['NN']) / r['NN']
    F['NN'] = V0 * np.exp(-u0 * r['NN']) / (r['NN'])**2 * (u0 * r['NN'] + 1)

    descr['NP'] = "# Pandha {0} potential for different species".format(pot)
    V['NP'] = Vr * np.exp(-ur * r['NP']) / r['NP'] -\
              Va * np.exp(-ua * r['NP']) / r['NP']
    F['NP'] = Vr * np.exp(-ur * r['NP']) / (r['NP'])**2 * (ur * r['NP'] + 1) -\
              Va * np.exp(-ua * r['NP']) / (r['NP'])**2 * (ua * r['NP'] + 1)

    descr['PP'] = "# Pandha {0} potential for same species \
    with Coulomb interaction lambda = {1}".format(pot, l)
    V['PP'] = V0 * np.exp(-u0 * r['PP']) / r['PP'] +\
              Vc * np.exp(-uc * r['PP']) / r['PP']
    F['PP'] = V0 * np.exp(-u0 * r['PP']) / (r['PP'])**2 * (u0 * r['PP'] + 1) +\
              Vc * np.exp(-uc * r['PP']) / (r['PP'])**2 * (uc * r['PP'] + 1)

  elif pot == "horowitz":
    raise AttributeError("Horowitz potential not yet implemented! Sorry :(")

  else:
    raise AttributeError("Option {0} for potential not found".format(pot))

  with open("potential.table", 'w') as fp:
    for p in pairs:
      print>>fp, descr[p]+"\n"
      print>>fp, p
      print>>fp, "N {0}\n".format(N)
      for i in xrange(N):
        print>>fp, i+1, r[p][i], V[p][i], F[p][i]



class MDSys(object):
  def __init__(self, gpu=False, silent=True, root='./data', mste=True):
    """
    Constructor: Instantiate the lammps class, so the system is always
    aware of the object it has.

    Pass dummy info to lammps and set root directory

    gpu = Use LAMMPS gpu package
    silent = Run in silent mode (no output on screen)
    root = Root directory for file hierarchy
    """

    _args = []

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
    )

    for cmd in script:
      self.lmp.command(cmd)

    if mste:
      self.lmp.command("compute mste all mste/atom 5.4")

    self.root = root
    self.parameters = {}
    self.collective = {}
    self.npairs = 4
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
      command = (('change_box all x final 0 {s} '
                  'y final 0 {s} '
                  'z final 0 {s} remap').format(s=_size),)

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

  def read_dump(self, fname):
    """
    Read information from the last snapshot of a dump file
    """

    _t = None
    with open(fname, 'r') as fp:
      for line in fp:
        if line == "ITEM: TIMESTEP":
          _t = fp.next()[:-1]

    if not _t:
      raise IOError("File {0} does not look correct.")

    cmd = ('read_dump {0} {1} x y z vx vy vz purge yes'
           'add yes replace no'.format(fname, _t))

    self.lmp.command(cmd)

  def run(self, Nsteps):
    """
    Wrapper for the "run" command in lammps
    """
    self.lmp.command("run {ns}".format(ns=Nsteps))


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
    while True:
      i = i + 1
      self.run(nfreq)
      # Extract thermo values
      [temp, ke, epair, etot, press] = self.thermo()
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

  def rdf(self, nbin, rmax):
    """
    Wrapper to the analysis rdf method
    """
    return A.rdf(self.lmp.lmp, nbin, rmax, self.npairs)

  def minkowski(self, rpart, rcell):
    """
    Wrapper to the analysis minkowski method
    """
    return A.minkowski(self.lmp.lmp, rpart, rcell)[:]

  def mste(self):
    """
    This is simply a wrapper to create a compute mste and extract
    it. LAMMPS returns simply an array with the cluster ID of each
    tag, so inside we do some calculations to return the mass
    distribution. The compute name *is* and *must be* mste [although
    this is meant to be obscured to the user in the setup() method].
    There is also a problem with setting the cutoff radius. It cannot
    be set from within the mste method, but it isn't very important: it
    should always be the cutoff of the pandha potential. Anyway, this is
    very unstable, so proceed with care when handling.

    There is a resize with mste, to make sure we can add different
    arrays afterwards. We add a lot of zeros and lose the sparsity
    of the mste, but early optimization...
    """
    _N = self.parameters['N']
    ext = self.lmp.extract_compute("mste", 1, 1)
    tmp = np.fromiter(ext, dtype=np.int, count=_N)
    # First bincount to count repeated indices => cluster size
    clust = np.bincount(tmp)
    # Filter out clusters with 0 mass
    clust = clust[clust > 0]
    mean = np.mean(clust)
    std = np.std(clust)
    # Second to histogram over sizes
    mste = np.bincount(clust)
    mste.resize(_N + 1)
    return mste, mean, std

  def lind(self):
    print "I'm inside lind and my path is {0}".format(self.this_path)
    raise AttributeError("Don't know what to do with Lindemann :(")

  def thermo(self):
    """
    Wrapper to LAMMPS internal computes.
    To avoid adding unnecesary computes to LAMMPS, we just reference
    to the default computes created for the LAMMPS inner thermo output.

    We have an advantage here: every time LAMMPS ends a run,
    calculates again thermo_temp, etc if they are in the thermo_style
    """

    temp = self.lmp.extract_compute("thermo_temp", 0, 0)
    epair = self.lmp.extract_compute("thermo_pe", 0, 0)/self.parameters['N']
    press = self.lmp.extract_compute("thermo_press", 0, 0)
    ke = 3.0/2.0 * temp
    etot = epair + ke
    return temp, ke, epair, etot, press

  def structure(self, rdf):
    """
    Calculate structure factor given the radial distribution function.
    Returns an array structured like rdf, with only one column per
    pair.

    Also return the first peak of the neutron-neutron scattering.
    """
    _d = self.parameters['d']
    r = rdf[:, 0]
    # Assume evenly spaced
    dr = r[1] - r[0]
    # Wave vectors
    n = len(r)
    q = np.linspace(0, 2*np.pi/dr, n)
    S = np.zeros((n, self.npairs+1))
    S[:, 0] = q
    for i in range(self.npairs):
      #Integrand in the fourier transform
      ker = (rdf[:, 2*i + 1] - 1) * r
      #Imaginary (sin) part of the Fourier transform
      ft = np.imag(np.fft.fft(ker)) * dr
      #Structure factor
      #We split the q = 0 case, since it is ill-defined
      S[0, i+1] = 1
      S[1:, i+1] = 1 - (ft[1:] / q[1:]) * (4 * np.pi * _d)

    data = S[:, 4]
    try:
      _change = np.sign(np.diff(data))
      c = np.diff(_change) < 0
      c = c.nonzero()[0][0] + 1 # first local max
      kmax = q[c]
      Smax = S[c, 4]
    except IndexError:
      kmax = float('nan')
      Smax = float('nan')
    return S, kmax, Smax

  def fit(self, rdf):
    """
    Fit the very-long correlations of neutron-neutron (near 15 fm)
    with a sine and returns the height of the fit and the wavelength
    """
    g = rdf[:, 7]
    r = rdf[:, 0]
    #Filter out initial zeros
    for i in range(len(g)):
      if g[i] > 0: break
      init = i
    g = g[init:]
    r = r[init:]
    sig = (g-1) * r
    sm_r, sm_sig = self.smooth(sig, r, dr=3.0)
    func = lambda x, A, l, ph: A * np.sin((2.0*np.pi*x) / l + ph)
    try:
      sm_popt, sm_pcov = curve_fit(func, sm_r, sm_sig, p0=[3.0, 15.0, 0.0])
      popt, pcov = curve_fit(func, r, sig, p0=sm_popt)
    except RuntimeError:
      return float("nan"), float("nan"), float("nan"), float("nan")
    tval = t.ppf(1.0-0.34/2, 3)
    h = sm_popt[0]
    dh = sm_pcov[0][0]**0.5 * tval
    l = sm_popt[1]
    dl = sm_pcov[1][1]**0.5 * tval
    return h, dh, l, dl

  def smooth(self, sig, r, dr=3.0):
    """
    Smooth signals
    """
    win_size = int(dr/(r[1] - r[0]))
    if win_size % 2 == 0: 
      win_size += 1
    init = (win_size -1)/2
    fin_size = len(r) - win_size
    sm_sig = np.zeros(fin_size)
    sm_r = np.zeros(fin_size)
    for i in range(fin_size):
      for j in range(win_size):
        sm_sig[i] += sig[i+j-init]
      sm_sig[i] /= win_size
      sm_r[i] = r[i+init]
    return sm_r, sm_sig


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
      t_r = self.rdf(nbins, rmax)
      self.c_rdf *= (n-1)/n
      self.c_rdf += t_r/n

    if "ssf" in self.computes:
      [t_s, a, b] = self.structure(t_r)
      self.c_ssf *= (n-1)/n
      self.c_ssf += t_s/n
      self.collective['k_absorption'] = a
      self.collective['S_absorption'] = b

    if "mste" in self.computes:
      [t_c, a, b] = self.mste()
      self.c_mste *= (n-1)/n
      self.c_mste += t_c/n
      self.collective['size_avg'] = a
      self.collective['size_std'] = b

    if "mink" in self.computes:
      [a, b, c, d] = self.minkowski(r_mink, r_cell)
      self.collective['volume'] = a
      self.collective['surface'] = b
      self.collective['breadth'] = c
      self.collective['euler'] = d


    if "thermo" in self.computes:
      [a, b, c, d, e] = self.thermo()
      self.collective['temperature'] = a
      self.collective['kinetic'] = b
      self.collective['potential'] = c
      self.collective['energy'] = d
      self.collective['pressure'] = e

    if "fit" in self.computes:
      [a, b, c, d] = self.fit(t_r)
      self.collective['height'] = a
      self.collective['del_height'] = b
      self.collective['lambda'] = c
      self.collective['del_lambda'] = d

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
    self.lmp.command("run 0 post no")
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
