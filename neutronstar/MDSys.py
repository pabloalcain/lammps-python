import random as R
import numpy as np
import pylab as pl
from lammps import lammps
from os import makedirs
import analysis as A

class MDSys(object):
  def __init__(self, T, l, N, x, d, V, gpu=False):
    """
    Constructor: so far we need to pass all the variables. Eventually
    some default values can be discussed, and some of these variables
    could be optional arguments. Some sanity checks could be here as
    well, so T is always passed as float, V as string and so on (maybe
    as a lambda function in a future?). We also instantiate the lammps
    class here, so the system is always aware of the object it has.
    """
    
    self.T = T
    self.l = l
    self.N = N
    self.x = x
    self.d = d
    self.V = V
    self.gpu = gpu
    self.build_table(V, l)
    self.lmp = lammps("")
    #Number of pairs to be considered in the interaction: 
    #*v*, 1v1, 1v2, 2v2
    self.npairs = 4
    self.init_variables()
    
  def init_variables(self):
    self.n_tally = 0
    self.c_mste = 0
    self.c_rdf = 0
    self.c_ssf = 0
    self.computes = {}
    self.variables = []

  def build_table(self, V_tag, l, N=5000, rc_nuc=5.4, rc_cou=20.0,
		  fname="potential.table"):
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
    
    self.table_fname = fname
    pairs = ['NN', 'NP', 'PP']
    r = {}
    V = {}
    F = {}
    descr = {}
    r['NN'] = np.linspace(0,rc_nuc,N+1)[1:]
    r['NP'] = np.linspace(0,rc_nuc,N+1)[1:]
    r['PP'] = np.linspace(0,rc_cou,N+1)[1:]

    Vc = 1.44
    uc = 1.0/l

    if V_tag == "medium" or V_tag == "stiff":
      if V_tag == "medium":
	V0=373.118
	u0=1.5

	Va=2666.647
	ua=1.6

	Vr=3088.118
	ur=1.7468

      if V_tag == "stiff":
	V0=17630.256
	u0=3.25

	Va=2834.338
	ua=2.0

	Vr=3601.482
	ur=2.2395
      
      descr['NN'] = "# Pandha {0} potential for same species".format(V_tag)
      V['NN'] = V0 * np.exp(-u0 * r['NN']) / r['NN']
      F['NN'] = V0 * np.exp(-u0 * r['NN']) / (r['NN'])**2 * (u0 * r['NN'] + 1)

      descr['NP'] = "# Pandha {0} potential for different species".format(V_tag)
      V['NP'] = Vr * np.exp(-ur * r['NP']) / r['NP'] -\
		Va * np.exp(-ua * r['NP']) / r['NP']
      F['NP'] = Vr * np.exp(-ur * r['NP']) / (r['NP'])**2 * (ur * r['NP'] + 1) -\
		Va * np.exp(-ua * r['NP']) / (r['NP'])**2 * (ua * r['NP'] + 1)

      descr['PP'] = "# Pandha {0} potential for same species \
  with Coulomb interaction lambda = {1}".format(V_tag, l)
      V['PP'] = V0 * np.exp(-u0 * r['PP']) / r['PP'] +\
		Vc * np.exp(-uc * r['PP']) / r['PP']
      F['PP'] = V0 * np.exp(-u0 * r['PP']) / (r['PP'])**2 * (u0 * r['PP'] + 1) +\
		Vc * np.exp(-uc * r['PP']) / (r['PP'])**2 * (uc * r['PP'] + 1)
		

    elif V_tag == "horowitz":
      raise AttributeError("Horowitz potential not yet implemented! Sorry :(")

    else:
      raise AttributeError("Option {0} for potential not found".format(V_tag))

    with open(self.table_fname, 'w') as fp:
      for p in pairs:
	print>>fp, descr[p]+"\n"
	print>>fp, p
	print>>fp, "N {0}\n".format(N)
	for i in xrange(N):
	  print>>fp, i+1, r[p][i], V[p][i], F[p][i]

  def build_script(self, fname = "lammps.inp", dump = None):
    """
    We get all the information from this object and build an input
    file. This input file does the following:
    
    - insert masses and number of particles
    - create velocities from a Gaussian at the set temperature
    - set the interaction according to the table
    - minimize system to remove the potential energy
    - set the fix nvt at the required temperature
    - set thermo style and frequency

    Thermo style and frequency is just for following up in the screen.

    We can set the initial position and velocities from a dump file
    as well, via the optional argument dump.

    As a file so we can properly debug if needed. We can also
    check if some kwargs should be added to specify when we want a
    lattice, but since this seems outdated and hard to extend, it is
    not implemented yet.

    read_dump does NOT support changing the number of particles!
    """

    self.input_fname = fname
    if self.gpu:
      package = "package    gpu force/neigh 0 1 -1"
      style = "table/gpu"
    else: 
      package = ""
      style = "table"
      
    #The read_dump can override previous configurations, so if we
    #are meant to read a dump_file, we don't minimize (takes some time)

    if (dump == None):
	config = "minimize    0 1.0 1000 100000"
    else:
	config = "read_dump {0} 0 x y z vx vy vz purge yes add yes replace no".format(dump)


    inp = """#Nuclear model
units		lj
{package}
atom_style	atomic
timestep	0.10

region		box block 0 {size} 0 {size} 0 {size}
create_box	2 box
create_atoms	1 random {nprot} {seed1} box
create_atoms	2 random {nneut} {seed2} box
mass		1 938.0
mass		2 938.0
velocity	all create {T} {seed3}

pair_style	{style} linear {ninter}
pair_coeff	1 1 {table_fname} NN 5.4
pair_coeff	1 2 {table_fname} NP 5.4
pair_coeff	2 2 {table_fname} NN 5.4

neighbor	1.2 bin
neigh_modify	every 1 delay 0 check yes one 8000 page 80000

thermo_style	custom step temp ke epair etotal press
thermo		1000

min_style	hftn

dump            1 all custom 1000 minim.lammpstrj type id x y z vx vy vz 
{config}
undump          1

pair_coeff	1 1 {table_fname} PP {cutoff}
fix		1 all nvt temp {T} {T} {tdamp}

reset_timestep  0
  """.format(package=package,
	   style=style,
	   size=(self.N/self.d)**(1.0/3.0),
	   nprot=int(self.x * self.N),
	   nneut=self.N - int(self.x * self.N),
	   T=self.T,
	   ninter=5000,
	   table_fname=self.table_fname,
	   cutoff=max(5.4,self.l),
	   tdamp=10.0,
	   config=config,
	   seed1=R.randint(0,10000),
	   seed2=R.randint(0,10000),
	   seed3=R.randint(0,10000),
	   )
    with open(self.input_fname, 'w') as fp:
      print>>fp, inp

  def setup(self, path='./data',
	    rdf = True, ssf = True,
	    mink = True, mste = True,
	    lind = True, thermo = True):
    """
    This method sets up the run in a specific lmp object according to
    the inputfile. We also set here the computes.
    """
    if ssf and not rdf:
      raise AttributeError("Cannot calculate structure factor without rdf")
    
    with open(self.input_fname) as fp:
      lines = fp.readlines()
      for line in lines:
	self.lmp.command(line)
    
    self.path = path
    self.is_rdf = rdf
    self.is_ssf = ssf
    self.is_mink = mink
    self.is_mste = mste
    self.is_lind = lind
    self.is_thermo = thermo
     
  def run(self, Nsteps):
    """
    Wrapper for the "run" command in lammps
    """
    self.lmp.command("run {ns}".format(ns=Nsteps))

  def set_T(self, T, tdamp = 10.0):
    """
    Set value of temperature and change all related values
    in the lammps script.
    """
    self.T = T
    self.update_path()
    self.lmp.command("unfix 1")
    temp = "fix 1 all nvt temp {T} {T} {tdamp}"
    self.lmp.command(temp.format(T=self.T,tdamp=10.0))

  def set_l(self, l):
    """
    Set value of lambda and change all related values in the lammps
    script.
    """
    self.l = l
    self.update_path()
    self.build_table(self.V, self.l, fname = self.table_fname)
    cmd = "pair_coeff 1 1 {t} PP {co}".format(t=self.table_fname,
					      co=max(5.4,self.l),
					      )
    self.lmp.command(cmd)

  def set_V(self, V):
    """
    Set value of the potential and change all related values in the
    lammps script.
    """
    self.V = V
    self.update_path()
    self.build_table(self.V, self.l, fname = self.table_fname)
    cmd = """pair_coeff 1 1 {t} PP {co}
    pair_coeff	1 2 {t} NP 5.4
    pair_coeff	2 2 {t} NN 5.4
    """.format(t=self.table_fname,
	       co=max(5.4,self.l),
	       )
    self.lmp.command(cmd)

  def set_N(self, N):
    """
    This should add particles randomly inside the already given 
    configuration. Not implemented yet
    """
    self.update_path()
    self.N = N
    raise AttributeError("This is not implemented yet :(")

  def set_x(self, x):
    """
    This should change randomly particles type to fulfill the x
    requirement. Not implemented yet
    """
    self.update_path()
    self.x = x
    raise AttributeError("This is not implemented yet :(")

  def set_d(self, d):
    """
    This method changes the box size in order to get the density
    requirement. The particles are remapped, so their relative
    position change.
    """
    self.update_path()
    self.d = d
    cmd = ("change_box all "
	   "x final 0 {size} "
	   "y final 0 {size} "
	   "z final 0 {size} "
	   "remap").format(size=(self.N/self.d)**(1.0/3.0))

  def update_path(self):
    """
    Update this_path according to value of parameters
    """
    self.this_path = self.path + "/{V}/l{l}/x{x}/N{N}/d{d}/T{T}/"
    self.this_path = self.this_path.format(V=self.V, l=self.l, x=self.x,
					   N=self.N, d=self.d, T=self.T)
    try:
      makedirs(self.this_path)
    except OSError as err:
      msg = "Directory {0} already exists: rename base path or delete old files"
      raise OSError(msg.format(self.this_path))
    
  def equilibrate(self, tdamp = 10.0, nfreq = 100, wind = 100):
    """
    This method takes care of the thermalization, with a berendsen
    thermostat. It has a quite unstable part, in which it unfixes a
    previous "fix 1 all nvt" lammps command. This works with no
    problem, because the script building already sets a temperature,
    but be VERY careful with respect to this.

    nfreq is how many timesteps to take between runs, and wind is the
    size of the window. The criterion for stability is that the
    average temperature of the last 100 steps is close to the set
    temperature by a standard deviation, while the energy stops
    decreasing (slope > 0). This works only when setting temperature
    from high to low (while for going from low to high, one would
    expect that just setting the temperatures should be enough).

    Thermalization doesn't write any log or dump file.
    """
    self.lmp.command("unfix 1")
    brd = "fix 1 all temp/berendsen {T} {T} {tdamp}"
    self.lmp.command(brd.format(T = self.T, tdamp = tdamp))
    self.lmp.command("fix 2 all nve")
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
      if i < wind: continue
      # Slow, we are calculating things again. Could be updated
      # and deleted each time. Check profiling!
      [slope, aux] = np.polyfit(step, energy, 1)
      diff = abs(self.T - np.mean(temperature))
      std = np.std(temperature)
      if slope > 0 and diff < std: break
    self.lmp.command("unfix 2")
    
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
    ext = self.lmp.extract_compute("mste", 1, 1)
    tmp = np.fromiter(ext, dtype = np.int, count = self.N)
    # First bincount to count repeated indices => cluster size
    clust = np.bincount(tmp)
    # Filter out clusters with 0 mass
    clust = clust[clust > 0] 
    mean = np.mean(clust)
    std = np.std(clust)
    # Second to histogram over sizes
    mste = np.bincount(clust)
    mste.resize(self.N + 1) 
    return mste, mean, std
    
  def lind(self):
    raise AttributeError("Don't know what to do with Lindemann :(")
    print "I'm inside lind and my path is {0}".format(self.this_path)
    pass

  def thermo(self):
    """
    Wrapper to LAMMPS internal computes.
    To avoid adding unnecesary computes to LAMMPS, we just reference
    to the default computes created for the LAMMPS inner thermo output.

    We have an advantage here: every time LAMMPS ends a run,
    calculates again thermo_temp, etc if they are in the thermo_style
    """

    temp = self.lmp.extract_compute("thermo_temp", 0, 0)
    epair = self.lmp.extract_compute("thermo_pe", 0, 0)/self.N
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
    r = rdf[:,0]
    # Assume evenly spaced
    dr = r[1] - r[0]
    # Wave vectors
    n = len(r)
    q = np.linspace(0,2*np.pi/dr,n)
    S = np.zeros((n,self.npairs+1))
    S[:,0] = q
    for i in range(self.npairs):
      #Integrand in the fourier transform
      ker = (rdf[:, 2*i + 1] - 1) * r
      #Imaginary (sin) part of the Fourier transform
      ft = np.imag(np.fft.fft(ker)) * dr
      #Structure factor
      #We split the q = 0 case, since it is ill-defined
      S[0,i+1] = 1
      S[1:,i+1] = 1 - ( ft[1:] / q[1:] ) * ( 4 * np.pi * self.d )
      
    data = S[1:,4]
    try: 
      c = (np.diff(np.sign(np.diff(data))) < 0).nonzero()[0][0] + 1 # first local max
      kmax = q[c]
      Smax = S[c, 4]
    except IndexError:
      kmax = float('nan')
      Smax = float('nan')
    return S, kmax, Smax

  def results(self, 
	      r_mink=1.8, r_cell=0.5,
	      nbins = 200, rmax = None):
    """
    Method to take all the results that have been set in the setup()
    method
    """
    
    if rmax == None: rmax = (float(self.N)/self.d)**(1.0/3)*0.5
    
    self.n_tally+=1
    n = float(self.n_tally)
    
    if self.is_rdf:
      t_r = self.rdf(nbins, rmax)
      self.c_rdf *= (n-1)/n
      self.c_rdf += t_r/n

    if self.is_ssf:
      [t_s, a, b] = self.structure(t_r)
      self.c_ssf *= (n-1)/n
      self.c_ssf += t_s/n
      self.computes['k_absorption'] = a
      self.computes['S_absorption'] = b

    if self.is_mste:
      [t_c, a, b] = self.mste()
      self.c_mste *= (n-1)/n
      self.c_mste += t_c/n
      self.computes['size_avg'] = a
      self.computes['size_std'] = b

    if self.is_mink:
      [a, b, c, d] = self.minkowski(r_mink, r_cell)
      self.computes['volume'] = a
      self.computes['surface'] = b
      self.computes['breadth'] = c
      self.computes['euler'] = d
      

    if self.is_thermo:
      [a, b, c, d, e] = self.thermo()
      self.computes['temperature'] = a
      self.computes['kinetic'] = b
      self.computes['potential'] = c
      self.computes['energy'] = d
      self.computes['pressure'] = e

    if self.is_lind:
      self.lind()

  def dump(self, prefix = ''):
    """
    Wrapper to dump positions
    """
    path = self.this_path + prefix
    dump_fname = path + 'dump.lammpstrj'
    tmp = "dump myDUMP all custom 1 {dumpfile} id x y z vx vy vz"
    self.lmp.command(tmp.format(dumpfile = dump_fname))
    self.lmp.command("dump_modify myDUMP sort id")
    self.lmp.command("run 0 post no")
    self.lmp.command("undump myDUMP")

  def flush(self, prefix = ''):
    """
    Write log in the data file and plot
    """
    path = self.this_path + prefix
    if self.is_mste:
      # To add cluster size to file
      temp = np.vstack((range(len(self.c_mste)),self.c_mste))
      mste_fname = path + 'cluster.dat'
      np.savetxt(mste_fname, temp.T, header = 'size, number', fmt='%6i + %1.4e')
      
    if self.is_rdf:
      rdf_fname = path + 'rdf.dat'
      h = 'r, a-a, ia-a, 1-1, i1-1, 1-2, i1-2, 2-2, i2-2'
      np.savetxt(rdf_fname, self.c_rdf, header = h, fmt = '%1.4e')

    if self.is_ssf:
      ssf_fname = path + 'ssf.dat'
      h = 'r, a-a, 1-1, 1-2, 2-2'
      np.savetxt(ssf_fname, self.c_ssf, header = h)

    if self.is_thermo:
      thermo_fname = path + 'thermo.dat'
      h = ', '.join(self.computes.keys())
      np.savetxt(thermo_fname, self.variables, header = h, fmt = '%1.4e')
    
    self.init_variables()
    self.lmp.command("reset_timestep 0")

  def log(self):
    """
    Append new data to variables. This has a problem with memory
    usage, and the fact that append may render slow for big
    arrays. Should be checked afterwards for long runs.

    Thoughts: log might upgrade the mean of each column + std dev on
    the fly.
    """
    self.variables.append(self.computes.values())
