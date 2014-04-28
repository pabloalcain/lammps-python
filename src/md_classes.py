import random as R
from numpy import linspace, exp
from lammps import lammps
from os import makedirs, path

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
    r['NN'] = linspace(0,rc_nuc,N+1)[1:]
    r['NP'] = linspace(0,rc_nuc,N+1)[1:]
    r['PP'] = linspace(0,rc_cou,N+1)[1:]

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
      V['NN'] = V0 * exp(-u0 * r['NN']) / r['NN']
      F['NN'] = V0 * exp(-u0 * r['NN']) / (r['NN'])**2 * (u0 * r['NN'] + 1)

      descr['NP'] = "# Pandha {0} potential for different species".format(V_tag)
      V['NP'] = Vr * exp(-ur * r['NP']) / r['NP'] -\
	        Va * exp(-ua * r['NP']) / r['NP']
      F['NP'] = Vr * exp(-ur * r['NP']) / (r['NP'])**2 * (ur * r['NP'] + 1) -\
                Va * exp(-ua * r['NP']) / (r['NP'])**2 * (ua * r['NP'] + 1)

      descr['PP'] = "# Pandha {0} potential for same species \
with Coulomb interaction lambda = {1}".format(V_tag, l)
      V['PP'] = V0 * exp(-u0 * r['PP']) / r['PP'] +\
	        Vc * exp(-uc * r['PP']) / r['PP']
      F['PP'] = V0 * exp(-u0 * r['PP']) / (r['PP'])**2 * (u0 * r['PP'] + 1) +\
                Vc * exp(-uc * r['PP']) / (r['PP'])**2 * (uc * r['PP'] + 1)
	        

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
    - set thermo style and frequency (can be overridden later?)
    - set dump style

    We can set the initial position and velocities from a dump file
    as well, via the optional argument dump.

    After this, to make a proper run, the dump frequency should be 
    done in the setup method.
    
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
        config = "read_dump {0} 0 x y z vx vy vz".format(dump)


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
thermo		100

min_style	hftn

dump            1 all custom 1000 minim.lammpstrj type id x y z vx vy vz 
{config}

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

  def setup(self, path='./data', ndump=1000, nthermo=1000):
    """
    This method sets up the run in a specific lmp object according to
    the inputfile.
    """
    with open(self.input_fname) as fp:
      lines = fp.readlines()
      for line in lines:
	self.lmp.command(line)
    
    self.path = path
    self.ndump = ndump
    self.nthermo = nthermo

  def run(self, Nsteps):
    """
    Wrapper for the "run" command in lammps
    """
    self.lmp.command("run {ns}".format(ns=Nsteps))

  def set_T(self, T, tdamp = 10.0, therm = False):
    """
    Set value of temperature and change all related values
    in the lammps script.
    """
    self.T = T
    if (therm): self.thermalize()
    self.update_files()
    self.lmp.command("unfix 1")
    self.lmp.command("fix 1 all nvt temp {T} {T} {tdamp}".format(T=self.T,
                                                         tdamp=10.0))

  def set_l(self, l):
    self.l = l
    self.update_files()
    self.build_table(self.V, self.l, fname = self.table_fname)
    cmd = "pair_coeff 1 1 {t} PP {co}".format(t=self.table_fname,
                                              co=max(5.4,self.l),
                                              )
    self.lmp.command(cmd)

  def set_V(self, V):
    self.V = V
    self.update_files()
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
    self.update_files()
    self.N = N
    raise AttributeError("This is not implemented yet :(")
  
  def set_x(self, x):
    """
    This should change randomly particles type to fulfill the x
    requirement. Not implemented yet
    """
    self.update_files()
    self.x = x
    raise AttributeError("This is not implemented yet :(")

  def set_d(self, d):
    """
    This method changes the box size in order to get the density
    requirement. The particles are remapped, so their relative
    position change.
    """
    
    self.d = d
    self.update_files()
    cmd = ("change_box all "
           "x final 0 {size} "
           "y final 0 {size} "
           "z final 0 {size} "
           "remap").format(size=(self.N/self.d)**(1.0/3.0))

  def update_files(self, therm = False):
    self.lmp.command("reset_timestep 0")
    self.lmp.command("undump 1")
    self.this_path = self.path + "/{V}/l{l}/x{x}/N{N}/d{d}/T{T}".format(V=self.V,
                                                                        l=self.l,
                                                                        x=self.x,
                                                                        N=self.N,
                                                                        d=self.d,
                                                                        T=self.T,
                                                                        )
    
    if (therm): fname = "thermalization"
    else: fname = "evolution"
    
    
    try:
      makedirs(self.this_path)
    except OSError as err:
      if (path.isfile(self.this_path+'/'+fname+'.log')):
        msg = "File {0} already exists: rename base path or delete old files"
        raise OSError(msg.format(self.this_path+'/'+fname+'.log'))

    dmp = ("dump 1 all custom {ndump} "
           "{d}/{fname}.lammpstrj "
           "type id x y z vx vy vz".format(d=self.this_path,
                                           ndump=self.ndump,
                                           fname=fname
                                           )
           )

    self.lmp.command("log {d}/{fname}.log".format(d=self.this_path,
                                                  fname=fname))
    self.lmp.command("thermo {nthermo}".format(nthermo=self.nthermo))
    self.lmp.command(dmp)
    self.lmp.command("dump_modify 1 sort id")

  def thermalize(self, tdamp = 1000, nsteps = 80000):
    """
    This method takes care of the thermalization, with a berendsen
    thermostat. It has a quite unstable part, in which it unfixes a
    previous "fix 1 all nvt" lammps command. This works with no
    problem, because the script building already sets a temperature,
    but be VERY careful with respect to this. 
    """
    
    self.update_files(therm = True)
    self.lmp.command("unfix 1")
    brd = ("fix 1 all temp/berendsen"
           " {T} {T} {tdamp}".format(T=self.T,
                                     tdamp=tdamp
                                     )
           )
    self.lmp.command(brd)
    self.lmp.command("fix 2 all nve")
    self.lmp.command("run {nsteps}".format(nsteps = nsteps))
    self.lmp.command("unfix 2")
    
    
  def results(self, rdf = True, mink = True, mste = True, lind = True, thermo = True):
    if (rdf): self.rdf()
    if (mink): self.mink()
    if (mste): self.mste()
    if (lind): self.lind()
    if (thermo): self.thermo()
    
  def rdf(self):
    print "I'm inside rdf and my path is {0}".format(self.this_path)
    pass
  
  def mink(self):
    print "I'm inside mink and my path is {0}".format(self.this_path)
    pass
  
  def mste(self):
    print "I'm inside mste and my path is {0}".format(self.this_path)
    pass
  
  def lind(self):
    print "I'm inside lind and my path is {0}".format(self.this_path)
    pass
  
  def thermo(self):
    print "I'm inside thermo and my path is {0}".format(self.this_path)
    pass

