import random as R
from numpy import linspace, exp

class MDSys(object):
  def __init__(self, T, l, N, x, d, V):
    """
    Constructor: so far we need to pass all the variables. Eventually
    some default values can be discussed, and some of these variables
    could be optional arguments. Some sanity checks could be here as
    well, so T is always passed as float, V as string and so on (maybe
    as a lambda function in a future?).
    """
    
    self.T = T
    self.l = l
    self.N = N
    self.x = x
    self.d = d
    self.V = V
    self.build_table(V, l)

  def build_table(self, V_tag, l, N=5000, rc_nuc=5.4, rc_cou=20.0,
		  fname="potential.table"):
    """
    This method builds the actual potential table that will be read
    in lammps. So far, three different potentials are going to be 
    supported:

    - Pandha medium
    - Pandha stiff (not yet implemented)
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
	raise ValueError("Stiff potential not yet implemented! Sorry :(")
      
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
      V['PP'] = V0 * exp(-u0 * r['PP']) / r['PP'] -\
	        Vc * exp(-uc * r['PP']) / r['PP']
      F['PP'] = V0 * exp(-u0 * r['PP']) / (r['PP'])**2 * (u0 * r['PP'] + 1) -\
                Vc * exp(-uc * r['PP']) / (r['PP'])**2 * (uc * r['PP'] + 1)
	        

    elif V_tag == "horowitz":
      raise ValueError("Horowitz potential not yet implemented! Sorry :(")

    else:
      raise ValueError("Option {0} for potential not found".format(V_tag))

    with open(self.table_fname, 'w') as fp:
      for p in pairs:
	print>>fp, descr[p]+"\n"
	print>>fp, p
	print>>fp, "N {0}\n".format(N)
	for i in xrange(N):
	  print>>fp, i+1, r[p][i], V[p][i], F[p][i]

  def build_script(self, fname = "lammps.inp"):
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

    After this, to make a proper run, the dump frequency should be 
    manually set *outside* the molecular dynamics system.
    
    As a file so we can properly debug if needed. We can also
    check if some kwargs should be added to specify when we want a
    lattice or a given data/restart/dump file, but since this seems
    outdated and hard to extend, it is not implemented yet.
    """

    self.input_fname = fname
    inp = """#Nuclear model
units		lj
atom_style	atomic
timestep	0.10

region		box block 0 {size} 0 {size} 0 {size}
create_box	2 box
create_atoms	1 random {nprot} {seed1} box
create_atoms	2 random {nneut} {seed2} box
mass		1 938.0
mass		2 938.0
velocity	all create {T} {seed3}

pair_style	table linear {ninter}
pair_coeff	1 1 {table_fname} NN 5.4
pair_coeff	1 2 {table_fname} NP 5.4
pair_coeff	2 2 {table_fname} NN 5.4

neighbor	1.2 bin
neigh_modify	every 1 delay 0 check yes one 8000 page 80000

thermo_style	custom step temp ke epair etotal press
thermo		100

min_style	hftn
minimize	0 1.0 1000 10000

thermo		1
pair_coeff	1 1 {table_fname} PP {cutoff}
fix		1 all nvt temp {T} {T} {tdamp}
""".format(size=(self.N/self.d)**(1.0/3.0),
	   nprot=int(self.x * self.N),
	   nneut=self.N - int(self.x * self.N),
	   T=self.T,
	   ninter=5000,
	   table_fname=self.table_fname,
	   cutoff=max(5.4,self.l),
	   tdamp=10.0,
	   seed1=R.randint(0,10000),
	   seed2=R.randint(0,10000),
	   seed3=R.randint(0,10000),
	   )
    with open(self.input_fname, 'w') as fp:
      print>>fp, inp

  def setup(self, lmp, path='./data', ndump='1000', thermo='1000'):
    """
    This method sets up the run in a specific lmp object according to
    the inputfile.

    Also sets path for the log and the dump file, according to the
    chosen parameters, with default values for ndump and thermo frequency.
    """
    with open(self.input_fname) as fp:
      lines = fp.readlines()
      for line in lines:
	lmp.command(line)
    
    this_path = path + "{V}/l{l}/x{x}/N{N}/d{d}/T{T}".format(V=self.V,
							l=self.l,
							x=self.x,
							N=self.N,
							d=self.d,
							T=self.T,
							)
						
    lmp.command("log {d}/thermo.log".format(d=this_path))
    lmp.command("dump 1 all custom {ndump} {d}/dump.lammpstraj\
	         type id x y z vx vy vz".format(d=this_path))

  def set_T(self, T):
    self.T = T

  def set_l(self, l):
    self.l = l

  def set_N(self, N):
    self.N = N

  def set_x(self, x):
    self.x = x

  def set_d(self, d):
    self.d = d

  def set_V(self, V):
    self.V = V
