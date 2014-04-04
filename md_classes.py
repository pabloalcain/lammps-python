import random as R

class mdsys(object):
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

  def build_table(self, V, l, N=5000, fname="table.dat"):
    """
    This method builds the actual potential table that will be read
    in lammps.
    """
    self.table_fname = fname
    """ BUNCH OF DIRTY FUNCTIONS """
    pass

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

pair_style	table linear 
pair_coeff	1 1 {table_fname} NN 5.4
pair_coeff	1 2 {table_fname} NP 5.4
pair_coeff	2 2 {table_fname} NN 5.4

neighbor	1.2 bin
neigh_modify	every 1 delay 0 check yes one 8000 page 80000

thermo_style	custom step temp epair etotal press
thermo		1000

min_style	hftn
minimize	0 1.0 1000 1000000

pair_coeff	1 1 {table_fname} PP {cutoff}
fix		1 all nvt temp {T} {T} {tdamp}
""".format(size=(self.N/self.d)**(1.0/3.0),
	             nprot=self.x * self.N,
		     nneut=(1-self.x) * self.N,
		     T=self.T,
		     table_fname=self.table_fname,
		     cutoff=max(5.4,self.l),
		     tdamp=10.0,
		     seed1=R.randint(0,10000),
		     seed2=R.randint(0,10000),
		     seed3=R.randint(0,10000),
		     )
    with open(self.input_fname, 'w') as fp:
      print>>fp, inp

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
