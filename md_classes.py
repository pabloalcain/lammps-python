class mdsys(Object):
  self.T = 1.0
  self.l = 20
  self.N = 5000
  self.x = 0.5
  self.d = 0.04
  self.V = "med"
  self.filename = "table.dat"

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
    We get all the information from this object and build the input
    file. As a file so we can properly debug if needed. We can also
    check if some kwargs should be added to specify when we want a
    lattice, but since this seems outdated and hard to extend, it is
    not implemented yet.
    """
    self.input_fname = fname
    """ BUNCH OF DIRTY FUNCTIONS """
    pass
  


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
