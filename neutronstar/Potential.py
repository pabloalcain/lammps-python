"""
Module for the calculation of the potential energy in Neutron Star
systems.
"""
import numpy as np
import random as rn

def _interaction(pot="medium", lda=20):
  """
  Calculate interaction function for a nuclear potential.

  Parameters
  ----------

  pot : {"medium", "stiff", "newmed", "horowitz"}
      Potential style to calculate

  lda : float
      Debye screening length

  Returns
  -------

  potfunc : dict
      A dictionary with the function of the potential for each pair
  """
  Vc = 1.44
  uc = 1.0/lda
  if pot in ["medium", "stiff", "newmed"]:
    if pot == "medium":
      (Vr, Va, V0) = (3088.118, 2666.647, 373.118)
      (ur, ua, u0) = (1.7468, 1.6, 1.5)

    if pot == "newmed": # This doesn't look correct
      (Vr, Va, V0) = (3097.0, 2696.0, 379.5)
      (ur, ua, u0) = (1.648, 1.528, 1.628)

    if pot == "stiff":
      (Vr, Va, V0) = (3601.482, 2834.338, 17630.256)
      (ur, ua, u0) = (2.2395, 2.0, 3.25)

    vnn = lambda r: V0 * np.exp(-u0 * r) / r
    vnp = lambda r: Vr * np.exp(-ur * r) / r - Va * np.exp(-ua * r) / r
    vpp = lambda r: vnn(r) + Vc * np.exp(-uc * r) / r

  elif pot == "horowitz":
    a = 110.0
    b = -26.0
    c = 24.0
    L = 1.25
    vnn = lambda r: a * np.exp(-r**2 / L) + (b + c) * np.exp(-r**2 / (2*L))
    vnp = lambda r: a * np.exp(-r**2 / L) + (b - c) * np.exp(-r**2 / (2*L))
    vpp = lambda r: vnn(r) + Vc * np.exp(-uc * r) / r

  else:
    raise AttributeError("Option {0} for potential not found".format(pot))

  return {'NN': vnn, 'NP': vnp, 'PP': vpp}

def build_table(pot, lda, name=None):
  """
  This method builds the actual potential table that will be read
  in lammps. So far, three different potentials are going to be
  supported:

  - Pandha medium
  - Pandha stiff
  - Horowitz
  """

  if name == None:
    name = '{0}-{1}-{2}.table'.format(pot, lda, rn.randint(1, 1000))
  rc_nuc = 5.4
  npoints = 5000
  potential = {}
  descr = {}
  positions = {}
  positions['NN'] = np.linspace(0, rc_nuc, npoints+1)[1:]
  positions['NP'] = np.linspace(0, rc_nuc, npoints+1)[1:]
  positions['PP'] = np.linspace(0, max(rc_nuc, lda), npoints+1)[1:]

  potfunction = _interaction(pot, lda)
  for i in positions.keys():
    potential[i] = potfunction[i](positions[i])

  descr['NN'] = "# {0} potential for same species".format(pot)
  descr['NP'] = "# {0} potential for different species".format(pot)
  descr['PP'] = ("# {0} potential for same species "
                 "with Coulomb interaction lambda = {1}").format(pot, lda)


  with open(name, 'w') as potfile:
    for pair in positions.keys():
      print>>potfile, descr[pair]+"\n"
      print>>potfile, pair
      print>>potfile, "N {0}\n".format(npoints)
      _pot = potential[pair]
      _pos = positions[pair]
      _force = - np.diff(_pot) / np.diff(_pos)
      _force.resize(npoints)
      for i, info in enumerate(zip(_pos, _pot, _force)):
        print>>potfile, i+1, info[0], info[1], info[2]

  return name
