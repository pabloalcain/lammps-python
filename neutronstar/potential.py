import numpy as np

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

