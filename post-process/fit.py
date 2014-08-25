from numpy import sin, pi, loadtxt, zeros
from scipy.optimize import curve_fit

def wavelength(fname):
  # This method returns the wavelength near lambda = 15 fm
  # and the amplitude for the usually formatted file "rdf.dat"
  # in neutron-neutron interaction
  A = loadtxt(fname)
  g = A[:,7]
  r = A[:,0]
  #Filter out initial zeros
  for i in range(len(g)):
    if g[i] > 0: break

  init = i
  g = g[init:]
  r = r[init:]
  sig = (g-1) * r
  sm_r, sm_sig = smooth(sig, r, dr = 1.8)
  func = lambda x, A, l, ph: A * sin((2.0*pi*x) / l + ph)
  try:
    sm_popt, sm_pcov = curve_fit(func, sm_r, sm_sig, p0 = [3.0, 15.0, 0.0])
    popt, pcov = curve_fit(func, r, sig, p0 = sm_popt)
  except RuntimeError:
    return [0, 0, 0], [[0, 0, 0],[0, 0, 0],[0,0,0]]
  return popt, pcov
  

def smooth(sig, r, dr = 1.8):
  # Get the signal and smooth variations of order "dr"

  win_size = int(dr/(r[1] - r[0]))
  if win_size % 2 == 0: win_size += 1
  init = (win_size -1)/2
  fin_size = len(r) - win_size
  sm_sig = zeros(fin_size)
  sm_r = zeros(fin_size)
  for i in range(fin_size):
    for j in range(win_size):
      sm_sig[i] += sig[i+j-init]
    sm_sig[i] /= win_size
    sm_r[i] = r[i+init]
  return sm_r, sm_sig
