from plotme import Plotter
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as pl
import sys
from matplotlib.backends.backend_pdf import PdfPages
#data="therm_data_1.dat"

font = {'size'   : 18}
p = Plotter('/home/pablo/lammps-neutron-stars-input/examples/data_explor')

mag = sys.argv[1]
pp = PdfPages('{0}.pdf'.format(mag))

x_l = [0.4, 0.46, 0.48, 0.5]
d_l = np.linspace(0.01,0.08,8)
T_l = np.linspace(0.5, 2.0, 31)
#x_l = [0.4]
#d_l = [0.06]
for x in x_l:
  for d in d_l:
    med = []
    newmed = []
    for T in T_l:
      out = p.extract(mag, V="medium", x=x, d=d, T=T)
      med.append(out)
      out = p.extract(mag, V="newmed", x=x, d=d, T=T)
      newmed.append(out)
    med = np.array(med)
    newmed = np.array(newmed)

    pl.figure()
    ax = pl.subplot(1,1,1)
    ax.errorbar(T_l, med[:, 0], marker="o",
                yerr=med[:, 1], label="Medium Illinois")
    ax.errorbar(T_l, newmed[:, 0], marker="^",
                yerr=newmed[:, 1], label="Medium Actual")
    ax.legend(numpoints=1)
    pl.xlabel('Temperature [MeV]')
#    pl.ylabel('Absorption wavelength [fm]')
    pl.ylabel(mag)
    pl.title('x = {x}, d = {d}'.format(x=x,d=d))
    pl.savefig('{mag}-x{x}-d{d}.pdf'.format(mag=mag,x=x,d=d))
    pp.savefig()
pp.close()
