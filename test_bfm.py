from postprocess import extract, plotter
from analysis import ECRA, MSTE
import pylab
ext = extract.Extraction();

for i in range(0, 300, 30):
  x = ext.x('dump2.lammpstrj', i)
  v = ext.v('dump2.lammpstrj', i)
  t = ext.t('dump2.lammpstrj', i)
  box = ext.box('dump2.lammpstrj', i)
  value_mste, _ = MSTE.MSTE.mste(x, v, t, box, True, 0.0001, model='newmed')
  plotter.mste(value_mste)
  pylab.savefig('mste_{0}.pdf'.format(i))
  pylab.savefig('mste_{0}.png'.format(i))
  pylab.close()
print "Done MSTE"
for i in range(0, 300, 30):
  print i
  x = ext.x('dump2.lammpstrj', i)
  v = ext.v('dump2.lammpstrj', i)
  t = ext.t('dump2.lammpstrj', i)
  box = ext.box('dump2.lammpstrj', i)
  value, _ = ECRA.ECRA.ecra_puente(x, v, t, box, 0.0001)
  plotter.mste(value)
  pylab.savefig('puente_{0}.pdf'.format(i))
  pylab.savefig('puente_{0}.png'.format(i))
  pylab.close()
