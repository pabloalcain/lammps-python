"""
Test refactor system. The way we want it to work, at least
"""
import pylammps as plmp
import numpy as np

sys = plmp.NeutronStarSystem(silent=False)

k = np.linspace(0.2, 0.5, 200)
pairs = (((0,), (0,),),)
computes = {'rdf': plmp.Computes.RDF(200, pairs),
            'ssf': plmp.Computes.StructureFactor(k, pairs),
            'thermo': plmp.Computes.Thermo(),
            'mste': plmp.Computes.MSTE()}
analyzer = plmp.Analyzer(computes)

sys['potential'] = 'medium'
sys['lambda'] = 20
sys['x'] = 0.5
sys['N'] = 50
sys['T'] = 2.0
sys['d'] = 0.05
sys.minimize(0, 1.0, 1000, 100000)
for T in np.linspace(2.0, 0.5, 31):
  sys['T'] = T
  log = plmp.Logger(sys)
  log.dump(sys, 'image')
  sys.thermalize(200, 30)
  analyzer.zero()
  for i in range(5):
    sys.run(100)
    log.dump(sys, 'text')
    analyzer.update(sys)
  log.dump(sys, 'image')
  log.log(analyzer)
  log.plot(analyzer)
