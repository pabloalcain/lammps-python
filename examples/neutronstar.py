"""
An example that runs a LAMMPS simulation for Neutron Star Matter.
"""
import pylammps as plmp
import numpy as np

#Instantiate the Neutron Star System, that inherits from a LAMMPS
#System.
sys = plmp.NeutronStarSystem(silent=False)

#Variables we we use in the computes
k = np.linspace(0.2, 0.5, 200)
pairs = (((0,), (0,),),)

#Compute dictionary to instantiate the analyzer. See that we can add
#many times the same style of computes with differemt parameters as
#long we give them different keys in the dict.
computes = {'rdf': plmp.Computes.RDF(200, pairs),
            'ssf': plmp.Computes.StructureFactor(k, pairs),
            'thermo': plmp.Computes.Thermo(),
            'mste': plmp.Computes.MSTE()}

analyzer = plmp.Analyzer(computes)

#Set values for the system
sys['potential'] = 'medium'
sys['lambda'] = 20
sys['x'] = 0.5
sys['N'] = 50
sys['T'] = 2.0
sys['d'] = 0.05

#Remove some of the initial energy
sys.minimize(0, 1.0, 1000, 100000)
for T in np.linspace(2.0, 0.5, 31):
  sys['T'] = T
  #We instantiate the logger in every case, since want to change the
  #directory for different temperatures
  log = plmp.Logger(sys)
  log.dump(sys, 'image')
  #Thermalize the system
  sys.thermalize(200, 30)
  #Zero out every compute in the analyzer
  analyzer.zero()
  for i in range(5):
    sys.run(100)
    log.dump(sys, 'text')
    #Update the computes in the analyzer
    analyzer.update(sys)
  log.dump(sys, 'image')
  #Log and plot the results of the analyzer
  log.log(analyzer)
  log.plot(analyzer)
