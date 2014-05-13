TODO LIST
=========

Global Issues
-------------

- [ ] What is the lindemann coefficient? How can we calculate it?
- [ ] Reorder code
  * [ ] Move files and folders
  * [ ] Put python files inside a folder, maybe a whole package?
- [ ] Analysis: do we _really_ need to add those horrible methods
      inside MDSys? I think not, should be worked out a bit.


LAMMPS wrapper
--------------

- [ ] Thermo and dump frequency should be none. Files written from 
      Python itself. 
  * [X] Extract value of computes and write to file whatever
        we want, to replace thermo
  * [X] Create a wrapper for the dump LAMMPS command
- [ ] Thermalization: improve naïve Berendsen implementation, look in
      old_files example
- [ ] Throw error if dump file number of particles doesnt mach with N
- [ ] Need compute mste/atom to be shipped with the distro?
- [ ] Add makes/installs/pip/whatever
- [ ] Devise a method for bookkeeping
- [ ] PostProcess:
  * [ ] Work out clean way to work with computes within LAMMPS (MSTE, for example)

- [X] When setting temperature, add flag for thermalization
- [X] Set initial position and velocities from file
- [X] Create classes
- [X] Put table creation inside the classes
- [X] Copy philosphy from CoffeeMD

Analysis wrapper
----------------

- [ ] Add warning to minkowski method related to lattice size
- [X] Decide what to do with the npair in rdf. Should be taken
      outside from analysis.py
- [X] Add S(k) to g(r)
- [ ] Check analysis
  * [ ] MSTE
  * [ ] Minkowski
  * [ ] Lindemann
  * [ ] RDF

- [X] Clean&Check the old analysis.py routines, c-style them
- [X] Insert computes inside the actual run as a library




- Finish Horowitz potential [O]
