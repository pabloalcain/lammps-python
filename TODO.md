TODO LIST
=========

Global Issues
-------------

- [ ] What is the lindemann coefficient? How can we calculate it?
- [ ] Reorder code
  * [ ] Move files and folders
  * [ ] Put python files inside a folder, maybe a whole package?
  * [ ] Attention! Methods should be a bit more encapsulated, use
        advantage of OOP and have different objects, for system,
        results, etc (if needed)
- [ ] Analysis: do we *really* need to add those horrible methods
      inside MDSys? I think not, should be worked out a bit.


LAMMPS wrapper
--------------

- [ ] Write mean and std of computes in a global file inside path
- [ ] Need to work a bit about when to add the files, when to 
      update_path to further check if the files *do* exist.
      Idea: makedirs without problem, but throw error when creating
      the files in flush()?
- [ ] Thermalization: improve naïve Berendsen implementation, look in
      old_files example
- [ ] Throw error if dump file number of particles doesnt mach with N
- [ ] Need compute mste/atom to be shipped with the distro?
- [ ] Add makes/installs/pip/whatever
- [ ] Devise a method for bookkeeping

Analysis wrapper
----------------

- [ ] Plot computes
- [ ] Add warning to minkowski method related to lattice size
- [ ] Check analysis
  * [ ] MSTE
  * [ ] Minkowski
  * [ ] Lindemann
  * [ ] RDF

LAMMPS wrapper
--------------

- [X] When setting temperature, add flag for thermalization
- [X] Set initial position and velocities from file
- [X] Create classes
- [X] Put table creation inside the classes
- [X] Copy philosphy from CoffeeMD
- [X] Thermo and dump frequency should be none. Files written from 
      Python itself. 
  * [X] Extract value of computes and write to file whatever
        we want, to replace thermo
  * [X] Create a wrapper for the dump LAMMPS command
- [X] PostProcess:
  * [X] Work out clean way to work with computes within LAMMPS (MSTE, 
        for example)

Analysis wrapper
----------------

- [X] Clean&Check the old analysis.py routines, c-style them
- [X] Insert computes inside the actual run as a library
- [X] Decide what to do with the npair in rdf. Should be taken
      outside from analysis.py
- [X] Add S(k) to g(r)
- [X] Return value of peak of S(k) for low momenta


- [O] Finish Horowitz potential
