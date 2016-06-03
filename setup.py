from setuptools import setup, find_packages
import os
import sys
#This is a list of files to install, and where
#(relative to the 'root' dir, where setup.py is)
#You could be more specific.

shr_libraries = []
for root, _, files in os.walk('.'):
  for f in files:
    if f.endswith('.so'):
      shr_libraries.append(f)

packages = find_packages()
long_description = """
The idea of this project is to create a good framework to call
LAMMPS from Python, mostly based on the Neutron Star typical
simulations. The core idea is that we have a 6 dimensional parameters
space, which are:

- V: Type of potential (Medium, Stiff, Horowitz)

- l: lambda, screening length of the Coulomb-Debye interaction

- N: Number of particles

- x: fraction of protons

- d: number density of particles

- T: temperature

It will have *no* input file, at least in its core. Python flexibility
allows an easy extension to support this behaviour, but the idea is to
create python files that instantiate the classes that will be
developed in this code."""

setup(name="LAMMPS Neutron Star Wrapper",
      version="0.1",
      description="LAMMPS from Python in Neutron Star simulations",
      author="Pablo Alcain",
      author_email="pabloalcain@gmail.com",
      url="none yet (soon!)",
      packages=packages,
      package_data={'' : shr_libraries},
      include_package_data=True,
      scripts=["examples/neutronstar.py"],
      long_description=long_description
)
