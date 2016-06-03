"""
Build and installation script
"""
from setuptools import setup, find_packages
import os

shr_libraries = []
for root, _, files in os.walk('.'):
  for f in files:
    if f.endswith('.so'):
      shr_libraries.append(f)

packages = find_packages()
setup(name="lammps-python",
      version="0.1",
      description="A Python wrapper to LAMMPS",
      author="Pablo Alcain",
      author_email="pabloalcain@gmail.com",
      url="none yet (soon!)",
      packages=packages,
      package_data={'' : shr_libraries},
      include_package_data=True,
      scripts=["examples/neutronstar.py"],)
