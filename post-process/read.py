import numpy as np
import pylab as pl
import os

class PostProcess(object):
  def load_file(self):
    file = open('thermo.dat')
    self.dic = {}
    keys = file.readline().rstrip()[2:].split(', ')
    for i in range(len(keys)):
      self.dic[keys[i]] = i
    file.close()
    
    self.values = np.loadtxt('thermo.dat')
    self.mean = [np.mean(self.values[:,i]) for i in range(len(keys))]
    self.std = [np.std(self.values[:,i]) for i in range(len(keys))]
  
  def hash():
    d = os.walk('data')
    while True:
      try:
        info = d.next()
      except StopIteration:
        break
      if 'thermo.dat' in info[2]:
        print info[0]
      
