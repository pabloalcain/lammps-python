import numpy as np
import pylab as pl
import os


class Processor(object):
  def __init__(self, path = 'data/'):
    self.path = path

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
    
  def hash(self):
    self.data = {}
    data_dir = os.walk(self.path)
    offset = len(self.path)
    while True:
      try:
        info = data_dir.next()
      except StopIteration:
        break
      if 'thermo.dat' in info[2]:
        pot, lam, x, N, d, T = info[0][offset:].split('/')
        lam, x, N, d, T = map(float, (lam[1:], x[1:], N[1:], d[1:], T[1:]))
        self.data[(pot, lam, x, N, d, T)] = info[0]+'/thermo.dat'
      
      
