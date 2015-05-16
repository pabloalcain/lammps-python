"""
Tally all the results given from analysis of the system
"""

import neutronstar.analysis as analysis

def rdf(computes, collective, rdf, n):
    computes['rdf'] *= (n-1)/n
    computes['rdf'] += rdf/n
    
def ssf(computes, collective, (ssf, k, s), n):
    computes['ssf'] *= (n-1)/n
    computes['ssf'] += ssf/n
    collective['k_absorption'].append(k)
    collective['S_absorption'].append(s)

def fit(computes, collective, (hei, dh, lam, dl), n):
    collective['height'].append(hei)
    collective['del_height'].append(dh)
    collective['lambda'].append(lam)
    collective['del_lambda'].append(dl)

def mste(computes, collective, (occ, frac, avg, std), n):
    new_occ = (n-1) * computes['mste'][0] + occ
    # If the cluster never occurs, its proton fraction is meaningless
    new_occ[new_occ == 0] = 1
    
    # Normalization of the proton fraction
    computes['mste'][1] *= (n-1) * computes['mste'][0]
    computes['mste'][1] += frac * occ
    computes['mste'][1] /= new_occ
    
    # Normalization of the occupancy
    computes['mste'][0] *= (n-1)/n
    computes['mste'][0] += occ/n
    
    collective['size_avg'].append(avg)
    collective['size_std'].append(std)
    
def mink(computes, collective, (vol, sur, bre, eul), n):
    collective['volume'].append(vol)
    collective['surface'].append(sur)
    collective['breadth'].append(bre)
    collective['euler'].append(eul)
    
def thermo(computes, collective, (tem, kin, pot, ene, pre), n):
    collective['temperature'].append(tem)
    collective['kinetic'].append(kin)
    collective['potential'].append(pot)
    collective['energy'].append(ene)
    collective['pressure'].append(pre)
    
def lind(computes, collective, lind, n):
    pass
