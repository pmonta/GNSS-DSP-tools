import numpy as np

NT = 2048
nco_table = np.exp(2*(np.pi)*(1j)*np.arange(NT)*(1.0/NT))

def nco(f,p,n):
  idx = p + f*np.arange(n)
  idx = np.floor(idx*NT).astype('int')
  idx = np.mod(idx,NT)
  return nco_table[idx]
