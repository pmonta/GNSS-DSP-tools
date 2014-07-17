import numpy as np

NT = 1024
nco_table = np.exp(2*(np.pi)*(1j)*np.arange(NT)*(1.0/NT))

def nco(f,p,n):
  idx = p + f*np.arange(n)
  idx = np.floor(idx*NT).astype('int')
  idx = np.mod(idx,NT)
  return nco_table[idx]

def boc11(chips,frac,incr,n):
  c = np.array([-1,1])
  boc11_length = 2
  idx = (chips%boc11_length) + frac + incr*np.arange(n)
  idx = idx*2
  idx = np.floor(idx).astype('int')
  idx = np.mod(idx,boc11_length)
  return c[idx]

from numba import jit

@jit(nopython=True)
def mix(x,f,p,tab):
  n = len(x)
  dp = int(p*NT*65536)
  df = int(f*NT*65536)
  for i in range(n):
    idx = dp>>16
    x[i] *= tab[idx&(NT-1)]
    dp += df
