import numpy as np

def nco(f,p,n):
  return np.exp(2*(np.pi)*(1j)*(p+f*np.arange(n)))
