import numpy as np

def get_samples_complex(fp,n):
  z = fp.read(2*n)
  if len(z)!=2*n:
    return None
  s = np.fromstring(z,dtype='int8')
  s.shape = (n,2)
  x = np.empty(n,dtype='c8')
  x.real = s[:,0]
  x.imag = s[:,1]
  return x
