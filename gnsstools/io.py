import numpy as np

def get_samples_complex(fp,n):
  s = np.fromfile(fp,'b',2*n)
  if len(s)!=2*n:
    return None
  s.shape = (n,2)
  x = np.empty(n,dtype='c8')
  x.real = s[:,0]
  x.imag = s[:,1]
  return x
