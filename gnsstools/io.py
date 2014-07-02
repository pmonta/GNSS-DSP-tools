import numpy as np

def get_samples_complex(fp,n):
  s = np.fromfile(fp,'b',2*n)
  if len(s)!=2*n:
    return None
  x = s[0:2*n:2] + (1j)*s[1:2*n:2]
  return x
