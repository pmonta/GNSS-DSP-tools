import numpy as np

try:
  from numba import jit
except:
  def jit(**kwargs):
    return lambda x: x

#
# decimate by n (boxcar filter), square, do m incoherent sums
#

@jit(nopython=True)
def squaring(x,r,n,m):
  blocks = len(x)//(n*m)
  for b in range(blocks):
    r[b] = 0j
    q = b*n*m
    for k in range(m):
      s = 0j
      for l in range(n):
        s += x[q+k*n+l]
      r[b] += s*s/n
