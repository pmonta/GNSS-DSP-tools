# GLONASS P code construction
#
# Copyright 2014 Peter Monta

import numpy as np

chip_rate = 5110000
code_length = 5110000

def glonass_p_shift(x):
  return [x[24]^x[2]] + x[0:24]

def make_glonass_p():
  n = code_length
  x = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]
  y = np.zeros(n)
  for i in range(n):
    y[i] = x[9]
    x = glonass_p_shift(x)
  return y

c = make_glonass_p()

def p_code():
  return c

def code(chips,frac,incr,n):
  idx = (chips%code_length) + frac + incr*np.arange(n)
  idx = np.floor(idx).astype('int')
  idx = np.mod(idx,code_length)
  x = c[idx]
  return 1.0 - 2.0*x

try:
  from numba import jit
except:
  def jit(**kwargs):
    return lambda x: x

@jit(nopython=True)
def correlate(x,chips,frac,incr,c):
  n = len(x)
  p = 0.0j
  cp = (chips+frac)%code_length
  for i in range(n):
    p += x[i]*(1.0-2.0*c[int(cp)])
    cp += incr
    if cp>=code_length:
      cp -= code_length
  return p

#
# testing: print out a small sample of the code
#

if __name__=='__main__':
  print(c[0:100])
