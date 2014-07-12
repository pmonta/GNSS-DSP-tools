# GLONASS C/A code construction
#
# Copyright 2014 Peter Monta

import numpy as np

chip_rate = 511000
code_length = 511

def glonass_ca_shift(x):
  return [x[8]^x[4]] + x[0:8]

def make_glonass_ca():
  n = code_length
  x = [1,1,1,1,1,1,1,1,1]
  y = np.zeros(n)
  for i in range(n):
    y[i] = x[6]
    x = glonass_ca_shift(x)
  return y

c = make_glonass_ca()

def code(chips,frac,incr,n):
  idx = (chips%code_length) + frac + incr*np.arange(n)
  idx = np.floor(idx).astype('int')
  idx = np.mod(idx,code_length)
  x = c[idx]
  return 1.0 - 2.0*x

#
# testing: print out a small sample of the code
#

if __name__=='__main__':
  print c[0:20]
