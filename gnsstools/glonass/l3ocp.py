# GLONASS L3OCp (pilot signal) code construction
#
# Copyright 2014 Peter Monta

import numpy as np

chip_rate = 10230000
code_length = 10230

secondary_code = np.array([0,0,0,0,1,1,0,1,0,1])
secondary_code = 1.0 - 2.0*secondary_code

def g2_shift(x):
  return [x[13]^x[12]^x[7]^x[3]] + x[0:13]

def g1_shift(x):
  return [x[6]^x[5]] + x[0:6]

def seq(n):
  s = []
  for i in range(7):
    s = s + [(n>>(6-i))&1]
  return s

def make_l3ocp(n):
  g1 = seq(n+64)
  g2 = [0,0,1,1,0,1,0,0,1,1,1,0,0,0]
  x = np.zeros(code_length)
  for i in range(code_length):
    x[i] = g1[6]^g2[13]
    g1 = g1_shift(g1)
    g2 = g2_shift(g2)
  return x

codes = {}

def l3ocp_code(n):
  if n not in codes:
    codes[n] = make_l3ocp(n)
  return codes[n]

def code(prn,chips,frac,incr,n):
  c = l3ocp_code(prn)
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
def correlate(x,prn,chips,frac,incr,c):
  n = len(x)
  p = 0.0j
  cp = (chips+frac)%code_length
  for i in range(n):
    p += x[i]*(1.0-2.0*c[int(cp)])
    cp = (cp+incr)%code_length
  return p

# print some codes (no ICD document or test vectors available, apparently)

if __name__=='__main__':
  import sys
  c = l3ocp_code(30)
  for i in range(200):
    sys.stdout.write('%d'%c[i])
  sys.stdout.write('\n')
