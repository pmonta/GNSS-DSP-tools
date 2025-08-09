# Xona X1 pilot code
#
# Copyright 2025 Peter Monta

import numpy as np

from .x1p_strings import *

chip_rate = 1023000
code_length = 1023

secondary_code = np.array([0,1,0,0,0,1,0,0,1,1,1,1,0,1,1,1,0,0,1,0,1,1,1,0,0,1,1,0,0,1,0,0,1,0,1,1,1,1,1,1,0,1,1,1,0,0,0,1,1,0,1,0,1,1,1,0,1,0,0,0,1,1,0,0,1,1,1,1,1,1,1,0,1,1,1,1,1,0,0,1,0,1,0,0,1,0,0,0,0,1,1,1,1,0,0,0,0,0,0,1])
secondary_code = 1.0 - 2.0*secondary_code

def x1p_parse_hex(prn):
  s = x1p_strings[prn]
  n = code_length
  y = np.zeros(n)
  for i in range(n):
    nib = i//4
    bit = 3-(i%4)
    y[i] = (int(s[nib],16)>>bit)&1
  return y

codes = {}

def x1p_code(prn):
  if prn not in codes:
    codes[prn] = x1p_parse_hex(prn)
  return codes[prn]

def code(prn,chips,frac,incr,n):
  c = x1p_code(prn)
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

if __name__=='__main__':
  print(x1p_code('pulsar-0')[0:20])
