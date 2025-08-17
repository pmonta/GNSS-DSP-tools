# Xona X5 data code
#
# Copyright 2025 Peter Monta

import numpy as np

from .x5d_strings import *

chip_rate = 10230000
code_length = 10230

def x5d_parse_hex(prn):
  s = x5d_strings[prn]
  n = code_length
  y = np.zeros(n)
  for i in range(n):
    nib = i//4
    bit = 3-(i%4)
    y[i] = (int(s[nib],16)>>bit)&1
  return y

codes = {}

def x5d_code(prn):
  if prn not in codes:
    codes[prn] = x5d_parse_hex(prn)
  return codes[prn]

def code(prn,chips,frac,incr,n):
  c = x5d_code(prn)
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
  print(x5d_code(0)[0:20])
