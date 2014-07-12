# Galileo E1-C code construction
#
# Copyright 2014 Peter Monta

import numpy as np

from e1c_strings import *

chip_rate = 1023000
code_length = 4092

secondary_code = np.array([0,0,1,1,1,0,0,0,0,0,0,0,1,0,1,0,1,1,0,1,1,0,0,1,0])
secondary_code = 1.0 - 2.0*secondary_code

def e1c_parse_hex(prn):
  s = e1c_strings[prn]
  n = code_length
  y = np.zeros(n)
  for i in range(n):
    nib = i/4
    bit = 3-(i%4)
    y[i] = (int(s[nib],16)>>bit)&1
  return y

codes = {}

def e1c_code(prn):
  if not codes.has_key(prn):
    codes[prn] = e1c_parse_hex(prn)
  return codes[prn]

def code(prn,chips,frac,incr,n):
  c = e1c_code(prn)
  idx = (chips%code_length) + frac + incr*np.arange(n)
  idx = np.floor(idx).astype('int')
  idx = np.mod(idx,code_length)
  x = c[idx]
  return 1.0 - 2.0*x

# test

if __name__=='__main__':
  print e1c_code(1)[0:20]
  print e1c_code(2)[0:20]
