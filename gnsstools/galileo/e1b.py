# Galileo E1-B code construction
#
# Copyright 2014 Peter Monta

import numpy as np

from e1b_strings import *

chip_rate = 1023000
code_length = 4092

def e1b_parse_hex(prn):
  s = e1b_strings[prn]
  n = code_length
  y = np.zeros(n)
  for i in range(n):
    nib = i/4
    bit = 3-(i%4)
    y[i] = (int(s[nib],16)>>bit)&1
  return y

codes = {}

def e1b_code(prn):
  if not codes.has_key(prn):
    codes[prn] = e1b_parse_hex(prn)
  return codes[prn]

def code(prn,chips,frac,incr,n):
  c = e1b_code(prn)
  idx = (chips%code_length) + frac + incr*np.arange(n)
  idx = np.floor(idx).astype('int')
  idx = np.mod(idx,code_length)
  x = c[idx]
  return 1.0 - 2.0*x

# test

if __name__=='__main__':
  print e1b_code(1)[0:20]
  print e1b_code(2)[0:20]
