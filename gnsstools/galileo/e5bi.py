# Galileo E5b-I code construction
#
# Copyright 2014 Peter Monta

import numpy as np

chip_rate = 10230000
code_length = 10230

# secondary code from Table 19 (page 18) of Galileo ICD SIS (2014)

secondary_code = np.array([1,1,1,0])   # table 19, CS4_1, hex string 'E'
secondary_code = 1.0 - 2.0*secondary_code

# initial-state table from Table 17 (page 15) of Galileo ICD SIS (2014)
# index is PRN

e5bi_init = {
   1: 0o07220,    2: 0o26047,    3: 0o00252,    4: 0o17166,
   5: 0o14161,    6: 0o02540,    7: 0o01537,    8: 0o26023,
   9: 0o01725,   10: 0o20637,   11: 0o02364,   12: 0o27731,
  13: 0o30640,   14: 0o34174,   15: 0o06464,   16: 0o07676,
  17: 0o32231,   18: 0o10353,   19: 0o00755,   20: 0o26077,
  21: 0o11644,   22: 0o11537,   23: 0o35115,   24: 0o20452,
  25: 0o34645,   26: 0o25664,   27: 0o21403,   28: 0o32253,
  29: 0o02337,   30: 0o30777,   31: 0o27122,   32: 0o22377,
  33: 0o36175,   34: 0o33075,   35: 0o33151,   36: 0o13134,
  37: 0o07433,   38: 0o10216,   39: 0o35466,   40: 0o02533,
  41: 0o05351,   42: 0o30121,   43: 0o14010,   44: 0o32576,
  45: 0o30326,   46: 0o37433,   47: 0o26022,   48: 0o35770,
  49: 0o06670,   50: 0o12017
}

def e5bi_reg1_shift(x):
  return [x[13]^x[12]^x[10]^x[3]] + x[0:13]

def e5bi_reg2_shift(x):
  return [x[13]^x[11]^x[8]^x[7]^x[4]^x[1]] + x[0:13]

def make_e5bi_reg1():
  x = [1,1,1,1,1,1,1,1,1,1,1,1,1,1]
  y = np.zeros(code_length)
  for i in range(code_length):
    y[i] = x[13]
    x = e5bi_reg1_shift(x)
  return y

def make_e5bi_reg2(start):
  x = start
  y = np.zeros(code_length)
  for i in range(code_length):
    y[i] = x[13]
    x = e5bi_reg2_shift(x)
  return y

r1 = make_e5bi_reg1()

def seq(a):
  s = []
  for i in range(14):
    s = s + [(a>>i)&1]
  return s

codes = {}

def make_e5bi(prn):
  start = seq(e5bi_init[prn])
  r2 = make_e5bi_reg2(start)
  return np.logical_xor(r1,r2)

def e5bi_code(prn):
  if prn not in codes:
    codes[prn] = make_e5bi(prn)
  return codes[prn]

def code(prn,chips,frac,incr,n):
  c = e5bi_code(prn)
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

# test vectors in Galileo ICD

if __name__=='__main__':
  for prn in [1,2,3,4]:
    x = e5bi_code(prn)
    print(x[0:24].astype('int'))
