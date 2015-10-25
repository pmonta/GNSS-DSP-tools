# Galileo E5a-I code construction
#
# Copyright 2014 Peter Monta

import numpy as np

chip_rate = 10230000
code_length = 10230

# secondary code from Table 19 (page 18) of Galileo ICD SIS (2014)

secondary_code = '842E9'

# parse the hex string, return a 20-bit secondary-code array

def secondary_seq(s):
  x = np.zeros(20)
  for i in range(20):
    nib = i//4
    bit = 3-(i%4)
    x[i] = (int(s[nib],16)>>bit)&1
    x[i] = 1.0 - 2.0*x[i]
  return x

# transform the secondary-code table entries to sequences

secondary_code = secondary_seq(secondary_code)

# initial-state table from Table 15 (page 15) of Galileo ICD SIS (2014)
# index is PRN

e5ai_init = {
   1: 0o30305,    2: 0o14234,    3: 0o27213,    4: 0o20577,
   5: 0o23312,    6: 0o33463,    7: 0o15614,    8: 0o12537,    
   9: 0o01527,   10: 0o30236,   11: 0o27344,   12: 0o07272,    
  13: 0o36377,   14: 0o17046,   15: 0o06434,   16: 0o15405,    
  17: 0o24252,   18: 0o11631,   19: 0o24776,   20: 0o00630,    
  21: 0o11560,   22: 0o17272,   23: 0o27445,   24: 0o31702,    
  25: 0o13012,   26: 0o14401,   27: 0o34727,   28: 0o22627,    
  29: 0o30623,   30: 0o27256,   31: 0o01520,   32: 0o14211,    
  33: 0o31465,   34: 0o22164,   35: 0o33516,   36: 0o02737,    
  37: 0o21316,   38: 0o35425,   39: 0o35633,   40: 0o24655,    
  41: 0o14054,   42: 0o27027,   43: 0o06604,   44: 0o31455,    
  45: 0o34465,   46: 0o25273,   47: 0o20763,   48: 0o31721,    
  49: 0o17312,   50: 0o13277
}

def e5ai_reg1_shift(x):
  return [x[13]^x[7]^x[5]^x[0]] + x[0:13]

def e5ai_reg2_shift(x):
  return [x[13]^x[11]^x[7]^x[6]^x[4]^x[3]] + x[0:13]

def make_e5ai_reg1():
  x = [1,1,1,1,1,1,1,1,1,1,1,1,1,1]
  y = np.zeros(code_length)
  for i in range(code_length):
    y[i] = x[13]
    x = e5ai_reg1_shift(x)
  return y

def make_e5ai_reg2(start):
  x = start
  y = np.zeros(code_length)
  for i in range(code_length):
    y[i] = x[13]
    x = e5ai_reg2_shift(x)
  return y

r1 = make_e5ai_reg1()

def seq(a):
  s = []
  for i in range(14):
    s = s + [(a>>i)&1]
  return s

codes = {}

def make_e5ai(prn):
  start = seq(e5ai_init[prn])
  r2 = make_e5ai_reg2(start)
  return np.logical_xor(r1,r2)

def e5ai_code(prn):
  if prn not in codes:
    codes[prn] = make_e5ai(prn)
  return codes[prn]

def code(prn,chips,frac,incr,n):
  c = e5ai_code(prn)
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
    x = e5ai_code(prn)
    print(x[0:24].astype('int'))
