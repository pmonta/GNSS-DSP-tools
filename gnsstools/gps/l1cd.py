# GPS L1Cd code construction
#
# Copyright 2018 Peter Monta

import numpy as np
from sympy.ntheory import legendre_symbol

chip_rate = 1023000
code_length = 10230

l1cd_params = {
   1: (5097,181),     2: (5110,359),    3: (5079,72),      4: (4403,1110),
   5: (4121,1480),    6: (5043,5034),   7: (5042,4622),    8: (5104,1),
   9: (4940,4547),   10: (5035,826),   11: (4372,6284),   12: (5064,4195),
  13: (5084,368),    14: (5048,1),     15: (4950,4796),   16: (5019,523),
  17: (5076,151),    18: (3736,713),   19: (4993,9850),   20: (5060,5734),
  21: (5061,34),     22: (5096,6142),  23: (4983,190),    24: (4783,644),
  25: (4991,467),    26: (4815,5384),  27: (4443,801),    28: (4769,594),
  29: (4879,4450),   30: (4894,9437),  31: (4985,4307),   32: (5056,5906),
  33: (4921,378),    34: (5036,9448),  35: (4812,9432),   36: (4838,5849),
  37: (4855,5547),   38: (4904,9546),  39: (4753,9132),   40: (4483,403),
  41: (4942,3766),   42: (4813,3),     43: (4957,684),    44: (4618,9711),
  45: (4669,333),    46: (4969,6124),  47: (5031,10216),  48: (5038,4251),
  49: (4740,9893),   50: (4073,9884),  51: (4843,4627),   52: (4979,4449),
  53: (4867,9798),   54: (4964,985),   55: (5025,4272),   56: (4579,126),
  57: (4390,10024),  58: (4763,434),   59: (4612,1029),   60: (4784,561),
  61: (3716,289),    62: (4703,638),   63: (4851,4353),
}

N = 10223
L = np.array([legendre_symbol(i,N) for i in range(N)])
L[L==-1] = 0
L[0] = 0

def l1cd(prn):
  w,p = l1cd_params[prn]
  W = np.array([L[k]^L[(k+w)%N] for k in range(N)])
  expansion = np.array([0,1,1,0,1,0,0])
  c = np.concatenate((W[0:p-1],expansion,W[p-1:N]))
  return c

codes = {}

def l1cd_code(prn):
  if prn not in codes:
    codes[prn] = l1cd(prn)
  return codes[prn]

def code(prn,chips,frac,incr,n):
  c = l1cd_code(prn)
  idx = (chips%code_length) + frac + incr*np.arange(n)
  idx = np.floor(idx).astype('int')
  idx = np.mod(idx,code_length)
  x = c[idx]
  return 1.0 - 2.0*x

boc11 = np.array([1.0,-1.0])

try:
  from numba import jit
except:
  def jit(**kwargs):
    return lambda x: x

@jit(nopython=True)
def correlate(x,prn,chips,frac,incr,c,boc11):
  n = len(x)
  p = 0.0j
  cp = (chips+frac)%code_length
  bp = (2*(chips+frac))%2
  for i in range(n):
    boc = boc11[int(bp)]
    p += x[i]*(1.0-2.0*c[int(cp)])*boc
    cp = (cp+incr)%code_length
    bp = (bp+2*incr)%2
  return p

# test

def chips2octal(x):
  s = ''
  for i in range(len(x)//3):
    d = 4*x[3*i] + 2*x[3*i+1] + x[3*i+2]
    s = s + '%o'%int(d)
  return s

if __name__=='__main__':
  for prn in range(1,64):
    c = l1cd_code(prn)
    s1 = chips2octal(c[0:24])
    s2 = chips2octal(c[-24:])
    print("%d %s %s"%(prn,s1,s2))
