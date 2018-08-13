# Beidou B1Cd code construction
#
# Copyright 2018 Peter Monta

import numpy as np
from sympy.ntheory import legendre_symbol

chip_rate = 1023000
code_length = 10230

b1cd_params = {
   1: (2678,699),    2: (4802,694),    3: (958,7318),    4: (859,2127),
   5: (3843,715),    6: (2232,6682),   7: (124,7850),    8: (4352,5495),
   9: (1816,1162),  10: (1126,7682),  11: (1860,6792),  12: (4800,9973),
  13: (2267,6596),  14: (424,2092),   15: (4192,19),    16: (4333,10151),
  17: (2656,6297),  18: (4148,5766),  19: (243,2359),   20: (1330,7136),
  21: (1593,1706),  22: (1470,2128),  23: (882,6827),   24: (3202,693),
  25: (5095,9729),  26: (2546,1620),  27: (1733,6805),  28: (4795,534),
  29: (4577,712),   30: (1627,1929),  31: (3638,5355),  32: (2553,6139),
  33: (3646,6339),  34: (1087,1470),  35: (1843,6867),  36: (216,7851),
  37: (2245,1162),  38: (726,7659),   39: (1966,1156),  40: (670,2672),
  41: (4130,6043),  42: (53,2862),    43: (4830,180),   44: (182,2663),
  45: (2181,6940),  46: (2006,1645),  47: (1080,1582),  48: (2288,951),
  49: (2027,6878),  50: (271,7701),   51: (915,1823),   52: (497,2391),
  53: (139,2606),   54: (3693,822),   55: (2054,6403),  56: (4342,239),
  57: (3342,442),   58: (2592,6769),  59: (1007,2560),  60: (310,2502),
  61: (4203,5072),  62: (455,7268),   63: (4318,341),
}

N = 10243
L = np.array([legendre_symbol(i,N) for i in range(N)])
L[L==-1] = 0
L[0] = 0

def b1cd(prn):
  w,p = b1cd_params[prn]
  W = np.array([L[k]^L[(k+w)%N] for k in range(N)])
  c = np.array([W[(n+p-1)%N] for n in range(code_length)])
  return c

codes = {}

def b1cd_code(prn):
  if prn not in codes:
    codes[prn] = b1cd(prn)
  return codes[prn]

def code(prn,chips,frac,incr,n):
  c = b1cd_code(prn)
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
    s = s + '%o'%d
  return s

if __name__=='__main__':
  for prn in range(1,64):
    c = b1cd_code(prn)
    s1 = chips2octal(c[0:24])
    s2 = chips2octal(c[-24:])
    print("%d %s %s"%(prn,s1,s2))
