# GPS L1Cp code construction
#
# Copyright 2018 Peter Monta

import numpy as np
from sympy.ntheory import legendre_symbol

chip_rate = 1023000
code_length = 10230

l1cp_params = {
   1: (5111,412),    2: (5109,161),     3: (5108,1),      4: (5106,303),
   5: (5103,207),    6: (5101,4971),    7: (5100,4496),   8: (5098,5),
   9: (5095,4557),  10: (5094,485),    11: (5093,253),   12: (5091,4676),
  13: (5090,1),     14: (5081,66),     15: (5080,4485),  16: (5069,282),
  17: (5068,193),   18: (5054,5211),   19: (5044,729),   20: (5027,4848),
  21: (5026,982),   22: (5014,5955),   23: (5004,9805),  24: (4980,670),
  25: (4915,464),   26: (4909,29),     27: (4893,429),   28: (4885,394),
  29: (4832,616),   30: (4824,9457),   31: (4591,4429),  32: (3706,4771),
  33: (5092,365),   34: (4986,9705),   35: (4965,9489),  36: (4920,4193),
  37: (4917,9947),  38: (4858,824),    39: (4847,864),   40: (4790,347),
  41: (4770,677),   42: (4318,6544),   43: (4126,6312),  44: (3961,9804),
  45: (3790,278),   46: (4911,9461),   47: (4881,444),   48: (4827,4839),
  49: (4795,4144),  50: (4789,9875),   51: (4725,197),   52: (4675,1156),
  53: (4539,4674),  54: (4535,10035),  55: (4458,4504),  56: (4197,5),
  57: (4096,9937),  58: (3484,430),    59: (3481,5),     60: (3393,355),
  61: (3175,909),   62: (2360,1622),   63: (1852,6284),
  193: (4311,9864),  194: (5024,9753),  195: (4352,9859),  199: (3646,164),
}

N = 10223
L = np.array([legendre_symbol(i,N) for i in range(N)])
L[L==-1] = 0
L[0] = 0

def l1cp(prn):
  w,p = l1cp_params[prn]
  W = np.array([L[k]^L[(k+w)%N] for k in range(N)])
  expansion = np.array([0,1,1,0,1,0,0])
  c = np.concatenate((W[0:p-1],expansion,W[p-1:N]))
  return c

codes = {}

def l1cp_code(prn):
  if prn not in codes:
    codes[prn] = l1cp(prn)
  return codes[prn]

def code(prn,chips,frac,incr,n):
  c = l1cp_code(prn)
  idx = (chips%code_length) + frac + incr*np.arange(n)
  idx = np.floor(idx).astype('int')
  idx = np.mod(idx,code_length)
  x = c[idx]
  return 1.0 - 2.0*x

l1cp_secondary_params = {
   1: (0o5111,0o3266),   2: (0o5421,0o2040),   3: (0o5501,0o1527),   4: (0o5403,0o3307),
   5: (0o6417,0o3756),   6: (0o6141,0o3026),   7: (0o6351,0o0562),   8: (0o6501,0o0420),
   9: (0o6205,0o3415),  10: (0o6235,0o0337),  11: (0o7751,0o0265),  12: (0o6623,0o1230),
  13: (0o6733,0o2204),  14: (0o7627,0o1440),  15: (0o5667,0o2412),  16: (0o5051,0o3516),
  17: (0o7665,0o2761),  18: (0o6325,0o3750),  19: (0o4365,0o2701),  20: (0o4745,0o1206),
  21: (0o7633,0o1544),  22: (0o6747,0o1774),  23: (0o4475,0o0546),  24: (0o4225,0o2213),
  25: (0o7063,0o3707),  26: (0o4423,0o2051),  27: (0o6651,0o3650),  28: (0o4161,0o1777),
  29: (0o7237,0o3203),  30: (0o4473,0o1762),  31: (0o5477,0o2100),  32: (0o6163,0o0571),
  33: (0o7223,0o3710),  34: (0o6323,0o3535),  35: (0o7125,0o3110),  36: (0o7035,0o1426),
  37: (0o4341,0o0255),  38: (0o4353,0o0321),  39: (0o4107,0o3124),  40: (0o5735,0o0572),
  41: (0o6741,0o1736),  42: (0o7071,0o3306),  43: (0o4563,0o1307),  44: (0o5755,0o3763),
  45: (0o6127,0o1604),  46: (0o4671,0o1021),  47: (0o4511,0o2624),  48: (0o4533,0o0406),
  49: (0o5357,0o0114),  50: (0o5607,0o0077),  51: (0o6673,0o3477),  52: (0o6153,0o1000),
  53: (0o7565,0o3460),  54: (0o7107,0o2607),  55: (0o6211,0o2057),  56: (0o4321,0o3467),
  57: (0o7201,0o0706),  58: (0o4451,0o2032),  59: (0o5411,0o1464),  60: (0o5141,0o0520),
  61: (0o7041,0o1766),  62: (0o6637,0o3270),  63: (0o4577,0o0341),
}

sec_code_length = 1800

def int2list(x,n):
  y = []
  for i in range(n):
    y.append((x>>i)&1)
  return y

def xorprod(a,b):
  t = 0
  for x,y in zip(a,b):
    t = t ^ (x*y)
  return t

def s1_shift(x,p):
  return [xorprod(x,p)] + x[0:-1]

def sec_l1cp(prn):
  p,init = l1cp_secondary_params[prn]
  p = int2list(p//2,11)
  x = int2list(init,11)
  c = np.zeros(sec_code_length)
  for i in range(sec_code_length):
    c[i] = x[10]
    x = s1_shift(x,p)
  return c

secondary_codes = {}

def secondary_code(prn):
  if prn not in secondary_codes:
    secondary_codes[prn] = sec_l1cp(prn)
  return secondary_codes[prn]

boc11 = np.array([1.0,-1.0])
tmboc_pattern = np.array([1,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0])

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
  bp6 = (12*(chips+frac))%2
  u = int(cp%33)
  for i in range(n):
    if tmboc_pattern[u]:
      boc = boc11[int(bp6)]
    else:
      boc = boc11[int(bp)]
    p += x[i]*(1.0-2.0*c[int(cp)])*boc
    cp = (cp+incr)%code_length
    bp = (bp+2*incr)%2
    bp6 = (bp6+12*incr)%2
    u = int(cp%33)
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
    c = l1cp_code(prn)
    s1 = chips2octal(c[0:24])
    s2 = chips2octal(c[-24:])
    print("%d %s %s"%(prn,s1,s2))
  print("secondary:")
  for prn in range(1,64):
    c = secondary_code(prn)
    print('%d %s'%(prn,chips2octal(np.concatenate((np.array([0]),c[-11:])))))
