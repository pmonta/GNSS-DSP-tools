# Beidou B1Cp code construction
#
# Copyright 2018 Peter Monta

import numpy as np
from sympy.ntheory import legendre_symbol

chip_rate = 1023000
code_length = 10230

b1cp_params = {
   1: (796,7575),     2: (156,2369),     3: (4198,5688),    4: (3941,539),
   5: (1374,2270),    6: (1338,7306),    7: (1833,6457),    8: (2521,6254),
   9: (3175,5644),   10: (168,7119),    11: (2715,1402),   12: (4408,5557),
  13: (3160,5764),   14: (2796,1073),   15: (459,7001),    16: (3594,5910),
  17: (4813,10060),  18: (586,2710),    19: (1428,1546),   20: (2371,6887),
  21: (2285,1883),   22: (3377,5613),   23: (4965,5062),   24: (3779,1038),
  25: (4547,10170),  26: (1646,6484),   27: (1430,1718),   28: (607,2535),
  29: (2118,1158),   30: (4709,526),    31: (1149,7331),   32: (3283,5844),
  33: (2473,6423),   34: (1006,6968),   35: (3670,1280),   36: (1817,1838),
  37: (771,1989),    38: (2173,6468),   39: (740,2091),    40: (1433,1581),
  41: (2458,1453),   42: (3459,6252),   43: (2155,7122),   44: (1205,7711),
  45: (413,7216),    46: (874,2113),    47: (2463,1095),   48: (1106,1628),
  49: (1590,1713),   50: (3873,6102),   51: (4026,6123),   52: (4272,6070),
  53: (3556,1115),   54: (128,8047),    55: (1200,6795),   56: (130,2575),
  57: (4494,53),     58: (1871,1729),   59: (3073,6388),   60: (4386,682),
  61: (4098,5565),   62: (1923,7160),   63: (1176,2277),
}

N = 10243
L = np.array([legendre_symbol(i,N) for i in range(N)])
L[L==-1] = 0
L[0] = 0

def b1cp(prn):
  w,p = b1cp_params[prn]
  W = np.array([L[k]^L[(k+w)%N] for k in range(N)])
  c = np.array([W[(n+p-1)%N] for n in range(code_length)])
  return c

codes = {}

def b1cp_code(prn):
  if prn not in codes:
    codes[prn] = b1cp(prn)
  return codes[prn]

def code(prn,chips,frac,incr,n):
  c = b1cp_code(prn)
  idx = (chips%code_length) + frac + incr*np.arange(n)
  idx = np.floor(idx).astype('int')
  idx = np.mod(idx,code_length)
  x = c[idx]
  return 1.0 - 2.0*x

b1cp_secondary_params = {
   1: (269,1889),    2: (1448,1268),   3: (1028,1593),   4: (1324,1186),
   5: (822,1239),    6: (5,1930),      7: (155,176),     8: (458,1696),
   9: (310,26),     10: (959,1344),   11: (1238,1271),  12: (1180,1182),
  13: (1288,1381),  14: (334,1604),   15: (885,1333),   16: (1362,1185),
  17: (181,31),     18: (1648,704),   19: (838,1190),   20: (313,1646),
  21: (750,1385),   22: (225,113),    23: (1477,860),   24: (309,1656),
  25: (108,1921),   26: (1457,1173),  27: (149,1928),   28: (322,57),
  29: (271,150),    30: (576,1214),   31: (1103,1148),  32: (450,1458),
  33: (399,1519),   34: (241,1635),   35: (1045,1257),  36: (164,1687),
  37: (513,1382),   38: (687,1514),   39: (422,1),      40: (303,1583),
  41: (324,1806),   42: (495,1664),   43: (725,1338),   44: (780,1111),
  45: (367,1706),   46: (882,1543),   47: (631,1813),   48: (37,228),
  49: (647,2871),   50: (1043,2884),  51: (24,1823),    52: (120,75),
  53: (134,11),     54: (136,63),     55: (158,1937),   56: (214,22),
  57: (335,1768),   58: (340,1526),   59: (661,1402),   60: (889,1445),
  61: (929,1680),   62: (1002,1290),  63: (1149,1245),
}

sec_N = 3607
sec_L = np.array([legendre_symbol(i,sec_N) for i in range(sec_N)])
sec_L[sec_L==-1] = 0
sec_L[0] = 0

sec_code_length = 1800

def sec_b1cp(prn):
  w,p = b1cp_secondary_params[prn]
  W = np.array([sec_L[k]^sec_L[(k+w)%sec_N] for k in range(sec_N)])
  c = np.array([W[(n+p-1)%sec_N] for n in range(sec_code_length)])
  return c

secondary_codes = {}

def secondary_code(prn):
  if prn not in secondary_codes:
    secondary_codes[prn] = sec_b1cp(prn)
  return secondary_codes[prn]

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
    c = b1cp_code(prn)
    s1 = chips2octal(c[0:24])
    s2 = chips2octal(c[-24:])
    print("%d %s %s"%(prn,s1,s2))
  print("secondary:")
  for prn in range(1,64):
    c = secondary_code(prn)
    s1 = chips2octal(c[0:24])
    s2 = chips2octal(c[-24:])
    print("%d %s %s"%(prn,s1,s2))
