# GPS L5Q code construction
#
# Copyright 2014 Peter Monta

import numpy as np

chip_rate = 10230000
code_length = 10230
nh_code = [0,0,0,0,0,1,0,0,1,1,0,1,0,1,0,0,1,1,1,0]

# initial-state table from pages 5--7 and pages 29--33 of IS-GPS-705D
# index is PRN

l5q_init = {
    1: 1701,     2:  323,     3: 5292,     4: 2020,
    5: 5429,     6: 7136,     7: 1041,     8: 5947,
    9: 4315,    10:  148,    11:  535,    12: 1939,
   13: 5206,    14: 5910,    15: 3595,    16: 5135,
   17: 6082,    18: 6990,    19: 3546,    20: 1523,
   21: 4548,    22: 4484,    23: 1893,    24: 3961,
   25: 7106,    26: 5299,    27: 4660,    28:  276,
   29: 4389,    30: 3783,    31: 1591,    32: 1601,
   33:  749,    34: 1387,    35: 1661,    36: 3210,
   37:  708,
   38: 4226,    39: 5604,    40: 6375,    41: 3056,
   42: 1772,    43: 3662,    44: 4401,    45: 5218,
   46: 2838,    47: 6913,    48: 1685,    49: 1194,
   50: 6963,    51: 5001,    52: 6694,    53:  991,
   54: 7489,    55: 2441,    56:  639,    57: 2097,
   58: 2498,    59: 6470,    60: 2399,    61:  242,
   62: 3768,    63: 1186,
   64: 5246,    65: 4259,    66: 5907,    67: 3870,
   68: 3262,    69: 7387,    70: 3069,    71: 2999,
   72: 7993,    73: 7849,    74: 4157,    75: 5031,
   76: 5986,    77: 4833,    78: 5739,    79: 7846,
   80:  898,    81: 2022,    82: 7446,    83: 6404,
   84:  155,    85: 7862,    86: 7795,    87: 6121,
   88: 4840,    89: 6585,    90:  429,    91: 6020,
   92:  200,    93: 1664,    94: 1499,    95: 7298,
   96: 1305,    97: 7323,    98: 7544,    99: 4438,
  100: 2485,   101: 3387,   102: 7319,   103: 1853,
  104: 5781,   105: 1874,   106: 7555,   107: 2132,
  108: 6441,   109: 6722,   110: 1192,   111: 2588,
  112: 2188,   113:  297,   114: 1540,   115: 4138,
  116: 5231,   117: 4789,   118:  659,   119:  871,
  120: 6837,   121: 1393,   122: 7383,   123:  611,
  124: 4920,   125: 5416,   126: 1611,   127: 2474,
  128:  118,   129: 1382,   130: 1092,   131: 7950,
  132: 7223,   133: 1769,   134: 4721,   135: 1252,
  136: 5147,   137: 2165,   138: 7897,   139: 4054,
  140: 3498,   141: 6571,   142: 2858,   143: 8126,
  144: 7017,   145: 1901,   146:  181,   147: 1114,
  148: 5195,   149: 7479,   150: 4186,   151: 3904,
  152: 7128,   153: 1396,   154: 4513,   155: 5967,
  156: 2580,   157: 2575,   158: 7961,   159: 2598,
  160: 4508,   161: 2090,   162: 3685,   163: 7748,
  164:  684,   165:  913,   166: 5558,   167: 2894,
  168: 5858,   169: 6432,   170: 3813,   171: 3573,
  172: 7523,   173: 5280,   174: 3376,   175: 7424,
  176: 2918,   177: 5793,   178: 1747,   179: 7079,
  180: 2921,   181: 2490,   182: 4119,   183: 3373,
  184:  977,   185:  681,   186: 4273,   187: 5419,
  188: 5626,   189: 1266,   190: 5804,   191: 2414,
  192: 6444,   193: 4757,   194:  427,   195: 5452,
  196: 5182,   197: 6606,   198: 6531,   199: 4268,
  200: 3115,   201: 6835,   202:  862,   203: 4856,
  204: 2765,   205:   37,   206: 1943,   207: 7977,
  208: 2512,   209: 4451,   210: 4071
}

def l5q_xa_shift(xa):
  if xa==[1,1,1,1,1,1,1,1,1,1,1,0,1]:
    return [1,1,1,1,1,1,1,1,1,1,1,1,1]
  else:
    return [xa[12]^xa[11]^xa[9]^xa[8]] + xa[0:12]

def l5q_xb_shift(xb):
  return [xb[12]^xb[11]^xb[7]^xb[6]^xb[5]^xb[3]^xb[2]^xb[0]] + xb[0:12]

def make_l5q_xa():
  xa = [1,1,1,1,1,1,1,1,1,1,1,1,1]
  y = np.zeros(code_length)
  for i in range(code_length):
    y[i] = xa[12]
    xa = l5q_xa_shift(xa)
  return y

def make_l5q_xb():
  xb = [1,1,1,1,1,1,1,1,1,1,1,1,1]
  y = np.zeros(code_length)
  for i in range(code_length):
    y[i] = xb[12]
    xb = l5q_xb_shift(xb)
  return y
  
xa = make_l5q_xa()
xb = make_l5q_xb()

codes = {}

def make_l5q(prn):
  xb_offset = l5q_init[prn]
  n = code_length
  xb_shift = xb[np.mod(np.arange(xb_offset,xb_offset+n),n)]
  return np.logical_xor(xa,xb_shift)

def l5q_code(prn):
  if prn not in codes:
    codes[prn] = make_l5q(prn)
  return codes[prn]

def code(prn,chips,frac,incr,n):
  c = l5q_code(prn)
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

# test vectors in IS-GPS-705D (for XB only unfortunately)

xb_start_state = {
  1: [1,0,0,1,0,1,1,0,0,1,1,0,0],
  2: [0,1,0,0,0,1,1,1,1,0,1,1,0]
}

def test_xb_start_state(prn):
  x = l5q_code(prn)
  y = []
  for i in range(13):
    bit = x[12-i].astype('int')
    y.append(1-bit)
  return y

if __name__=='__main__':
  for prn in list(xb_start_state.keys()):
    s = test_xb_start_state(prn)
    t = xb_start_state[prn]
    if s!=t:
      print("prn %d: ***mismatch*** %s %s" % (prn,s,t))
    else:
#      print("prn %d: %s %s" % (prn,s,t))
      pass
