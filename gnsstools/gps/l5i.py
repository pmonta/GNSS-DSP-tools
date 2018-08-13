# GPS L5I code construction
#
# Copyright 2014 Peter Monta

import numpy as np

chip_rate = 10230000
code_length = 10230

secondary_code = np.array([0,0,0,0,1,1,0,1,0,1])
secondary_code = 1.0 - 2.0*secondary_code

# initial-state table from pages 5--7 and pages 29--33 of IS-GPS-705D
# index is PRN

l5i_init = {
    1:  266,     2:  365,     3:  804,     4: 1138,
    5: 1509,     6: 1559,     7: 1756,     8: 2084,
    9: 2170,    10: 2303,    11: 2527,    12: 2687,
   13: 2930,    14: 3471,    15: 3940,    16: 4132,
   17: 4332,    18: 4924,    19: 5343,    20: 5443,
   21: 5641,    22: 5816,    23: 5898,    24: 5918,
   25: 5955,    26: 6243,    27: 6345,    28: 6477,
   29: 6518,    30: 6875,    31: 7168,    32: 7187,
   33: 7329,    34: 7577,    35: 7720,    36: 7777,
   37: 8057,
   38: 5358,    39: 3550,    40: 3412,    41:  819,
   42: 4608,    43: 3698,    44:  962,    45: 3001,
   46: 4441,    47: 4937,    48: 3717,    49: 4730,
   50: 7291,    51: 2279,    52: 7613,    53: 5723,
   54: 7030,    55: 1475,    56: 2593,    57: 2904,
   58: 2056,    59: 2757,    60: 3756,    61: 6205,
   62: 5053,    63: 6437,
   64: 7789,    65: 2311,    66: 7432,    67: 5155,
   68: 1593,    69: 5841,    70: 5014,    71: 1545,
   72: 3016,    73: 4875,    74: 2119,    75:  229,
   76: 7634,    77: 1406,    78: 4506,    79: 1819,
   80: 7580,    81: 5446,    82: 6053,    83: 7958,
   84: 5267,    85: 2956,    86: 3544,    87: 1277,
   88: 2996,    89: 1758,    90: 3360,    91: 2718,
   92: 3754,    93: 7440,    94: 2781,    95: 6756,
   96: 7314,    97:  208,    98: 5252,    99:  696,
  100:  527,   101: 1399,   102: 5879,   103: 6868,
  104:  217,   105: 7681,   106: 3788,   107: 1337,
  108: 2424,   109: 4243,   110: 5686,   111: 1955,
  112: 4791,   113:  492,   114: 1518,   115: 6566,
  116: 5349,   117:  506,   118:  113,   119: 1953,
  120: 2797,   121:  934,   122: 3023,   123: 3632,
  124: 1330,   125: 4909,   126: 4867,   127: 1183,
  128: 3990,   129: 6217,   130: 1224,   131: 1733,
  132: 2319,   133: 3928,   134: 2380,   135:  841,
  136: 5049,   137: 7027,   138: 1197,   139: 7208,
  140: 8000,   141:  152,   142: 6762,   143: 3745,
  144: 4723,   145: 5502,   146: 4796,   147:  123,
  148: 8142,   149: 5091,   150: 7875,   151:  330,
  152: 5272,   153: 4912,   154:  374,   155: 2045,
  156: 6616,   157: 6321,   158: 7605,   159: 2570,
  160: 2419,   161: 1234,   162: 1922,   163: 4317,
  164: 5110,   165:  825,   166:  958,   167: 1089,
  168: 7813,   169: 6058,   170: 7703,   171: 6702,
  172: 1714,   173: 6371,   174: 2281,   175: 1986,
  176: 6282,   177: 3201,   178: 3760,   179: 1056,
  180: 6233,   181: 1150,   182: 2823,   183: 6250,
  184:  645,   185: 2401,   186: 1639,   187: 2946,
  188: 7091,   189:  923,   190: 7045,   191: 6493,
  192: 1706,   193: 5836,   194:  926,   195: 6086,
  196:  950,   197: 5905,   198: 3240,   199: 6675,
  200: 3197,   201: 1555,   202: 3589,   203: 4555,
  204: 5671,   205: 6948,   206: 4664,   207: 2086,
  208: 5950,   209: 5521,   210: 1515
}

def l5i_xa_shift(xa):
  if xa==[1,1,1,1,1,1,1,1,1,1,1,0,1]:
    return [1,1,1,1,1,1,1,1,1,1,1,1,1]
  else:
    return [xa[12]^xa[11]^xa[9]^xa[8]] + xa[0:12]

def l5i_xb_shift(xb):
  return [xb[12]^xb[11]^xb[7]^xb[6]^xb[5]^xb[3]^xb[2]^xb[0]] + xb[0:12]

def make_l5i_xa():
  xa = [1,1,1,1,1,1,1,1,1,1,1,1,1]
  y = np.zeros(code_length)
  for i in range(code_length):
    y[i] = xa[12]
    xa = l5i_xa_shift(xa)
  return y

def make_l5i_xb():
  xb = [1,1,1,1,1,1,1,1,1,1,1,1,1]
  y = np.zeros(code_length)
  for i in range(code_length):
    y[i] = xb[12]
    xb = l5i_xb_shift(xb)
  return y
  
xa = make_l5i_xa()
xb = make_l5i_xb()

codes = {}

def make_l5i(prn):
  xb_offset = l5i_init[prn]
  n = code_length
  xb_shift = xb[np.mod(np.arange(xb_offset,xb_offset+n),n)]
  return np.logical_xor(xa,xb_shift)

def l5i_code(prn):
  if prn not in codes:
    codes[prn] = make_l5i(prn)
  return codes[prn]

def code(prn,chips,frac,incr,n):
  c = l5i_code(prn)
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
  1: [0,1,0,1,0,1,1,1,0,0,1,0,0],
  2: [1,1,0,0,0,0,0,1,1,0,1,0,1]
}

def test_xb_start_state(prn):
  x = l5i_code(prn)
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
