# GPS C/A code construction
#
# Copyright 2012-2014 Peter Monta

import numpy as np

chip_rate = 1023000
code_length = 1023

# G2 delay table from pages 6--8 and pages 57--61 of IS-GPS-200H
# index is PRN, value is G2 delay

g2_delay = {
    1:    5,    2:    6,    3:    7,    4:    8,
    5:   17,    6:   18,    7:  139,    8:  140,
    9:  141,   10:  251,   11:  252,   12:  254,
   13:  255,   14:  256,   15:  257,   16:  258,
   17:  469,   18:  470,   19:  471,   20:  472,
   21:  473,   22:  474,   23:  509,   24:  512,
   25:  513,   26:  514,   27:  515,   28:  516,
   29:  859,   30:  860,   31:  861,   32:  862,
   33:  863,   34:  950,   35:  947,   36:  948,
   37:  950,
   38:   67,   39:  103,   40:   91,   41:   19,
   42:  679,   43:  225,   44:  625,   45:  946,
   46:  638,   47:  161,   48: 1001,   49:  554,
   50:  280,   51:  710,   52:  709,   53:  775,
   54:  864,   55:  558,   56:  220,   57:  397,
   58:   55,   59:  898,   60:  759,   61:  367,
   62:  299,   63: 1018,
   64:  729,   65:  695,   66:  780,   67:  801,
   68:  788,   69:  732,   70:   34,   71:  320,
   72:  327,   73:  389,   74:  407,   75:  525,
   76:  405,   77:  221,   78:  761,   79:  260,
   80:  326,   81:  955,   82:  653,   83:  699,
   84:  422,   85:  188,   86:  438,   87:  959,
   88:  539,   89:  879,   90:  677,   91:  586,
   92:  153,   93:  792,   94:  814,   95:  446,
   96:  264,   97: 1015,   98:  278,   99:  536,
  100:  819,  101:  156,  102:  957,  103:  159,
  104:  712,  105:  885,  106:  461,  107:  248,
  108:  713,  109:  126,  110:  807,  111:  279,
  112:  122,  113:  197,  114:  693,  115:  632,
  116:  771,  117:  467,  118:  647,  119:  203,
  120:  145,  121:  175,  122:   52,  123:   21,
  124:  237,  125:  235,  126:  886,  127:  657,
  128:  634,  129:  762,  130:  355,  131: 1012,
  132:  176,  133:  603,  134:  130,  135:  359,
  136:  595,  137:   68,  138:  386,  139:  797,
  140:  456,  141:  499,  142:  883,  143:  307,
  144:  127,  145:  211,  146:  121,  147:  118,
  148:  163,  149:  628,  150:  853,  151:  484,
  152:  289,  153:  811,  154:  202,  155: 1021,
  156:  463,  157:  568,  158:  904,  159:  670,
  160:  230,  161:  911,  162:  684,  163:  309,
  164:  644,  165:  932,  166:   12,  167:  314,
  168:  891,  169:  212,  170:  185,  171:  675,
  172:  503,  173:  150,  174:  395,  175:  345,
  176:  846,  177:  798,  178:  992,  179:  357,
  180:  995,  181:  877,  182:  112,  183:  144,
  184:  476,  185:  193,  186:  109,  187:  445,
  188:  291,  189:   87,  190:  399,  191:  292,
  192:  901,  193:  339,  194:  208,  195:  711,
  196:  189,  197:  263,  198:  537,  199:  663,
  200:  942,  201:  173,  202:  900,  203:   30,
  204:  500,  205:  935,  206:  556,  207:  373,
  208:   85,  209:  652,  210:  310,
}

def g1_shift(x):
  return [x[9]^x[2]] + x[0:9]

def g2_shift(x):
  return [x[9]^x[8]^x[7]^x[5]^x[2]^x[1]] + x[0:9]

def make_g1():
  n = code_length
  x = [1,1,1,1,1,1,1,1,1,1]
  g1 = np.zeros(code_length)
  for i in range(n):
    g1[i] = x[9]
    x = g1_shift(x)
  return g1

def make_g2():
  n = code_length
  x = [1,1,1,1,1,1,1,1,1,1]
  g2 = np.zeros(code_length)
  for i in range(n):
    g2[i] = x[9]
    x = g2_shift(x)
  return g2

def circular_shift(g,n):
  return np.concatenate([g[code_length-n:code_length],g[0:code_length-n]])

g1 = make_g1()
g2 = make_g2()
codes = {}

def ca_code(prn):
  if prn not in codes:
    codes[prn] = np.logical_xor(g1,circular_shift(g2,g2_delay[prn]))
  return codes[prn]

def code(prn,chips,frac,incr,n):
  c = ca_code(prn)
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

def correlate_slow(x,prn,chips,frac,incr,c):
  n = len(x)
  q = code(prn,chips,frac,incr,n)
  return np.sum(x*q)

# test vectors in IS-GPS-200H

def first_10_chips(prn):
  c = ca_code(prn)
  r = 0
  for i in range(0,10):
    r = 2*r + c[i]
  return r

def print_first_10_chips():
  for prn in list(g2_delay.keys()):
    print('%d: %04o' % (prn,first_10_chips(prn)))

if __name__=='__main__':
  print_first_10_chips()
