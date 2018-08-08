# Beidou B2ad code construction
#
# Copyright 2018 Peter Monta

import numpy as np

chip_rate = 10230000
code_length = 10230

secondary_code = np.array([0,0,0,1,0])
secondary_code = 1.0 - 2.0*secondary_code

b2ad_g2_initial = {
   1: "1000000100101",   2: "1000000110100",   3: "1000010101101",   4: "1000101001111",
   5: "1000101010101",   6: "1000110101110",   7: "1000111101110",   8: "1000111111011",
   9: "1001100101001",  10: "1001111011010",  11: "1010000110101",  12: "1010001000100",
  13: "1010001010101",  14: "1010001011011",  15: "1010001011100",  16: "1010010100011",
  17: "1010011110111",  18: "1010100000001",  19: "1010100111110",  20: "1010110101011",
  21: "1010110110001",  22: "1011001010011",  23: "1011001100010",  24: "1011010011000",
  25: "1011010110110",  26: "1011011110010",  27: "1011011111111",  28: "1011100010010",
  29: "1011100111100",  30: "1011110100001",  31: "1011111001000",  32: "1011111010100",
  33: "1011111101011",  34: "1011111110011",  35: "1100001010001",  36: "1100010010100",
  37: "1100010110111",  38: "1100100010001",  39: "1100100011001",  40: "1100110101011",
  41: "1100110110001",  42: "1100111010010",  43: "1101001010101",  44: "1101001110100",
  45: "1101011001011",  46: "1101101010111",  47: "1110000110100",  48: "1110010000011",
  49: "1110010001011",  50: "1110010100011",  51: "1110010101000",  52: "1110100111011",
  53: "1110110010111",  54: "1111001001000",  55: "1111010010100",  56: "1111010011001",
  57: "1111011011010",  58: "1111011111000",  59: "1111011111111",  60: "1111110110101",
  61: "0010000000010",  62: "1101111110101",  63: "0001111010010",
}

def str2list(s):
  x = []
  for c in s:
    if c=='0':
      x.append(0)
    else:
      x.append(1)
  return x

def b2ad_g1_shift(x):
  return [x[0]^x[4]^x[10]^x[12]] + x[0:12]

def b2ad_g2_shift(x):
  return [x[2]^x[4]^x[8]^x[10]^x[11]^x[12]] + x[0:12]

def b2ad(prn):
  n = code_length
  g1 = [1,1,1,1,1,1,1,1,1,1,1,1,1]
  g2 = str2list(b2ad_g2_initial[prn])
  b2ad = np.zeros(n)
  for i in range(n):
    b2ad[i] = g1[12] ^ g2[12]
    if i==8189:
      g1 = [1,1,1,1,1,1,1,1,1,1,1,1,1]
    else:
      g1 = b2ad_g1_shift(g1)
    g2 = b2ad_g2_shift(g2)
  return b2ad

codes = {}

def b2ad_code(prn):
  if prn not in codes:
    codes[prn] = b2ad(prn)
  return codes[prn]

def code(prn,chips,frac,incr,n):
  c = b2ad_code(prn)
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

# test

def chips2octal(x):
  s = ''
  for i in range(len(x)//3):
    d = 4*x[3*i] + 2*x[3*i+1] + x[3*i+2]
    s = s + '%o'%int(d)
  return s

if __name__=='__main__':
  for prn in range(1,64):
    c = b2ad_code(prn)
    s1 = chips2octal(c[0:24])
    s2 = chips2octal(c[-24:])
    print("%d %s %s"%(prn,s1,s2))
