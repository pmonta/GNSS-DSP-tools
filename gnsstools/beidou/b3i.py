# Beidou B3I code construction
#
# Copyright 2018 Peter Monta

import numpy as np

chip_rate = 10230000
code_length = 10230

secondary_code = np.array([0,0,0,0,0,1,0,0,1,1,0,1,0,1,0,0,1,1,1,0])
secondary_code = 1.0 - 2.0*secondary_code

b3i_g2_initial = {
   1: "1010111111111",   2: "1111000101011",   3: "1011110001010",   4: "1111111111011",
   5: "1100100011111",   6: "1001001100100",   7: "1111111010010",   8: "1110111111101",
   9: "1010000000010",  10: "0010000011011",  11: "1110101110000",  12: "0010110011110",
  13: "0110010010101",  14: "0111000100110",  15: "1000110001001",  16: "1110001111100",
  17: "0010011000101",  18: "0000011101100",  19: "1000101010111",  20: "0001011011110",
  21: "0010000101101",  22: "0010110001010",  23: "0001011001111",  24: "0011001100010",
  25: "0011101001000",  26: "0100100101001",  27: "1011011010011",  28: "1010111100010",
  29: "0001011110101",  30: "0111111111111",  31: "0110110001111",  32: "1010110001001",
  33: "1001010101011",  34: "1100110100101",  35: "1101001011101",  36: "1111101110100",
  37: "0010101100111",  38: "1110100010000",  39: "1101110010000",  40: "1101011001110",
  41: "1000000110100",  42: "0101111011001",  43: "0110110111100",  44: "1101001110001",
  45: "0011100100010",  46: "0101011000101",  47: "1001111100110",  48: "1111101001000",
  49: "0000101001001",  50: "1000010101100",  51: "1111001001100",  52: "0100110001111",
  53: "0000000011000",  54: "1000000000100",  55: "0011010100110",  56: "1011001000110",
  57: "0111001111000",  58: "0010111001010",  59: "1100111110110",  60: "1001001000101",
  61: "0111000100000",  62: "0011001000010",  63: "0010001001110",
}

def str2list(s):
  x = []
  for c in s:
    if c=='0':
      x.append(0)
    else:
      x.append(1)
  return x

def b3i_g1_shift(x):
  if x==[1,1,1,1,1,1,1,1,1,1,1,0,0]:
    return [1,1,1,1,1,1,1,1,1,1,1,1,1]
  else:
    return [x[0]^x[2]^x[3]^x[12]] + x[0:12]

def b3i_g2_shift(x):
  return [x[0]^x[4]^x[5]^x[6]^x[8]^x[9]^x[11]^x[12]] + x[0:12]

def b3i(prn):
  n = code_length
  g1 = [1,1,1,1,1,1,1,1,1,1,1,1,1]
  g2 = str2list(b3i_g2_initial[prn])
  b3i = np.zeros(n)
  for i in range(n):
    b3i[i] = g1[12] ^ g2[12]
    g1 = b3i_g1_shift(g1)
    g2 = b3i_g2_shift(g2)
  return b3i

codes = {}

def b3i_code(prn):
  if prn not in codes:
    codes[prn] = b3i(prn)
  return codes[prn]

def code(prn,chips,frac,incr,n):
  c = b3i_code(prn)
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

if __name__=='__main__':
  print(b3i_code(1)[0:20])
  print(b3i_code(2)[0:20])
