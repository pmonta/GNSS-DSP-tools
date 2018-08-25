# Galileo E6-C code construction
#
# Copyright 2014 Peter Monta

import numpy as np

from .e6c_strings import *

chip_rate = 5115000
code_length = 5115

# secondary code table from Table 19 (page 18) of Galileo ICD SIS (2014)
# index is PRN

secondary_code = {
  1: '83F6F69D8F6E15411FB8C9B1C',   2: '66558BD3CE0C7792E83350525',
  3: '59A025A9C1AF0651B779A8381',   4: 'D3A32640782F7B18E4DF754B7',
  5: 'B91FCAD7760C218FA59348A93',   6: 'BAC77E933A779140F094FBF98',
  7: '537785DE280927C6B58BA6776',   8: 'EFCAB4B65F38531ECA22257E2',
  9: '79F8CAE838475EA5584BEFC9B',  10: 'CA5170FEA3A810EC606B66494',
 11: '1FC32410652A2C49BD845E567',  12: 'FE0A9A7AFDAC44E42CB95D261',
 13: 'B03062DC2B71995D5AD8B7DBE',  14: 'F6C398993F598E2DF4235D3D5',
 15: '1BB2FB8B5BF24395C2EF3C5A1',  16: '2F920687D238CC7046EF6AFC9',
 17: '34163886FC4ED7F2A92EFDBB8',  18: '66A872CE47833FB2DFD5625AD',
 19: '99D5A70162C920A4BB9DE1CA8',  20: '81D71BD6E069A7ACCBEDC66CA',
 21: 'A654524074A9E6780DB9D3EC6',  22: 'C3396A101BEDAF623CFC5BB37',
 23: 'C3D4AB211DF36F2111F2141CD',  24: '3DFF25EAE761739265AF145C1',
 25: '994909E0757D70CDE389102B5',  26: 'B938535522D119F40C25FDAEC',
 27: 'C71AB549C0491537026B390B7',  28: '0CDB8C9E7B53F55F5B0A0597B',
 29: '61C5FA252F1AF81144766494F',  30: '626027778FD3C6BB4BAA7A59D',
 31: 'E745412FF53DEBD03F1C9A633',  32: '3592AC083F3175FA724639098',
 33: '52284D941C3DCAF2721DDB1FD',  34: '73B3D8F0AD55DF4FE814ED890',
 35: '94BF16C83BD7462F6498E0282',  36: 'A8C3DE1AC668089B0B45B3579',
 37: 'E23FFC2DD2C14388AD8D6BEC8',  38: 'F2AC871CDF89DDC06B5960D2B',
 39: '06191EC1F622A77A526868BA1',  40: '22D6E2A768E5F35FFC8E01796',
 41: '25310A06675EB271F2A09EA1D',  42: '9F7993C621D4BEC81A0535703',
 43: 'D62999EACF1C99083C0B4A417',  44: 'F665A7EA441BAA4EA0D01078C',
 45: '46F3D3043F24CDEABD6F79543',  46: 'E2E3E8254616BD96CEFCA651A',
 47: 'E548231A82F9A01A19DB5E1B2',  48: '265C7F90A16F49EDE2AA706C8',
 49: '364A3A9EB0F0481DA0199D7EA',  50: '9810A7A898961263A0F749F56'
}

# parse the hex string, return a 100-bit secondary-code array

def secondary_seq(s):
  x = np.zeros(100)
  for i in range(100):
    nib = i//4
    bit = 3-(i%4)
    x[i] = (int(s[nib],16)>>bit)&1
    x[i] = 1.0 - 2.0*x[i]
  return x

# transform the secondary-code table entries to sequences

for i in range(1,51):
  secondary_code[i] = secondary_seq(secondary_code[i])

b64 = {
  'A':0, 'B':1, 'C':2, 'D':3, 'E':4, 'F':5, 'G':6, 'H':7,
  'I':8, 'J':9, 'K':10, 'L':11, 'M':12, 'N':13, 'O':14, 'P':15,
  'Q':16, 'R':17, 'S':18, 'T':19, 'U':20, 'V':21, 'W':22, 'X':23,
  'Y':24, 'Z':25, 'a':26, 'b':27, 'c':28, 'd':29, 'e':30, 'f':31,
  'g':32, 'h':33, 'i':34, 'j':35, 'k':36, 'l':37, 'm':38, 'n':39,
  'o':40, 'p':41, 'q':42, 'r':43, 's':44, 't':45, 'u':46, 'v':47,
  'w':48, 'x':49, 'y':50, 'z':51, '0':52, '1':53, '2':54, '3':55,
  '4':56, '5':57, '6':58, '7':59, '8':60, '9':61, '+':62, '/':63,
}

def e6c_parse_base64(prn):
  s = e6c_strings[prn]
  n = code_length
  y = np.zeros(n)
  for i in range(n):
    idx = i//6
    bit = 5-(i%6)
    x = b64[s[idx]]
    y[i] = (x>>bit)&1
  return y

codes = {}

def e6c_code(prn):
  if prn not in codes:
    codes[prn] = e6c_parse_base64(prn)
  return codes[prn]

def code(prn,chips,frac,incr,n):
  c = e6c_code(prn)
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
  print(e6c_code(1)[0:20])
