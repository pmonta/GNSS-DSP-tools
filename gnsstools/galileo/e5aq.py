# Galileo E5a-Q code construction
#
# Copyright 2014 Peter Monta

import numpy as np

chip_rate = 10230000
code_length = 10230

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

# initial-state table from Table 16 (page 15) of Galileo ICD SIS (2014)
# index is PRN

e5aq_init = {
   1: 0o25652,    2: 0o05142,    3: 0o24723,    4: 0o31751,
   5: 0o27366,    6: 0o24660,    7: 0o33655,    8: 0o27450,
   9: 0o07626,   10: 0o01705,   11: 0o12717,   12: 0o32122,
  13: 0o16075,   14: 0o16644,   15: 0o37556,   16: 0o02477,
  17: 0o02265,   18: 0o06430,   19: 0o25046,   20: 0o12735,
  21: 0o04262,   22: 0o11230,   23: 0o00037,   24: 0o06137,
  25: 0o04312,   26: 0o20606,   27: 0o11162,   28: 0o22252,
  29: 0o30533,   30: 0o24614,   31: 0o07767,   32: 0o32705,
  33: 0o05052,   34: 0o27553,   35: 0o03711,   36: 0o02041,
  37: 0o34775,   38: 0o05274,   39: 0o37356,   40: 0o16205,
  41: 0o36270,   42: 0o06600,   43: 0o26773,   44: 0o17375,
  45: 0o35267,   46: 0o36255,   47: 0o12044,   48: 0o26442,
  49: 0o21621,   50: 0o25411
}

def e5aq_reg1_shift(x):
  return [x[13]^x[7]^x[5]^x[0]] + x[0:13]

def e5aq_reg2_shift(x):
  return [x[13]^x[11]^x[7]^x[6]^x[4]^x[3]] + x[0:13]

def make_e5aq_reg1():
  x = [1,1,1,1,1,1,1,1,1,1,1,1,1,1]
  y = np.zeros(code_length)
  for i in range(code_length):
    y[i] = x[13]
    x = e5aq_reg1_shift(x)
  return y

def make_e5aq_reg2(start):
  x = start
  y = np.zeros(code_length)
  for i in range(code_length):
    y[i] = x[13]
    x = e5aq_reg2_shift(x)
  return y

r1 = make_e5aq_reg1()

def seq(a):
  s = []
  for i in range(14):
    s = s + [(a>>i)&1]
  return s

codes = {}

def make_e5aq(prn):
  start = seq(e5aq_init[prn])
  r2 = make_e5aq_reg2(start)
  return np.logical_xor(r1,r2)

def e5aq_code(prn):
  if prn not in codes:
    codes[prn] = make_e5aq(prn)
  return codes[prn]

def code(prn,chips,frac,incr,n):
  c = e5aq_code(prn)
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

# test vectors in Galileo ICD

if __name__=='__main__':
  for prn in [1,2,3,4]:
    x = e5aq_code(prn)
    print(x[0:24].astype('int'))
