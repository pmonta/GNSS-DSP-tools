# Galileo E5b-Q code construction
#
# Copyright 2014 Peter Monta

import numpy as np

chip_rate = 10230000
code_length = 10230

# secondary code table from Table 20 (page 20) of Galileo ICD SIS (2014)
# index is PRN

secondary_code = {
  1: 'CFF914EE3C6126A49FD5E5C94',   2: 'FC317C9A9BF8C6038B5CADAB3',
  3: 'A2EAD74B6F9866E414393F239',   4: '72F2B1180FA6B802CB84DF997',
  5: '13E3AE93BC52391D09E84A982',   6: '77C04202B91B22C6D3469768E',
  7: 'FEBC592DD7C69AB103D0BB29C',   8: '0B494077E7C66FB6C51942A77',
  9: 'DD0E321837A3D52169B7B577C',  10: '43DEA90EA6C483E7990C3223F',
 11: '0366AB33F0167B6FA979DAE18',  12: '99CCBBFAB1242CBE31E1BD52D',
 13: 'A3466923CEFDF451EC0FCED22',  14: '1A5271F22A6F9A8D76E79B7F0',
 15: '3204A6BB91B49D1A2D3857960',  16: '32F83ADD43B599CBFB8628E5B',
 17: '3871FB0D89DB77553EB613CC1',  18: '6A3CBDFF2D64D17E02773C645',
 19: '2BCD09889A1D7FC219F2EDE3B',  20: '3E49467F4D4280B9942CD6F8C',
 21: '658E336DCFD9809F86D54A501',  22: 'ED4284F345170CF77268C8584',
 23: '29ECCE910D832CAF15E3DF5D1',  24: '456CCF7FE9353D50E87A708FA',
 25: 'FB757CC9E18CBC02BF1B84B9A',  26: '5686229A8D98224BC426BC7FC',
 27: '700A2D325EA14C4B7B7AA8338',  28: '1210A330B4D3B507D854CBA3F',
 29: '438EE410BD2F7DBCDD85565BA',  30: '4B9764CC455AE1F61F7DA432B',
 31: 'BF1F45FDDA3594ACF3C4CC806',  32: 'DA425440FE8F6E2C11B8EC1A4',
 33: 'EE2C8057A7C16999AFA33FED1',  34: '2C8BD7D8395C61DFA96243491',
 35: '391E4BB6BC43E98150CDDCADA',  36: '399F72A9EADB42C90C3ECF7F0',
 37: '93031FDEA588F88E83951270C',  38: 'BA8061462D873705E95D5CB37',
 39: 'D24188F88544EB121E963FD34',  40: 'D5F6A8BB081D8F383825A4DCA',
 41: '0FA4A205F0D76088D08EAF267',  42: '272E909FAEBC65215E263E258',
 43: '3370F35A674922828465FC816',  44: '54EF96116D4A0C8DB0E07101F',
 45: 'DE347C7B27FADC48EF1826A2B',  46: '01B16ECA6FC343AE08C5B8944',
 47: '1854DB743500EE94D8FC768ED',  48: '28E40C684C87370CD0597FAB4',
 49: '5E42C19717093353BCAAF4033',  50: '64310BAD8EB5B36E38646AF01'
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

e5bq_init = {
   1: 0o03331,    2: 0o06143,    3: 0o25322,    4: 0o23371,
   5: 0o00413,    6: 0o36235,    7: 0o17750,    8: 0o04745,
   9: 0o13005,   10: 0o37140,   11: 0o30155,   12: 0o20237,
  13: 0o03461,   14: 0o31662,   15: 0o27146,   16: 0o05547,
  17: 0o02456,   18: 0o30013,   19: 0o00322,   20: 0o10761,
  21: 0o26767,   22: 0o36004,   23: 0o30713,   24: 0o07662,
  25: 0o21610,   26: 0o20134,   27: 0o11262,   28: 0o10706,
  29: 0o34143,   30: 0o11051,   31: 0o25460,   32: 0o17665,
  33: 0o32354,   34: 0o21230,   35: 0o20146,   36: 0o11362,
  37: 0o37246,   38: 0o16344,   39: 0o15034,   40: 0o25471,
  41: 0o25646,   42: 0o22157,   43: 0o04336,   44: 0o16356,
  45: 0o04075,   46: 0o02626,   47: 0o11706,   48: 0o37011,
  49: 0o27041,   50: 0o31024
}

def e5bq_reg1_shift(x):
  return [x[13]^x[12]^x[10]^x[3]] + x[0:13]

def e5bq_reg2_shift(x):
  return [x[13]^x[9]^x[8]^x[5]^x[4]^x[0]] + x[0:13]

def make_e5bq_reg1():
  x = [1,1,1,1,1,1,1,1,1,1,1,1,1,1]
  y = np.zeros(code_length)
  for i in range(code_length):
    y[i] = x[13]
    x = e5bq_reg1_shift(x)
  return y

def make_e5bq_reg2(start):
  x = start
  y = np.zeros(code_length)
  for i in range(code_length):
    y[i] = x[13]
    x = e5bq_reg2_shift(x)
  return y

r1 = make_e5bq_reg1()

def seq(a):
  s = []
  for i in range(14):
    s = s + [(a>>i)&1]
  return s

codes = {}

def make_e5bq(prn):
  start = seq(e5bq_init[prn])
  r2 = make_e5bq_reg2(start)
  return np.logical_xor(r1,r2)

def e5bq_code(prn):
  if prn not in codes:
    codes[prn] = make_e5bq(prn)
  return codes[prn]

def code(prn,chips,frac,incr,n):
  c = e5bq_code(prn)
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
    x = e5bq_code(prn)
    print(x[0:24].astype('int'))
