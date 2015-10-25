# GPS P code construction
#
# Copyright 2014 Peter Monta

import numpy as np

chip_rate = 10230000
code_length = chip_rate*86400*7

def x1a_shift(x):
  return [x[11]^x[10]^x[7]^x[5]] + x[0:11]

def x1b_shift(x):
  return [x[11]^x[10]^x[9]^x[8]^x[7]^x[4]^x[1]^x[0]] + x[0:11]

def x2a_shift(x):
  return [x[11]^x[10]^x[9]^x[8]^x[7]^x[6]^x[4]^x[3]^x[2]^x[0]] + x[0:11]

def x2b_shift(x):
  return [x[11]^x[8]^x[7]^x[3]^x[2]^x[1]] + x[0:11]

def make_12bit(reg_shift,reg_initial,n):
  x = reg_initial
  s = np.zeros(n)
  for i in range(n):
    s[i] = x[11]
    x = reg_shift(x)
  return s

x1a = make_12bit(x1a_shift,[0,0,0,1,0,0,1,0,0,1,0,0],4092)
x1b = make_12bit(x1b_shift,[0,0,1,0,1,0,1,0,1,0,1,0],4093)
x2a = make_12bit(x2a_shift,[1,0,1,0,0,1,0,0,1,0,0,1],4092)
x2b = make_12bit(x2b_shift,[0,0,1,0,1,0,1,0,1,0,1,0],4093)

def x1(prn,start,len):
  idx = start + np.arange(len)
  idx = idx % 15345000

  p_x1a = x1a[idx % 4092]

  hold = idx>=(15345037-343)
  idx[hold] = 4092
  p_x1b = x1b[idx % 4093]

  return np.logical_xor(p_x1a,p_x1b)

def x2(prn,start,len):
  idx = start + np.arange(len)
  idx = idx % 15345037

  idx_a = idx
  hold = idx_a>=(15345037-37)
  idx_a[hold] = 4091
  p_x2a = x2a[idx_a % 4092]

  idx_b = idx
  hold = idx_b>=(15345037-37-343)
  idx_b[hold] = 4092
  p_x2b = x2b[idx_b % 4093]

  return np.logical_xor(p_x2a,p_x2b)

def last_x2(prn,start,len):
  idx = start + np.arange(len)
  idx = idx % 15345037

  idx_a = idx
  hold = idx_a>=(15345037-1069)
  idx_a[hold] = 4091
  p_x2a = x2a[idx_a % 4092]

  idx_b = idx
  hold = idx_b>=(15345037-965)
  idx_b[hold] = 4092
  p_x2b = x2b[idx_b % 4093]

  return np.logical_xor(p_x2a,p_x2b)

def x2_delay(prn):
  day = (prn-1)//37
  prn = prn - 37*day
  return -prn + chip_rate*86400*day

def p_code(prn,start,len):
#  return x1(prn,start,len) ^ x2(prn,start+x2_delay(prn),len)

  p_x1 = x1(prn,start,len)

  p_x2 = x2(prn,start+x2_delay(prn),len)

  p_last_x2 = last_x2(prn,start+x2_delay(prn),len)
  idx_x2 = start + x2_delay(prn) + np.arange(len)
  idx_x2 = idx_x2 % code_length
  idx_last_x2 = idx_x2>=(code_length-4092)
  p_x2[idx_last_x2] = p_last_x2[idx_last_x2]

  return np.logical_xor(p_x1,p_x2)

def code(prn,chips,frac,incr,n):
  len = np.int(np.floor(n*incr)+5)
  c = p_code(prn,chips,len).astype('int')
  idx = frac + incr*np.arange(n)
  idx = np.floor(idx).astype('int')
  x = c[idx]
  return 1.0 - 2.0*x

# test vectors in IS-GPS-200H

def first_12_chips(prn):
  c = p_code(prn,0,12)
  r = 0
  for i in range(12):
    r = 2*r + c[i]
  return r

def print_first_12_chips():
  for prn in range(1,211):
    print('%d: %04o' % (prn,first_12_chips(prn)))

if __name__=='__main__':
  print_first_12_chips()
