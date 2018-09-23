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

  hold = idx>=(15345000-343)
  idx[hold] = 4092
  p_x1b = x1b[idx % 4093]

  return np.logical_xor(p_x1a,p_x1b)

def x2(prn,start,len):
  idx = start + np.arange(len)
  idx = idx % 15345037

  idx_a = idx.copy()
  hold = idx_a>=(15345037-37)
  idx_a[hold] = 4091
  p_x2a = x2a[idx_a % 4092]

  idx_b = idx.copy()
  hold = idx_b>=(15345037-37-343)
  idx_b[hold] = 4092
  p_x2b = x2b[idx_b % 4093]

  return np.logical_xor(p_x2a,p_x2b)

def last_x2(prn,start,len):
  idx = start + np.arange(len)
  idx_x2 = idx % 15345037

  idx_a = idx % 15345000
  hold = idx_a>=(15345000-1069)
  idx_x2_a = idx_x2.copy()
  idx_x2_a[hold] = 4091
  p_x2a = x2a[idx_x2_a % 4092]

  idx_b = idx % 15345000
  hold = idx_b>=(15345000-965)
  idx_x2_b = idx_x2.copy()
  idx_x2_b[hold] = 4092
  p_x2b = x2b[idx_x2_b % 4093]

  return np.logical_xor(p_x2a,p_x2b)

def p_code(prn,start,len):
  day = (prn-1)//37
  prn = prn - 37*day
  start += chip_rate*86400*day
  start = start%code_length

  p_x1 = x1(prn,start,len)
  p_x2 = x2(prn,start-prn,len)

  p_last_x2 = last_x2(prn,(start-prn)%code_length,len)
  idx_x2 = (start - prn + np.arange(len)) % code_length
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

# test vectors in IS-GPS-200J

def first_12_chips(prn):
  c = p_code(prn,0,12)
  r = 0
  for i in range(12):
    r = 2*r + c[i]
  return r

def binary(s):
  n = s.shape[0]
  r = 0
  for i in range(n):
    r += s[i]*2**((n-1)-i)
  return r

def hex_digit(s):
  return "%x"%binary(s)

def chips2hex(c):
  n = c.shape[0]//4
  r = c.shape[0]%4
  if r!=0:
    s = octal_digit(c[0:r])
  else:
    s = ""
  for i in range(n):
    s += hex_digit(c[r+4*i:r+4*(i+1)])
  return s

def first_256_chips_hex(prn):
  start = 0
  len = 256
  c = p_code(prn,start,len)
  return chips2hex(c)

def last_1024_chips_hex(prn):
  start = (7*86400*10230000) - 1024
  len = 1024
  c = p_code(prn,start,len)
  return chips2hex(c)

def print_first_12_chips():
  for prn in range(1,211):
    print('%3d: %04o' % (prn,first_12_chips(prn)))

def print_first_256_chips():
  for prn in [1,2,37,38,74,75,111,112,148,149,185,186,210]:
    print('%3d: %s' % (prn,first_256_chips_hex(prn)))

def print_last_1024_chips():
  for prn in [1,2,37,38,74,75,111,112,148,149,185,186,210]:
    print('%3d: %s' % (prn,last_1024_chips_hex(prn)))

if __name__=='__main__':
  print_first_12_chips()
  print_first_256_chips()
  print_last_1024_chips()
