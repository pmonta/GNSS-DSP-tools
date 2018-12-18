# Beidou B2bi code construction
#
# Copyright 2018 Peter Monta

import numpy as np

from .b2bi_strings import *

chip_rate = 10230000
code_length = 10230

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

def b2bi_parse_base64(prn):
  s = b2bi_strings[prn]
  n = code_length
  y = np.zeros(n)
  for i in range(n):
    idx = i//6
    bit = 5-(i%6)
    x = b64[s[idx]]
    y[i] = (x>>bit)&1
  return y

codes = {}

def b2bi_code(prn):
  if prn not in codes:
    codes[prn] = b2bi_parse_base64(prn)
  return codes[prn]

def code(prn,chips,frac,incr,n):
  c = b2bi_code(prn)
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

@jit(nopython=True)
def accum(x, cp, incr, a):
  n = len(x)
  for i in range(n):
    idx = int(cp)
    a[idx] += x[i]
    cp = (cp+incr)%code_length
