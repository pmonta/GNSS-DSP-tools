#!/usr/bin/env python

import sys
import numpy as np

import gnsstools.gps.ca as ca
import gnsstools.nco as nco
import gnsstools.io as io

class tracking_state:
  def __init__(self,fs,prn,code_p,code_f,code_i,carrier_p,carrier_f,carrier_i,mode):
    self.fs = fs
    self.prn = prn
    self.code_p = code_p
    self.code_f = code_f
    self.code_i = code_i
    self.carrier_p = carrier_p
    self.carrier_f = carrier_f
    self.carrier_i = carrier_i
    self.mode = mode
    self.prompt1 = 0 + 0*(1j)
    self.carrier_e1 = 0
    self.code_e1 = 0

def costas(x):
  if np.real(x)>0:
    return np.arctan2(np.imag(x),np.real(x))
  else:
    return np.arctan2(-np.imag(x),-np.real(x))

# tracking loops

def track(x,s):
  n = len(x)
  fs = s.fs

  x = x * nco.nco(-s.carrier_f/fs, s.carrier_p, n)

  c_early = ca.code(s.prn, s.code_p-0.5, 0, s.code_f/fs, n)
  p_early = np.sum(x*c_early)
  c_prompt = ca.code(s.prn, s.code_p, 0, s.code_f/fs, n)
  p_prompt = np.sum(x*c_prompt)
  c_late = ca.code(s.prn, s.code_p+0.5, 0, s.code_f/fs, n)
  p_late = np.sum(x*c_late)

  if s.mode=='FLL_WIDE':
    fll_k = 2.0
    a = p_prompt
    b = s.prompt1
    e = np.arctan2(np.imag(a)*np.real(b)-np.real(a)*np.imag(b),np.real(a)*np.real(b)+np.imag(a)*np.imag(b))
    s.carrier_f = s.carrier_f + fll_k*e
    s.prompt1 = p_prompt
  elif s.mode=='FLL_NARROW':
    fll_k = 0.3
    a = p_prompt
    b = s.prompt1
    e = np.arctan2(np.imag(a)*np.real(b)-np.real(a)*np.imag(b),np.real(a)*np.real(b)+np.imag(a)*np.imag(b))
    s.carrier_f = s.carrier_f + fll_k*e
    s.prompt1 = p_prompt
  elif s.mode=='PLL':
    pll_k1 = 0.15
    pll_k2 = 6.0
    e = costas(p_prompt)
    e1 = s.carrier_e1
    s.carrier_f = s.carrier_f + pll_k1*e + pll_k2*(e-e1)
    s.carrier_e1 = e

  s.carrier_p = s.carrier_p - n*s.carrier_f/fs
  s.carrier_p = np.mod(s.carrier_p,1)

# code loop

  dll_k1 = 0.005
  dll_k2 = 0.6
  d = np.absolute(p_prompt)
  if d==0:
    e = 0
  else:
    e = (np.absolute(p_late) - np.absolute(p_early))/d
  e1 = s.code_e1
  s.code_f = s.code_f + dll_k1*e + dll_k2*(e-e1)  # fixme: carrier aiding
  s.code_e1 = e

  s.code_p = s.code_p + n*s.code_f/fs
  s.code_p = np.mod(s.code_p,ca.code_length)

  return p_prompt,s

#
# main program
#

# parse command-line arguments
# example:
#   ./track-gps-l1.py data/gps-6002-l1_a.dat 68873142.857 -8662285.714 12 -400 781.2

filename = sys.argv[1]             # input data, raw file, i/q interleaved, 8 bit signed (two's complement)
fs = float(sys.argv[2])            # sampling rate, Hz
coffset = float(sys.argv[3])       # offset to L1 carrier, Hz (positive or negative)
prn = int(sys.argv[4])             # PRN code
doppler = float(sys.argv[5])       # initial doppler estimate from acquisition
code_offset = float(sys.argv[6])   # initial code offset from acquisition

n = int(round(0.001*fs))           # number of samples per block, approx 1 ms
fp = open(filename,"rb")

s = tracking_state(fs=fs, prn=prn,                    # initialize tracking state
  code_p=code_offset, code_f=ca.chip_rate, code_i=0,
  carrier_p=0, carrier_f=doppler, carrier_i=0,
  mode='FLL_WIDE')

block = 0
coffset_phase = 0

T=10000
r_p = np.zeros(T) + (1j)*np.zeros(T)
r_f = np.zeros(T)
r_cf = np.zeros(T)

while True:
  x = io.get_samples_complex(fp,n)
  if x==None:
    break

  w = nco.nco(-coffset/fs,coffset_phase,n)
  coffset_phase = coffset_phase - n*coffset/fs
  coffset_phase = np.mod(coffset_phase,1)
  x = x*w

  p_prompt,s = track(x,s)
#  print block,p_prompt,s.carrier_f,s.code_f
  r_p[block] = p_prompt
  r_f[block] = s.carrier_f
  r_cf[block] = s.code_f

  block = block + 1
  if (block%100)==0:
    print block
  if block==1000:
    s.mode = 'FLL_NARROW'
  if block==2000:
    s.mode = 'PLL'
  if block==T:
    break
