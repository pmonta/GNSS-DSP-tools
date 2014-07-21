#!/usr/bin/env python

import sys
import numpy as np

import gnsstools.galileo.e1b as e1b
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
    self.eml = 0

def costas(x):
  if np.real(x)>0:
    return np.arctan2(np.imag(x),np.real(x))
  else:
    return np.arctan2(-np.imag(x),-np.real(x))

# tracking loops

def track(x,s):
  n = len(x)
  fs = s.fs

  nco.mix(x,-s.carrier_f/fs, s.carrier_p, nco.nco_table)
  s.carrier_p = s.carrier_p - n*s.carrier_f/fs
  s.carrier_p = np.mod(s.carrier_p,1)

  cf = (s.code_f+s.carrier_f/1540.0)/fs

  p_early = e1b.correlate(x, s.prn, 0, s.code_p-0.2, cf, e1b.e1b_code(prn), e1b.boc11)
  p_prompt = e1b.correlate(x, s.prn, 0, s.code_p, cf, e1b.e1b_code(prn), e1b.boc11)
  p_late = e1b.correlate(x, s.prn, 0, s.code_p+0.2, cf, e1b.e1b_code(prn), e1b.boc11)

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
    pll_k1 = 0.03
    pll_k2 = 1.5
    e = costas(p_prompt)
    e1 = s.carrier_e1
    s.carrier_f = s.carrier_f + pll_k1*e + pll_k2*(e-e1)
    s.carrier_e1 = e

# code loop

  dll_k1 = 0.003
  dll_k2 = 1.0
  pwr_early = np.real(p_early*np.conj(p_early))
  pwr_late = np.real(p_late*np.conj(p_late))
  if (pwr_late+pwr_early)==0:
    e = 0
  else:
    e = (pwr_late-pwr_early)/(pwr_late+pwr_early)
  s.eml = e
  e1 = s.code_e1
  s.code_f = s.code_f + dll_k1*e + dll_k2*(e-e1)
  s.code_e1 = e

  s.code_p = s.code_p + n*cf
  s.code_p = np.mod(s.code_p,e1b.code_length)

  return p_prompt,s

#
# main program
#

# parse command-line arguments
# example:
#   ./track-galileo-e1b.py /dev/stdin 68873142.857 -8662285.714 12 -350 1274.2

filename = sys.argv[1]             # input data, raw file, i/q interleaved, 8 bit signed (two's complement)
fs = float(sys.argv[2])            # sampling rate, Hz
coffset = float(sys.argv[3])       # offset to L1 carrier, Hz (positive or negative)
prn = int(sys.argv[4])             # PRN code
doppler = float(sys.argv[5])       # initial doppler estimate from acquisition
code_offset = float(sys.argv[6])   # initial code offset from acquisition

n = int(round(0.001*fs))           # number of samples per block, approx 1 ms
fp = open(filename,"rb")

s = tracking_state(fs=fs, prn=prn,                    # initialize tracking state
  code_p=code_offset, code_f=e1b.chip_rate, code_i=0,
  carrier_p=0, carrier_f=doppler, carrier_i=0,
  mode='PLL')

block = 0
coffset_phase = 0.0

do_plots = False

if do_plots:
  from plotting import stripchart
  s1 = stripchart.stripchart(n=2000)
  s2 = stripchart.stripchart(n=2000)
  s3 = stripchart.stripchart(n=2000)
  s4 = stripchart.stripchart(n=2000)
  s5 = stripchart.stripchart(n=2000)

while True:
  x = io.get_samples_complex(fp,n)
  if x==None:
    break

  nco.mix(x,-coffset/fs,coffset_phase,nco.nco_table)
  coffset_phase = coffset_phase - n*coffset/fs
  coffset_phase = np.mod(coffset_phase,1)

  p_prompt,s = track(x,s)
  print block,np.real(p_prompt),np.imag(p_prompt),s.carrier_f,s.code_f-e1b.chip_rate,s.eml
  if do_plots:
    s1.point(s.carrier_f)
    s2.point(s.code_f)
    s3.point(np.real(p_prompt))
    s4.point(np.imag(p_prompt))
    s5.point(s.eml)

  block = block + 1
  if (block%100)==0:
    sys.stderr.write("%d\n"%block)
