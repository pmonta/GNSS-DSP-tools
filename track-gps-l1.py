#!/usr/bin/env python

import sys
import numpy as np

import gnsstools.gps.ca as ca
import gnsstools.nco as nco
import gnsstools.io as io
import gnsstools.discriminator as discriminator

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

# tracking loops

def track(x,s):
  n = len(x)
  fs = s.fs

  nco.mix(x,-s.carrier_f/fs, s.carrier_p, nco.nco_table)
  s.carrier_p = s.carrier_p - n*s.carrier_f/fs
  s.carrier_p = np.mod(s.carrier_p,1)

  cf = (s.code_f+s.carrier_f/1540.0)/fs

  p_early = ca.correlate(x, s.prn, 0, s.code_p-0.5, cf, ca.ca_code(prn))
  p_prompt = ca.correlate(x, s.prn, 0, s.code_p, cf, ca.ca_code(prn))
  p_late = ca.correlate(x, s.prn, 0, s.code_p+0.5, cf, ca.ca_code(prn))

  if s.mode=='FLL_WIDE':
    fll_k = 3.0
    a = p_prompt
    b = s.prompt1
    e = discriminator.fll_atan2(a,b)
    s.carrier_f = s.carrier_f + fll_k*e
    s.prompt1 = p_prompt
  elif s.mode=='FLL_NARROW':
    fll_k = 0.3
    a = p_prompt
    b = s.prompt1
    e = discriminator.fll_atan2(a,b)
    s.carrier_f = s.carrier_f + fll_k*e
    s.prompt1 = p_prompt
  elif s.mode=='PLL':
    pll_k1 = 0.1
    pll_k2 = 5.0
    e = discriminator.pll_costas(p_prompt)
    e1 = s.carrier_e1
    s.carrier_f = s.carrier_f + pll_k1*e + pll_k2*(e-e1)
    s.carrier_e1 = e

# code loop

  dll_k1 = 0.0005
  dll_k2 = 0.2
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

fp = open(filename,"rb")

n = int(fs*0.001*((ca.code_length-code_offset)/ca.code_length))  # align with 1 ms code boundary
x = io.get_samples_complex(fp,n)
code_offset += n*1000.0*ca.code_length/fs

s = tracking_state(fs=fs, prn=prn,                    # initialize tracking state
  code_p=code_offset, code_f=ca.chip_rate, code_i=0,
  carrier_p=0, carrier_f=doppler, carrier_i=0,
  mode='FLL_WIDE')

block = 0
coffset_phase = 0.0

while True:
  if s.code_p<ca.code_length/2:
    n = int(fs*0.001*(ca.code_length-s.code_p)/ca.code_length)
  else:
    n = int(fs*0.001*(2*ca.code_length-s.code_p)/ca.code_length)

  x = io.get_samples_complex(fp,n)
  if x==None:
    break

  nco.mix(x,-coffset/fs,coffset_phase,nco.nco_table)
  coffset_phase = coffset_phase - n*coffset/fs
  coffset_phase = np.mod(coffset_phase,1)

  p_prompt,s = track(x,s)
  print block,np.real(p_prompt),np.imag(p_prompt),s.carrier_f,s.code_f-ca.chip_rate,(180/np.pi)*np.angle(p_prompt)

  block = block + 1
  if (block%100)==0:
    sys.stderr.write("%d\n"%block)
  if block==500:
    s.mode = 'FLL_NARROW'
  if block==1000:
    s.mode = 'PLL'
