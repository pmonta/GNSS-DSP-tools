#!/usr/bin/env python

import sys
import numpy as np

import gnsstools.glonass.p as p
import gnsstools.nco as nco
import gnsstools.io as io

#
# Acquisition search
#

def search(x,chan,doppler,ca_code_phase):
  n = int(fs*0.004)
  w = nco.nco(-(437500*chan+doppler)/fs,0,n)
  m_metric,m_k = 0,0
  for k in range(1000):
    q = 0
    cp = 5110*k + 10*ca_code_phase
    for block in range(20):
      incr = 5110000.0/fs
      c = p.code(0,cp,incr,n)
      xp = x[n*block:n*(block+1)]*c*w
      q = q + np.absolute(np.sum(xp))
      cp += n*incr
    print('%f %f'%(k,q))
    if q>m_metric:
      m_metric = q
      m_k = k
  return m_metric,m_k

#
# main program
#

# parse command-line arguments
# example:
#   ./acquire-glonass-l2-p.py /dev/stdin 68873142.857 6283428.571 1 400 183.03

filename = sys.argv[1]        # input data, raw file, i/q interleaved, 8 bit signed (two's complement)
fs = float(sys.argv[2])       # sampling rate, Hz
coffset = float(sys.argv[3])  # offset to L1 GLONASS carrier channel 0, Hz (positive or negative)
chan = int(sys.argv[4])
doppler = float(sys.argv[5])
ca_code_phase = float(sys.argv[6])

# read first 85 ms of file

n = int(fs*0.085)
fp = open(filename,"rb")
x = io.get_samples_complex(fp,n)

nco.mix(x,-coffset/fs,0,nco.nco_table)

metric,k = search(x,chan,doppler,ca_code_phase)
print('%f %f'%(5110*k+10*ca_code_phase,metric))
