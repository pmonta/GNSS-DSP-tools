#!/usr/bin/env python

import sys
import numpy as np

import gnsstools.gps.l2cl as l2cl
import gnsstools.nco as nco
import gnsstools.io as io

#
# Acquisition search
#

def search(x,prn,doppler,l2cm_code_phase):
  n = int(fs*0.020)
  w = nco.nco(-doppler/fs,0,n)
  incr = l2cl.chip_rate/fs
  m_metric,m_k = 0,0
  for k in range(75):
    q = 0
    for block in range(4):
      c = l2cl.code(prn,(k+block)*10230+l2cm_code_phase,0,incr,n)
      p = x[n*block:n*(block+1)]*c*w
      q = q + np.absolute(np.sum(p))
    if q>m_metric:
      m_metric = q
      m_k = k
  return m_metric,m_k

#
# main program
#

# parse command-line arguments
# example:
#   ./acquire-gps-l2cl.py /dev/stdin 68873142.857 -12116571.429 31 -620 141.2

filename = sys.argv[1]        # input data, raw file, i/q interleaved, 8 bit signed (two's complement)
fs = float(sys.argv[2])       # sampling rate, Hz
coffset = float(sys.argv[3])  # offset to L1 carrier, Hz (positive or negative)
prn = int(sys.argv[4])
doppler = float(sys.argv[5])
l2cm_code_phase = float(sys.argv[6])

# read first 85 ms of file

n = int(fs*0.085)
fp = open(filename,"rb")
x = io.get_samples_complex(fp,n)

nco.mix(x,-coffset/fs,0)

metric,k = search(x,prn,doppler,l2cm_code_phase)
print('%f %f'%(10230*k+l2cm_code_phase,metric))
