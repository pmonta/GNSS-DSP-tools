#!/usr/bin/env python

import sys
import os
import numpy as np
import scipy.signal
import scipy.fftpack as fft

import gnsstools.glonass.ca as ca
import gnsstools.nco as nco
import gnsstools.io as io

#
# Acquisition search
#

def search(x,chan):
  fs = 68873142.857
  n = 68875
  incr = float(ca.code_length)/n
  c = ca.code(0,0,incr,n)                          # obtain samples of the C/A code
  c = fft.fft(c)
  m_metric,m_code,m_doppler = 0,0,0
  for doppler in np.arange(-5000,5000,200):        # doppler bins
    q = np.zeros(n)
    w = nco.nco(-(437500*chan+doppler)/fs,0,n)
    for block in range(20):                        # 20 incoherent sums
      b = x[(block*n):((block+1)*n)]
      b = b*w
      r = fft.ifft(c*np.conj(fft.fft(b)))
      q = q + np.absolute(r)
    idx = np.argmax(q)
    if q[idx]>m_metric:
      m_metric = q[idx]
      m_code = ca.code_length*(float(idx)/n)
      m_doppler = doppler
  return m_metric,m_code,m_doppler

#
# main program
#

# parse command-line arguments
# example:
#   ./acquire-glonass-l2.py data/gps-5001-l1_a.dat 68873142.857 6283428.571

filename = sys.argv[1]        # input data, raw file, i/q interleaved, 8 bit signed (two's complement)
fs = float(sys.argv[2])       # sampling rate, Hz
coffset = float(sys.argv[3])  # offset to L1 GLONASS carrier channel 0, Hz (positive or negative)

# read first 25 ms of file

n = int(fs*0.025)
fp = open(filename,"rb")
x = io.get_samples_complex(fp,n)

nco.mix(x,-coffset/fs,0,nco.nco_table)

# iterate over channels of interest

for chan in range(-7,7):
  metric,code,doppler = search(x,chan)
  if metric>0.0:    # fixme: need a proper metric and threshold; and estimate cn0
    print 'chan % 2d doppler % 7.1f metric %7.1f code_offset %6.1f' % (chan,doppler,metric,code)
