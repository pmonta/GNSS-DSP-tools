#!/usr/bin/env python

import sys
import os
import numpy as np
import scipy.signal
import scipy.fftpack as fft

import gnsstools.galileo.e1c as e1c
import gnsstools.nco as nco
import gnsstools.io as io

#
# Acquisition search
#

def search(x,prn):
  fs = 8184000.0
  n = 32736                                        # 4 ms coherent integration
  incr = float(e1c.code_length)/n
  c = e1c.code(prn,0,0,incr,n)                     # obtain samples of the E1-C code
  boc = nco.boc11(0,0,incr,n)
  c = fft.fft(c*boc)
  m_metric,m_code,m_doppler = 0,0,0
  for doppler in np.arange(-4000,4000,50):         # doppler bins
    q = np.zeros(n)
    w = nco.nco(-doppler/fs,0,n)
    for block in range(20):                        # 20 incoherent sums
      b = x[(block*n):((block+1)*n)]
      b = b*w
      r = fft.ifft(c*np.conj(fft.fft(b)))
      q = q + np.absolute(r)
    idx = np.argmax(q)
    if q[idx]>m_metric:
      m_metric = q[idx]
      m_code = e1c.code_length*(float(idx)/n)
      m_doppler = doppler
  return m_metric,m_code,m_doppler

#
# main program
#

# parse command-line arguments
# example:
#   ./acquire-galileo-e1c.py data/gps-5001-l1_a.dat 68873142.857 -8662285.714

filename = sys.argv[1]        # input data, raw file, i/q interleaved, 8 bit signed (two's complement)
fs = float(sys.argv[2])       # sampling rate, Hz
coffset = float(sys.argv[3])  # offset to E1 Galileo carrier, Hz (positive or negative)

# read first 90 ms of file

n = int(fs*0.090)
fp = open(filename,"rb")
x = io.get_samples_complex(fp,n)

# resample to 8.184 MHz

fsr = 8184000.0/fs
nco.mix(x,-coffset/fs,0,nco.nco_table)
h = scipy.signal.firwin(81,3e6/(fs/2),window='hanning')
x = scipy.signal.filtfilt(h,[1],x)
xr = np.interp((1/fsr)*np.arange(90*8184),np.arange(len(x)),np.real(x))
xi = np.interp((1/fsr)*np.arange(90*8184),np.arange(len(x)),np.imag(x))
x = xr+(1j)*xi

# iterate over channels of interest

#for prn in range(1,51):
for prn in [11,12,19,20]:
  metric,code,doppler = search(x,prn)
  if metric>0.0:    # fixme: need a proper metric and threshold; and estimate cn0
    print 'prn %2d doppler % 7.1f metric %7.1f code_offset %6.1f' % (prn,doppler,metric,code)
