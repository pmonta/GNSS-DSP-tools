#!/usr/bin/env python

import sys
import os
import numpy as np
import scipy.signal
import scipy.fftpack as fft

import gnsstools.gps.l2cm as l2cm
import gnsstools.nco as nco
import gnsstools.io as io

#
# Acquisition search
#

def search(x,prn):
  fs = 4096000.0
  n = 81920                                        # 20 ms coherent integration
  incr = float(l2cm.code_length)/n
  c = l2cm.code(prn,0,0,incr,n)                    # obtain samples of the L2CM code
  c = fft.fft(np.concatenate((c,np.zeros(n))))
  m_metric,m_code,m_doppler = 0,0,0
  for doppler in np.arange(-5000,5000,20):         # doppler bins
    q = np.zeros(2*n)
    w = nco.nco(-doppler/fs,0,2*n)
    for block in range(2):                         # 2 incoherent sums
      b = x[(block*n):((block+2)*n)]
      b = b*w
      r = fft.ifft(c*np.conj(fft.fft(b)))
      q = q + np.absolute(r)
    idx = np.argmax(q)
    if q[idx]>m_metric:
      m_metric = q[idx]
      m_code = l2cm.code_length*(float(idx)/n)
      m_doppler = doppler
  m_code = m_code%l2cm.code_length
  return m_metric,m_code,m_doppler

#
# main program
#

# parse command-line arguments
# example:
#   ./acquire-gps-l2cm.py /dev/stdin 69984000 -127126

filename = sys.argv[1]        # input data, raw file, i/q interleaved, 8 bit signed (two's complement)
fs = float(sys.argv[2])       # sampling rate, Hz
coffset = float(sys.argv[3])  # offset to L1 carrier, Hz (positive or negative)

# read first 75 ms of file

n = int(fs*0.075)
fp = open(filename,"rb")
x = io.get_samples_complex(fp,n)

# resample to 4.096 MHz

fsr = 4096000.0/fs
nco.mix(x,-coffset/fs,0)
h = scipy.signal.firwin(161,1.5e6/(fs/2),window='hanning')
x = scipy.signal.filtfilt(h,[1],x)
xr = np.interp((1/fsr)*np.arange(75*4096),np.arange(len(x)),np.real(x))
xi = np.interp((1/fsr)*np.arange(75*4096),np.arange(len(x)),np.imag(x))
x = xr+(1j)*xi

# iterate over PRNs of interest

for prn in range(1,33):
  metric,code,doppler = search(x,prn)
  if metric>0.0:    # fixme: need a proper metric and threshold; and estimate cn0
    print 'prn %3d doppler % 7.1f metric %7.1f code_offset %6.1f' % (prn,doppler,metric,code)
