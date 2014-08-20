#!/usr/bin/env python

import sys
import os
import numpy as np
import scipy.signal
import scipy.fftpack as fft

import gnsstools.glonass.l3q as l3q
import gnsstools.nco as nco
import gnsstools.io as io

#
# Acquisition search
#

def search(x,prn):
  fs = 3*10230000.0
  n = 3*10230                                       # 1 ms coherent integration
  incr = float(l3q.code_length)/n
  c = l3q.code(prn,0,0,incr,n)                     # obtain samples of the L3-I code
  c = fft.fft(np.concatenate((c,np.zeros(n))))
  m_metric,m_code,m_doppler = 0,0,0
  for doppler in np.arange(-4000,4000,200):         # doppler bins
    q = np.zeros(2*n)
    w = nco.nco(-doppler/fs,0,2*n)
    for block in range(80):                         # 80 incoherent sums
      b = x[(block*n):((block+2)*n)]
      b = b*w
      r = fft.ifft(c*np.conj(fft.fft(b)))
      q = q + np.absolute(r)
    idx = np.argmax(q)
    if q[idx]>m_metric:
      m_metric = q[idx]
      m_code = l3q.code_length*(float(idx)/n)
      m_doppler = doppler
  m_code = m_code%l3q.code_length
  return m_metric,m_code,m_doppler

#
# main program
#

# parse command-line arguments
# example:
#   ./acquire-glonass-l3q.py /dev/stdin 68873142.857 -3255000.000

filename = sys.argv[1]        # input data, raw file, i/q interleaved, 8 bit signed (two's complement)
fs = float(sys.argv[2])       # sampling rate, Hz
coffset = float(sys.argv[3])  # offset to L3 GLONASS carrier, Hz (positive or negative)

# read first 85 ms of file

n = int(fs*0.085)
fp = open(filename,"rb")
x = io.get_samples_complex(fp,n)

# resample to 3*10.230 MHz

fsr = 3*10230000.0/fs
nco.mix(x,-coffset/fs,0,nco.nco_table)
h = scipy.signal.firwin(161,12e6/(fs/2),window='hanning')
x = scipy.signal.filtfilt(h,[1],x)
xr = np.interp((1/fsr)*np.arange(85*3*10230),np.arange(len(x)),np.real(x))
xi = np.interp((1/fsr)*np.arange(85*3*10230),np.arange(len(x)),np.imag(x))
x = xr+(1j)*xi

# iterate over channels of interest

#for prn in range(1,33):
for prn in [28,29,30]:
  metric,code,doppler = search(x,prn)
  if metric>0.0:    # fixme: need a proper metric and threshold; and estimate cn0
    print 'prn %2d doppler % 7.1f metric %7.1f code_offset %6.1f' % (prn,doppler,metric,code)
