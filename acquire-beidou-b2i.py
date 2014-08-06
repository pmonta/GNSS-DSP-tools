#!/usr/bin/env python

import sys
import os
import numpy as np
import scipy.signal
import scipy.fftpack as fft

import gnsstools.beidou.b1i as b1i
import gnsstools.nco as nco

#
# Acquisition search
#

def search(x,prn):
  fs = 4092000.0
  n = 4092                                         # 1 ms coherent integration
  incr = float(b1i.code_length)/n
  c = b1i.code(prn,0,0,incr,n)                     # obtain samples of the B1I code
  c = fft.fft(c)
  m_metric,m_code,m_doppler = 0,0,0
  for doppler in np.arange(-5000,5000,200):        # doppler bins
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
      m_code = b1i.code_length*(float(idx)/n)
      m_doppler = doppler
  return m_metric,m_code,m_doppler

#
# main program
#

# parse command-line arguments
# example:
#   ./acquire-beidou-b2i.py /dev/stdin 68873142.857 -32576571.429

filename = sys.argv[1]        # input data, raw file, i/q interleaved, 8 bit signed (two's complement)
fs = float(sys.argv[2])       # sampling rate, Hz
coffset = float(sys.argv[3])  # offset to L1 BeiDou carrier, Hz (positive or negative)

# read first 25 ms of file

n = 2*int(round((2*fs*0.025)/2))
fp = open(filename,"rb")
#fp.seek(bytes, os.SEEK_SET)  # fixme: make command-line parameter
s = np.fromfile(fp,'b',n)
fp.close()
x = s[0:n:2] + (1j)*s[1:n:2]

# resample to 4.092 MHz

fsr = 4092000.0/fs
x = x * nco.nco(-coffset/fs,0,len(x))
h = scipy.signal.firwin(81,3e6/(fs/2),window='hanning')
x = scipy.signal.filtfilt(h,[1],x)
xr = np.interp((1/fsr)*np.arange(25*4092),np.arange(len(x)),np.real(x))
xi = np.interp((1/fsr)*np.arange(25*4092),np.arange(len(x)),np.imag(x))
x = xr+(1j)*xi

# iterate over channels of interest

for prn in range(1,38):
  metric,code,doppler = search(x,prn)
  if metric>0.0:    # fixme: need a proper metric and threshold; and estimate cn0
    print 'prn %2d doppler % 7.1f metric %7.1f code_offset %6.1f' % (prn,doppler,metric,code)
