#!/usr/bin/env python

import sys
import os
import numpy as np
import scipy.signal
import scipy.fftpack as fft

import gnsstools.gps.ca as ca
import gnsstools.nco as nco
import gnsstools.io as io

#
# Acquisition search
#

def search(x,prn):
  fs = 4096000.0
  n = 4096                                         # 1 ms coherent integration
  incr = float(ca.code_length)/n
  c = ca.code(prn,0,0,incr,n)                      # obtain samples of the C/A code
  c = fft.fft(c)
  m_metric,m_code,m_doppler = 0,0,0
  for doppler in np.arange(-7000,7000,200):        # doppler bins
    q = np.zeros(n)
    w = nco.nco(-doppler/fs,0,n)
    for block in range(80):                        # 80 incoherent sums
      b = x[(block*n):((block+1)*n)]
      b = b*w
      r = fft.ifft(c*np.conj(fft.fft(b)))
      q = q + np.absolute(r)
    idx = np.argmax(q)
    metric = q[idx]/np.mean(q)
    if metric>m_metric:
      m_metric = metric
      m_code = ca.code_length*(float(idx)/n)
      m_doppler = doppler
  return m_metric,m_code,m_doppler

#
# main program
#

filename = sys.argv[1]        # input data, raw file, i/q interleaved, 8 bit signed (two's complement)
fs = float(sys.argv[2])       # sampling rate, Hz
coffset = float(sys.argv[3])  # offset to L1 carrier, Hz (positive or negative)

# read first 85 ms of file

ms = 85

n = int(fs*0.001*ms)
fp = open(filename,"rb")
x = io.get_samples_complex(fp,n)

# wipe off nominal offset from channel center to GPS L1 carrier

nco.mix(x,-coffset/fs,0)

# resample to 4.096 MHz

fsr = 4096000.0/fs
h = scipy.signal.firwin(161,1.5e6/(fs/2),window='hanning')
x = scipy.signal.filtfilt(h,[1],x)
xr = np.interp((1/fsr)*np.arange(ms*4096),np.arange(len(x)),np.real(x))
xi = np.interp((1/fsr)*np.arange(ms*4096),np.arange(len(x)),np.imag(x))
x = xr+(1j)*xi

# iterate (in parallel) over PRNs of interest

def worker(p):
  x,prn = p
  metric,code,doppler = search(x,prn)
  return 'prn %3d doppler % 7.1f metric % 5.2f code_offset %6.1f' % (prn,doppler,metric,code)

import multiprocessing as mp

prns = list(range(1,33))+[131,133,135,138]
cpus = mp.cpu_count()
results = mp.Pool(cpus).map(worker, map(lambda prn: (x,prn),prns))

for r in results:
  print(r)
