#!/usr/bin/env python

import sys
import numpy as np
import scipy.fftpack as fft

import gnsstools.io as io

# ./spectrum.py /dev/stdin 1239716571.429 68873142.857 4096 1000

filename = sys.argv[1]             # input data, raw file, i/q interleaved, 8 bit signed (two's complement)
fc = float(sys.argv[2])            # center frequency, Hz
fs = float(sys.argv[3])            # sampling rate, Hz
n = int(sys.argv[4])               # FFT length
ns = int(sys.argv[5])              # number of blocks

fp = open(filename,"rb")

p = np.zeros(n)
w = np.hanning(n)

for k in range(ns):
  x = io.get_samples_complex(fp,n)
  z = fft.fft(x*w)
  p += np.real(z*np.conj(z))/ns

import pylab

f = fc + fs*((np.arange(n)-(n/2.0))/n)
pylab.plot(f/1e6,10*np.log10(fft.fftshift(p)))
pylab.xlabel('Frequency (MHz)')
pylab.ylabel('Power spectral density (dB)')
pylab.title('Spectrum')
pylab.grid('on')
pylab.show()
