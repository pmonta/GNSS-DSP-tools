#!/usr/bin/env python

import sys
import numpy as np

import gnsstools.nco as nco
import gnsstools.io as io
import gnsstools.squaring as squaring

# parse command-line arguments
# example:
#   ./squaring.py /dev/stdin 69984000 -9334875 | baudline -reset -stdin -fftsize 16384 -samplerate 43740 -channels 2 -format le16 -quadrature -flipcomplex

filename = sys.argv[1]             # input data, raw file, i/q interleaved, 8 bit signed (two's complement)
fs = float(sys.argv[2])            # sampling rate, Hz
coffset = float(sys.argv[3])       # offset to L1 carrier, Hz (positive or negative)

fp = open(filename,"rb")

coffset_phase = 0.0

b = 1000
n = 16
m = 100
r = np.zeros(b).astype('complex')
y = np.zeros(2*b).astype('int16')

while True:
  x = io.get_samples_complex(fp,b*n*m)
  if x is None:
    break

  nco.mix(x,-coffset/fs,coffset_phase)
  coffset_phase = coffset_phase - len(x)*coffset/fs
  coffset_phase = np.mod(coffset_phase,1)

  squaring.squaring(x,r,n,m)

  y[0:2*b:2] = np.round(20*np.real(r)).astype('int16')
  y[1:2*b:2] = np.round(20*np.imag(r)).astype('int16')

  y.tofile(sys.stdout)
