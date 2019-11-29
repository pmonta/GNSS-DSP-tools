#!/usr/bin/env python

import optparse

import numpy as np
import sys

def read_block(N):
  xi = np.zeros(N)
  xq = np.zeros(N)
  for i in range(N):
    t = sys.stdin.readline()
    if not t:
      return np.array([])
    t = t.split()
    xi[i] = float(t[1])
    xq[i] = float(t[2])
  return xi+(1j)*xq

def cn0(x):
  s = np.mean(np.abs(np.real(x)))
  r = np.sqrt(2)*np.std(np.imag(x))
  snr = 20*np.log10(s/r)
  cn0 = snr + 30
  return cn0  

#
# main program
#

parser = optparse.OptionParser(usage="""cn0.py [options]

Estimate carrier-to-noise (C/N_0) from a track sampled at 1 kHz

Examples:
  Print cn0 estimates with default options:
    cn0.py <track.dat
  Use two-second blocks, rather than the default:
    cn0.py --time 2000 <track.dat""")

parser.disable_interspersed_args()

parser.add_option("--time", default="300", help="integration time in milliseconds (default %default)")

(options, args) = parser.parse_args()

N = int(options.time)

while True:
  x = read_block(N)
  if len(x)==0:
    break
  print("%.2f" % cn0(x))
