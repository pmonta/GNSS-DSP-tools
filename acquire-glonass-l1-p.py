#!/usr/bin/env python

import optparse

import numpy as np

import gnsstools.glonass.p as p
import gnsstools.nco as nco
import gnsstools.io as io

#
# Acquisition search
#

def search(x,chan,doppler,ca_code_phase,ms):
  blocks = ms//4
  n = int(fs*0.004)
  w = nco.nco(-(562500*chan+doppler)/fs,0,n)
  m_metric,m_k = 0,0
  for k in range(1000):
    q = 0
    cp = 5110*k + 10*ca_code_phase
    for block in range(blocks):
      incr = 5110000.0/fs
      c = p.code(0,cp,incr,n)
      xp = x[n*block:n*(block+1)]*c*w
      q = q + np.absolute(np.sum(xp))
      cp += n*incr
#    print('%f %f'%(k,q))
    if q>m_metric:
      m_metric = q
      m_k = k
  return m_metric,m_k

#
# main program
#

parser = optparse.OptionParser(usage="""acquire-glonass-l1-p.py [options] input_filename sample_rate carrier_offset

Acquire GLONASS L1-P signals (given an L1-CA acquisition)

Examples:
  Acquire a GLONASS L1-P signal using standard input with sample rate 69.984 MHz, carrier (channel 0) offset 17.245125 MHz,
  RF channel -4, doppler 2600 Hz, and CA code phase 278.6 chips:
    acquire-glonass-l1-p.py /dev/stdin 69984000 17245125 -4 2600.0 278.6

Arguments:
  input_filename    input data file, i/q interleaved, 8 bit signed
  sample_rate       sampling rate in Hz
  carrier_offset    offset to GLONASS L1 carrier (channel 0) in Hz (positive or negative)
  channel           RF channel index
  doppler           Doppler from CA acquisition
  ca_code_phase     Code phase from CA acquisition""")

parser.disable_interspersed_args()

parser.add_option("--time", type="int", default=80, help="integration time in milliseconds (default %default)")

(options, args) = parser.parse_args()

filename = args[0]
fs = float(args[1])
coffset = float(args[2])
chan = int(args[3])
doppler = float(args[4])
ca_code_phase = float(args[5])

ms = options.time

# read first portion of file

ms_pad = ms + 5
n = int(fs*0.001*ms_pad)
fp = open(filename,"rb")
x = io.get_samples_complex(fp,n)

nco.mix(x,-coffset/fs,0)

metric,k = search(x,chan,doppler,ca_code_phase,ms)
print('%f %f'%(5110*k+10*ca_code_phase,metric))
