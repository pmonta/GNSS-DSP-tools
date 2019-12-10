#!/usr/bin/env python

import optparse

import numpy as np
import scipy.signal
import scipy.fftpack as fft

import gnsstools.glonass.ca as ca
import gnsstools.nco as nco
import gnsstools.io as io
import gnsstools.util as util

#
# Acquisition search
#

def search(x,chan,doppler_search,ms):
  fs = 16384000.0
  n = 16384
  doppler_min, doppler_max, doppler_incr = doppler_search
  incr = float(ca.code_length)/n
  c = ca.code(0,0,incr,n)                          # obtain samples of the C/A code
  c = fft.fft(c)
  m_metric,m_code,m_doppler = 0,0,0
  for doppler in np.arange(doppler_min,doppler_max,doppler_incr):        # doppler bins
    q = np.zeros(n)
    w = nco.nco(-(562500*chan+doppler)/fs,0,n)
    for block in range(ms):                        # incoherent sums
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

parser = optparse.OptionParser(usage="""acquire-glonass-l1.py [options] input_filename sample_rate carrier_offset

Acquire GLONASS L1 signals

Examples:
  Acquire all GLONASS channels using standard input with sample rate 69.984 MHz and carrier (channel 0) offset 17.245125 MHz:
    acquire-glonass-l1.py /dev/stdin 69984000 17245125

Arguments:
  input_filename    input data file, i/q interleaved, 8 bit signed
  sample_rate       sampling rate in Hz
  carrier_offset    offset to GLONASS L1 carrier (channel 0) in Hz (positive or negative)""")

parser.disable_interspersed_args()

parser.add_option("--channel", default="-7:7", help="channels to search, e.g. -6,-4,-1:2,7 (default %default)")
parser.add_option("--doppler-search", metavar="MIN,MAX,INCR", default="-7000,7000,200", help="Doppler search grid: min,max,increment (default %default)")
parser.add_option("--time", type="int", default=80, help="integration time in milliseconds (default %default)")

(options, args) = parser.parse_args()

filename = args[0]
fs = float(args[1])
coffset = float(args[2])
chans = util.parse_list_ranges(options.channel,sep=':')
doppler_search = util.parse_list_floats(options.doppler_search)
ms = options.time

# read first portion of file

ms_pad = ms + 5
n = int(fs*0.001*ms_pad)
fp = open(filename,"rb")
x = io.get_samples_complex(fp,n)

# wipe off nominal offset from channel center to GLONASS L1 carrier

nco.mix(x,-coffset/fs,0)

# resample to 16.384 MHz

fsr = 16384000.0/fs
h = scipy.signal.firwin(161,6e6/(fs/2),window='hanning')
x = scipy.signal.filtfilt(h,[1],x)
xr = np.interp((1/fsr)*np.arange(ms_pad*16384),np.arange(len(x)),np.real(x))
xi = np.interp((1/fsr)*np.arange(ms_pad*16384),np.arange(len(x)),np.imag(x))
x = xr+(1j)*xi

# iterate (in parallel) over channels of interest

def worker(p):
  x,chan = p
  metric,code,doppler = search(x,chan,doppler_search,ms)
  return 'chan % 2d doppler % 7.1f metric % 7.1f code_offset %7.2f' % (chan,doppler,metric,code)

import multiprocessing as mp

cpus = mp.cpu_count()
results = mp.Pool(cpus).map(worker, map(lambda chan: (x,chan),chans))

for r in results:
  print(r)
