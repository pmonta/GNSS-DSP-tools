#!/usr/bin/env python

import optparse

import numpy as np
import scipy.signal
import scipy.fftpack as fft

import gnsstools.xona.x5p as x5p
import gnsstools.nco as nco
import gnsstools.io as io
import gnsstools.util as util

#
# Acquisition search
#

def search(x,prn,doppler_search,ms):
  fs = 3*10230000.0
  n = 3*10230                                       # 3 ms coherent integration
  doppler_min, doppler_max, doppler_incr = doppler_search
  incr = float(x5p.code_length)/n
  c = x5p.code(prn,0,0,incr,n)                      # obtain samples of the X5P code
  c = fft.fft(c)
  m_metric,m_code,m_doppler = 0,0,0
  for doppler in np.arange(doppler_min,doppler_max,doppler_incr):        # doppler bins
    q = np.zeros(n)
    w = nco.nco(-doppler/fs,0,n)
    for block in range(ms):                        # incoherent sums
      b = x[(block*n):((block+1)*n)]
      b = b*w
      r = fft.ifft(c*np.conj(fft.fft(b)))
      q = q + np.absolute(r)
    idx = np.argmax(q)
    metric = q[idx]/np.mean(q)
    if metric>m_metric:
      m_metric = metric
      m_code = x5p.code_length*(float(idx)/n)
      m_doppler = doppler
  return m_metric,m_code,m_doppler

#
# main program
#

parser = optparse.OptionParser(usage="""acquire-xona-x5p.py [options] input_filename sample_rate carrier_offset

Acquire Xona X5 pilot signals

Examples:
  Acquire Xona using standard input with sample rate 69.984 MHz and carrier offset -1.125375 MHz:
    acquire-xona-x5p.py /dev/stdin 69984000 -1125375
  Acquire Xona with a custom doppler search grid and integration time of 20 ms:
    acquire-xona-x5p.py --doppler-search -50000,50000,100 --time 20 /dev/stdin 69984000 -1125375
  Acquire Xona from raw sample file "recording.iq":
    acquire-xona-x5p.py recording.iq 69984000 -1125375

Arguments:
  input_filename    input data file, i/q interleaved, 8 bit signed
  sample_rate       sampling rate in Hz
  carrier_offset    offset to X5 carrier in Hz (positive or negative)""")

parser.disable_interspersed_args()

parser.add_option("--prn", default="0", help="PRNs to search, e.g. 1,3,7-14,31 (default %default)")
parser.add_option("--doppler-search", metavar="MIN,MAX,INCR", default="-50000,50000,200", help="Doppler search grid: min,max,increment (default %default)")
parser.add_option("--time", type="int", default=80, help="integration time in milliseconds (default %default)")

(options, args) = parser.parse_args()

filename = args[0]
fs = float(args[1])
coffset = float(args[2])
prns = util.parse_list_ranges(options.prn)
doppler_search = util.parse_list_floats(options.doppler_search)
ms = options.time

# read first portion of file

ms_pad = ms + 5
n = int(fs*0.001*ms_pad)
fp = open(filename,"rb")
x = io.get_samples_complex(fp,n)

# wipe off nominal offset from channel center to Xona X5 carrier

nco.mix(x,-coffset/fs,0)

# resample to 3*10.230 MHz

fsr = 3*10230000.0/fs
h = scipy.signal.firwin(161,12e6/(fs/2),window='hann')
x = scipy.signal.filtfilt(h,[1],x)
xr = np.interp((1/fsr)*np.arange(ms_pad*3*10230),np.arange(len(x)),np.real(x))
xi = np.interp((1/fsr)*np.arange(ms_pad*3*10230),np.arange(len(x)),np.imag(x))
x = xr+(1j)*xi

# iterate (in parallel) over PRNs of interest

def worker(p):
  x,prn = p
  metric,code,doppler = search(x,prn,doppler_search,ms)
  return 'prn %3d doppler % 7.1f metric % 5.2f code_offset %6.1f' % (prn,doppler,metric,code)

import multiprocessing as mp

cpus = mp.cpu_count()
results = mp.Pool(cpus).map(worker, map(lambda prn: (x,prn),prns))

for r in results:
  print(r)
