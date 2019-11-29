#!/usr/bin/env python

import optparse

import numpy as np
import scipy.signal
import scipy.fftpack as fft

import gnsstools.beidou.b1cd as b1cd
import gnsstools.nco as nco
import gnsstools.io as io
import gnsstools.util as util

#
# Acquisition search
#

def search(x,prn,doppler_search,ms):
  blocks = ms//10
  fs = 8192000.0
  n = 81920                                        # 10 ms coherent integration
  doppler_min, doppler_max, doppler_incr = doppler_search
  incr = float(b1cd.code_length)/n
  c = b1cd.code(prn,0,0,incr,n)                     # obtain samples of the B1Cd code
  boc = nco.boc11(0,0,incr,n)
  c = fft.fft(c*boc)
  m_metric,m_code,m_doppler = 0,0,0
  for doppler in np.arange(doppler_min,doppler_max,doppler_incr):        # doppler bins
    q = np.zeros(n)
    w = nco.nco(-doppler/fs,0,n)
    for block in range(blocks):                    # incoherent sums
      b = x[(block*n):((block+1)*n)]
      b = b*w
      r = fft.ifft(c*np.conj(fft.fft(b)))
      q = q + np.absolute(r)
    idx = np.argmax(q)
    if q[idx]>m_metric:
      m_metric = q[idx]
      m_code = b1cd.code_length*(float(idx)/n)
      m_doppler = doppler
  m_code = m_code%b1cd.code_length
  return m_metric,m_code,m_doppler

#
# main program
#

parser = optparse.OptionParser(usage="""acquire-beidou-b1cd.py [options] input_filename sample_rate carrier_offset

Acquire BeiDou B1Cd signals

Examples:
  Acquire all BeiDou PRNs using standard input with sample rate 69.984 MHz and carrier offset -9.334875 MHz:
    acquire-beidou-b1i.py /dev/stdin 69984000 -9334875

Arguments:
  input_filename    input data file, i/q interleaved, 8 bit signed
  sample_rate       sampling rate in Hz
  carrier_offset    offset to B1C carrier in Hz (positive or negative)""")

parser.disable_interspersed_args()

parser.add_option("--prn", default="1-63", help="PRNs to search, e.g. 1,3,7-14,31 (default %default)")
parser.add_option("--doppler-search", metavar="MIN,MAX,INCR", default="-7000,7000,20", help="Doppler search grid: min,max,increment (default %default)")
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

# resample to 8.192 MHz

fsr = 8192000.0/fs
nco.mix(x,-coffset/fs,0)
h = scipy.signal.firwin(161,4e6/(fs/2),window='hanning')
x = scipy.signal.filtfilt(h,[1],x)
xr = np.interp((1/fsr)*np.arange(ms_pad*8192),np.arange(len(x)),np.real(x))
xi = np.interp((1/fsr)*np.arange(ms_pad*8192),np.arange(len(x)),np.imag(x))
x = xr+(1j)*xi

# iterate (in parallel) over PRNs of interest

def worker(p):
  x,prn = p
  metric,code,doppler = search(x,prn,doppler_search,ms)
  return 'prn %3d doppler % 7.1f metric % 7.1f code_offset %6.1f' % (prn,doppler,metric,code)

import multiprocessing as mp

cpus = mp.cpu_count()
results = mp.Pool(cpus).map(worker, map(lambda prn: (x,prn),prns))

for r in results:
  print(r)
