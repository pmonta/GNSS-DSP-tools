#!/usr/bin/env python

import optparse

import numpy as np
import scipy.signal
import scipy.fftpack as fft

import gnsstools.gps.ca as ca
import gnsstools.nco as nco
import gnsstools.io as io
import gnsstools.util as util

#
# Acquisition search
#

def search(x,prn,doppler_search,ms):
  fs = 4096000.0
  n = 4096                                         # 1 ms coherent integration
  doppler_min, doppler_max, doppler_incr = doppler_search
  incr = float(ca.code_length)/n
  c = ca.code(prn,0,0,incr,n)                      # obtain samples of the C/A code
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
      m_code = ca.code_length*(float(idx)/n)
      m_doppler = doppler
  return m_metric,m_code,m_doppler

#
# main program
#

parser = optparse.OptionParser(usage="""acquire-gps-l1.py [options] input_filename sample_rate carrier_offset

Acquire GPS L1 signals

Examples:
  Acquire all GPS PRNs using standard input with sample rate 69.984 MHz and carrier offset -9.334875 MHz:
    acquire-gps-l1.py /dev/stdin 69984000 -9334875
  Acquire all GPS and WAAS PRNs with a custom doppler search grid and integration time of 20 ms:
    acquire-gps-l1.py --prn 1-32,131,133,135,138 --doppler-search -6000,6000,500 --time 20 /dev/stdin 69984000 -9334875
  Acquire all GPS and QZSS PRNs from raw sample file "recording.iq":
    acquire-gps-l1.py --prn 1-32,193,194,195,199 recording.iq 69984000 -9334875

Arguments:
  input_filename    input data file, i/q interleaved, 8 bit signed
  sample_rate       sampling rate in Hz
  carrier_offset    offset to L1 carrier in Hz (positive or negative)""")

parser.disable_interspersed_args()

parser.add_option("--prn", default="1-32", help="PRNs to search, e.g. 1,3,7-14,31 (default %default)")
parser.add_option("--doppler-search", metavar="MIN,MAX,INCR", default="-7000,7000,200", help="Doppler search grid: min,max,increment (default %default)")
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

# wipe off nominal offset from channel center to GPS L1 carrier

nco.mix(x,-coffset/fs,0)

# resample to 4.096 MHz

fsr = 4096000.0/fs
h = scipy.signal.firwin(161,1.5e6/(fs/2),window='hann')
x = scipy.signal.filtfilt(h,[1],x)
xr = np.interp((1/fsr)*np.arange(ms_pad*4096),np.arange(len(x)),np.real(x))
xi = np.interp((1/fsr)*np.arange(ms_pad*4096),np.arange(len(x)),np.imag(x))
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
