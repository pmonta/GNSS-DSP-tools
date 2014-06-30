import sys
import os
import numpy as np
import scipy.signal
import scipy.fftpack as fft

import gnsstools.gps.ca as ca
import gnsstools.nco as nco

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
      m_code = idx
      m_doppler = doppler
  return m_metric,m_code,m_doppler

#
# main program
#

# read first 25 ms of file

fref = (38880000.0*62)/70
coffset = 1575420000 - 46*fref
fs = 2*fref
n = 2*int(round((2*fs*0.025)/2))
fp = open("data/gps-5001-l1_a.dat","rb")
#fp.seek(40*n, os.SEEK_SET)
s = np.fromfile(fp,'b',n)
fp.close()
x = s[0:n:2] + (1j)*s[1:n:2]

# resample to 4.096 MHz

fsr = 4096000.0/fs
x = x * nco.nco(-coffset/fs,0,len(x))
h = scipy.signal.firwin(41,3e6/(fs/2),window='hanning')
x = scipy.signal.filtfilt(h,[1],x)
xr = np.interp((1/fsr)*np.arange(25*4096),np.arange(len(x)),np.real(x))
xi = np.interp((1/fsr)*np.arange(25*4096),np.arange(len(x)),np.imag(x))
x = xr+(1j)*xi

# iterate over PRNs of interest

for prn in range(1,33)+[133,135,138]:
  metric,code,doppler = search(x,prn)
  if metric>2200.0:
    print 'gps_ca prn %3d doppler % 7.1f metric %7.1f code_offset %6.1f' % (prn,doppler,metric,code)
