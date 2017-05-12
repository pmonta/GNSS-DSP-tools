#!/usr/bin/env python

import sys
import numpy as np
import scipy.fftpack as fft
import matplotlib.pyplot as plt

import gnsstools.io as io

class myplot:
  def __init__(self,fc,fs,n):
    plt.ion()
    self.fc = fc
    self.fs = fs
    self.n = n
    self.fig = plt.figure()
    self.ax = self.fig.add_subplot(111)
    self.x = (fc + fs*((np.arange(n)-(n/2.0))/n))/1e6
    self.y = np.zeros(self.n)
    self.line, = self.ax.plot(self.x,self.y)
    self.ax.relim()
    self.ax.autoscale_view(True,True,True)
#    self.ax.axis([1255,1260,20,60])
    self.ax.set_xlabel('Frequency (MHz)')
    self.ax.set_ylabel('Power spectral density (dB)')
    self.ax.set_title('Spectrum')
    self.ax.grid('on')
    self.fig.canvas.draw()
  def update(self,new_y):
    self.line.set_ydata(new_y)
    self.ax.relim()
    self.ax.autoscale_view(True,True,True)
    self.fig.canvas.draw()
    plt.pause(0.1)

# ./spectrum.py /dev/stdin 1584754875 69984000 2048 1000
# ./spectrum.py /dev/stdin 1227727126 69984000 2048 1000

filename = sys.argv[1]             # input data, raw file, i/q interleaved, 8 bit signed (two's complement)
fc = float(sys.argv[2])            # center frequency, Hz
fs = float(sys.argv[3])            # sampling rate, Hz
n = int(sys.argv[4])               # FFT length
ns = int(sys.argv[5])              # number of blocks

fp = open(filename,"rb")
m = myplot(fc,fs,n)

while True:
  p = np.zeros(n)
  w = np.hanning(n)
  for k in range(ns):
    x = io.get_samples_complex(fp,n)
    if x is None:
      sys.exit()
    z = fft.fft(x*w)
    p += np.real(z*np.conj(z))/ns
  m.update(10*np.log10(fft.fftshift(p)))
