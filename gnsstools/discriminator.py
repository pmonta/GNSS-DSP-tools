import numpy as np

def pll_costas(x):
  if np.real(x)>0:
    return np.arctan2(np.imag(x),np.real(x))
  else:
    return np.arctan2(-np.imag(x),-np.real(x))

def fll_atan(x,x1):
  if np.real(x)==0:
    t = np.pi/2
  else:
    t = np.arctan(np.imag(x)/np.real(x))
  if np.real(x1)==0:
    t1 = np.pi/2
  else:
    t1 = np.arctan(np.imag(x1)/np.real(x1))
  d = t - t1
  if d>np.pi/2:
    d = np.pi-d
  if d<-np.pi/2:
    d = -np.pi-d
  return d

def fll_atan2(a,b):
  d = np.arctan2(np.imag(a)*np.real(b)-np.real(a)*np.imag(b),np.real(a)*np.real(b)+np.imag(a)*np.imag(b))
  return d
