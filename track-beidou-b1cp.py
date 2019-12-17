#!/usr/bin/env python

import optparse

import numpy as np

import gnsstools.beidou.b1cp as b1cp
import gnsstools.nco as nco
import gnsstools.io as io
import gnsstools.discriminator as discriminator
import gnsstools.util as util

class tracking_state:
  def __init__(self,fs,prn,code_p,code_f,code_i,carrier_p,carrier_f,carrier_i,mode):
    self.fs = fs
    self.prn = prn
    self.code_p = code_p
    self.code_f = code_f
    self.code_i = code_i
    self.carrier_p = carrier_p
    self.carrier_f = carrier_f
    self.carrier_i = carrier_i
    self.mode = mode
    self.prompt1 = 0 + 0*(1j)
    self.carrier_e1 = 0
    self.code_e1 = 0
    self.eml = 0

# tracking loops

def track(x,s):
  n = len(x)
  fs = s.fs

  nco.mix(x,-s.carrier_f/fs, s.carrier_p)
  s.carrier_p = s.carrier_p - n*s.carrier_f/fs
  s.carrier_p = np.mod(s.carrier_p,1)

  cf = (s.code_f+s.carrier_f/1540.0)/fs

  p_early = b1cp.correlate(x, s.prn, 0, s.code_p-0.2, cf, b1cp.b1cp_code(prn), b1cp.boc11)
  p_prompt = b1cp.correlate(x, s.prn, 0, s.code_p, cf, b1cp.b1cp_code(prn), b1cp.boc11)
  p_late = b1cp.correlate(x, s.prn, 0, s.code_p+0.2, cf, b1cp.b1cp_code(prn), b1cp.boc11)

  if s.mode=='FLL_WIDE':
    fll_k = 3.0
    a = p_prompt
    b = s.prompt1
    e = discriminator.fll_atan(a,b)
    s.carrier_f = s.carrier_f + fll_k*e
    s.prompt1 = p_prompt
  elif s.mode=='FLL_NARROW':
    fll_k = 0.8
    a = p_prompt
    b = s.prompt1
    e = discriminator.fll_atan(a,b)
    s.carrier_f = s.carrier_f + fll_k*e
    s.prompt1 = p_prompt
  elif s.mode=='PLL':
    pll_k1 = 0.1
    pll_k2 = 3.5
    e = discriminator.pll_costas(p_prompt)
    e1 = s.carrier_e1
    s.carrier_f = s.carrier_f + pll_k1*e + pll_k2*(e-e1)
    s.carrier_e1 = e

# code loop

  dll_k1 = 0.00002
  dll_k2 = 0.2
  s.early = np.absolute(p_early)
  s.prompt = np.absolute(p_prompt)
  s.late = np.absolute(p_late)
  if (s.late+s.early)==0:
    e = 0
  else:
    e = (s.late-s.early)/(s.late+s.early)
  s.eml = e
  e1 = s.code_e1
  s.code_f = s.code_f + dll_k1*e + dll_k2*(e-e1)
  s.code_e1 = e

  s.code_p = s.code_p + n*cf
  s.code_p = np.mod(s.code_p,b1cp.code_length)

  return p_prompt,s

#
# main program
#

parser = optparse.OptionParser(usage="""track-beidou-b1cp.py [options] input_filename sample_rate carrier_offset PRN doppler code_offset

Track BeiDou B1Cp signal

Examples:
  Track with default options:
    track-beidou-b1cp.py /dev/stdin 69984000 -9334875 31 1200.0 831.15
  Track with pure PLL (no FLL intervals at the start) and with a specified carrier phase:
    track-beidou-b1cp.py --carrier-phase 0.214 /dev/stdin 69984000 -9334875 31 1200.0 831.15

Arguments:
  input_filename    input data file, i/q interleaved, 8 bit signed
  sample_rate       sampling rate in Hz
  carrier_offset    offset to L1 carrier in Hz (positive or negative)
  PRN               PRN
  doppler           Doppler estimate from acquisition
  code_offset       Code-offset estimate from acquisition""")

parser.disable_interspersed_args()

parser.add_option("--loop-dwells", default="500,500", help="initial time intervals for wide FLL, then narrow FLL, in milliseconds (default %default)")
parser.add_option("--carrier-phase", help="initial carrier phase in cycles (disables FLL: uses PLL from the start)")

(options, args) = parser.parse_args()

filename = args[0]
fs = float(args[1])
coffset = float(args[2])
prn = int(args[3])
doppler = float(args[4])
code_offset = float(args[5])
loop_dwells = util.parse_list_floats(options.loop_dwells)
carrier_p = 0.0
if options.carrier_phase is not None:
  carrier_p = float(options.carrier_phase)
  loop_dwells = 0,0

fll_wide_time,fll_narrow_time = loop_dwells

fp = open(filename,"rb")

n = int(fs*0.01*((b1cp.code_length-code_offset)/b1cp.code_length))  # align with 10 ms code boundary
x = io.get_samples_complex(fp,n)
code_offset += n*100.0*b1cp.code_length/fs

s = tracking_state(fs=fs, prn=prn,                    # initialize tracking state
  code_p=code_offset, code_f=b1cp.chip_rate, code_i=0,
  carrier_p=carrier_p, carrier_f=doppler, carrier_i=0,
  mode='FLL_WIDE')

block = 0
coffset_phase = 0.0

while True:
  if block>=fll_wide_time:
    s.mode = 'FLL_NARROW'
  if block>=fll_wide_time+fll_narrow_time:
    s.mode = 'PLL'

  if s.code_p<b1cp.code_length/2:
    n = int(fs*0.01*(b1cp.code_length-s.code_p)/b1cp.code_length)
  else:
    n = int(fs*0.01*(2*b1cp.code_length-s.code_p)/b1cp.code_length)

  x = io.get_samples_complex(fp,n)
  if x is None:
    break

  nco.mix(x,-coffset/fs,coffset_phase)
  coffset_phase = coffset_phase - n*coffset/fs
  coffset_phase = np.mod(coffset_phase,1)

  for j in range(10):
    a,b = int(j*n/10),int((j+1)*n/10)
    p_prompt,s = track(x[a:b],s)
    vars = block, np.real(p_prompt), np.imag(p_prompt), s.carrier_f, s.code_f-b1cp.chip_rate, (180/np.pi)*np.angle(p_prompt), s.early, s.prompt, s.late
    print('%d %f %f %f %f %f %f %f %f' % vars)

    block = block + 1
