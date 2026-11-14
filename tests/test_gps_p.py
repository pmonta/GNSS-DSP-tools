import unittest

import numpy as np

from gnsstools.gps import p


class PCodeSamplingTest(unittest.TestCase):
  def test_one_sample_per_chip(self):
    expected = 1.0 - 2.0*p.p_code(1,0,12)
    np.testing.assert_array_equal(p.code(1,0,0.0,1.0,12),expected)

  def test_two_samples_per_chip_with_offset(self):
    expected = np.repeat(1.0 - 2.0*p.p_code(7,100,6),2)[1:]
    np.testing.assert_array_equal(p.code(7,100,0.5,0.5,11),expected)

  def test_zero_samples(self):
    self.assertEqual(p.code(1,0,0.0,1.0,0).shape,(0,))


if __name__=='__main__':
  unittest.main()
