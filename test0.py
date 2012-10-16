#!/usr/bin/python
import unittest

from fft import *
from conv import *

class conv_test(unittest.TestCase):
  def setUp(self):
    import random
    self.n = 64
    self.a = list(random.uniform(-100,+100) for _ in range(64))
    self.b = list(random.uniform(-100,+100) for _ in range(64))

  def test_conv_and_conv_simple_return_the_same(self):
    c1 = conv_simple(self.a, self.b)
    c2 = conv       (self.a, self.b)
    for i in range(self.n):
      self.assertAlmostEqual(c1[i], c2[i])

class fft_test(unittest.TestCase):
  def test_fft_and_dft_return_the_same(self):
     x = list(range(8))

  def test_fft_0(self):
     x = list(range(8))
     y_ref = [28.0000, -4.0000 + 9.6569j, -4.0000 + 4.0000j, -4.0000 + 1.6569j, -4.0000, -4.0000 - 1.6569j, -4.0000 - 4.0000j, -4.0000 - 9.6569j]
     y     = fft(x)
     for i in range(8):
       self.assertAlmostEqual(y[i], y_ref[i], 4)

if __name__ == "__main__":
  unittest.main()
