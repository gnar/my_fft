#!/usr/bin/python
import unittest

from fft import *
from fftn import fft2
from conv import *

def rnd_data(n=64):
  import random
  return list(random.uniform(-100,+100) for _ in range(n))

class conv_test(unittest.TestCase):
  def setUp(self):
    self.n = 64
    self.a = rnd_data(64)
    self.b = rnd_data(64)

  def test_conv_and_conv_simple_return_the_same(self):
    c1 = conv_simple(self.a, self.b)
    c2 = conv       (self.a, self.b)
    for i in range(self.n):
      self.assertAlmostEqual(c1[i], c2[i])

class fft_test(unittest.TestCase):

  def test_fft_and_dft_return_the_same(self):
    n = 64
    x = rnd_data(n)
    y0, y1 = dft(x), fft(x)
    for i in range(n): self.assertAlmostEqual(y0[i], y1[i], 5)

  ###################################################################
    
  x_ref = [0, 1, 2, 3, 4, 5, 6, 7]
  y_ref = [28.0000, -4.0000 + 9.6569j, -4.0000 + 4.0000j, -4.0000 + 1.6569j, -4.0000, -4.0000 - 1.6569j, -4.0000 - 4.0000j, -4.0000 - 9.6569j]

  def test_fft_0(self):
     y = fft(fft_test.x_ref)
     for i in range(8): self.assertAlmostEqual(y[i], fft_test.y_ref[i], 4)

  def test_dft_0(self):
     y = dft(fft_test.x_ref)
     for i in range(8): self.assertAlmostEqual(y[i], fft_test.y_ref[i], 4)

  a_ref = [[1,2,3,4],[5,6,7,8],[9,10,11,12],[13,14,15,16]]
  b_ref = [[ 136.00          , -008.00 + 008.00j, -008.00,           -008.00 - 008.00j],
           [-032.00 + 032.00j,       0          ,       0,                 0          ],
           [-032.00          ,       0          ,       0,                 0          ],
           [-032.00 - 032.00j,       0          ,       0,                 0          ]]

  def test_fft2(self):
    my_b = fft2(fft_test.a_ref)
    for x in range(4):
      for y in range(4): 
        self.assertAlmostEqual(my_b[y][x], fft_test.b_ref[y][x], 7)

if __name__ == "__main__":
  unittest.main()
