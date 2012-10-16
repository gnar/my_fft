from bitrev import bit_reverse
from cmath import exp
from math import pi, sqrt

def omega(n, k):
  return exp(2*pi*1j*k/n)

def dft(x):
  n = len(x)
  y = [0]*n
  for i in range(n):
    y[i] = sum(x[j] * omega(n, -i*j) for j in range(n))
  return y

def fft_r(x):
  n = len(x)
  if n == 1: 
    return x  
  else:
    y_e = fft_r(list(x[i*2+0] for i in range(n/2)))
    y_o = fft_r(list(x[i*2+1] for i in range(n/2)))
  
    y = [0] * n
    for i in range(n/2):
      y[i    ] = y_e[i] + y_o[i] * omega(n, -i)
      y[i+n/2] = y_e[i] - y_o[i] * omega(n, -i)
    return y

def fft(x, sign=+1):
  n = len(x)
  x = bit_reverse(x)

  fft_num = n/2
  fft_len = 2

  while fft_num > 0:
    for m in range(fft_len/2):
      W = omega(fft_len, sign * m)
      for i in range(m, n, fft_len):
        j = i + fft_len/2;
        a =     x[i];
        b = W * x[j];
        x[i] = a + b;
        x[j] = a - b;
    fft_num >>= 1
    fft_len <<= 1

  return x
  #return list(x[i] / sqrt(n) for i in range(n))

x = list(range(4))
print dft(x)
print fft(x)
