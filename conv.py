from fft import fft, ifft

def conv_simple(a, b):
  assert len(a) == len(b)
  n = len(a)

  c = [0] * n
  for i in range(n):
    c[i] = sum(a[j] * b[(i-j)%n] for j in range(n))
  return c

def conv(a, b):
  assert len(a) == len(b)
  n = len(a)

  a, b = fft(a), fft(b)
  return ifft(list(a[i] * b[i] / n for i in range(len(a))))
