import fft

def transpose(a):
  return map(list, zip(*a))

def fft2(a):
  ny = len(a); 
  if ny == 0: return a
  nx = len(a[0]); 
  if nx == 0: return a
  for y in range(ny): assert len(a[y]) == nx

  b = [None] * ny
  for y in range(ny):
    b[y] = fft.fft(a[y])

  b = transpose(b)

  for x in range(nx):
    b[x] = fft.fft(b[x])

  b = transpose(b)

  return b
