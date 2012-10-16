def bit_reverse(x):
  n = len(x)

  j = 0
  for i in range(0, n-1):
    # count i up in normal order

    # do the swapping
    if j > i: x[i], x[j] = x[j], x[i]
    
    # count j up in bit-reversed order
    k = n >> 1
    while k <= j:
    	j -= k
    	k >>= 1
    j += k

  return x
