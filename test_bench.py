#!/usr/bin/python
from fft import *
from conv import *
from test import rnd_data

import time

n = 2**11

t0 = time.time()
conv(rnd_data(n), rnd_data(n))
t1 = time.time()
print t1-t0

t0 = time.time()
conv_simple(rnd_data(n), rnd_data(n))
t1 = time.time()
print t1-t0
