import numpy as np 
import matplotlib.pyplot as plt
import random



sumsq = 1.1

while (sumsq >= 1.0):
	u = np.random.uniform(-1,1,1000)
	s = np.random.uniform(-1,1,1000)
	x1 = random.choice(u)
	x2 = random.choice(s)
	sumsq = x1**2 + x2**2

x = 2.0 * x1 * np.sqrt(1.0 - sumsq)
y = 2.0 * x2 * np.sqrt(1.0 - sumsq)
z = 1.0 - 2.0*(sumsq)
r2 = x**2 + y**2 + z**2

print x, y, z, r2