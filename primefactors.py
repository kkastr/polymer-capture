import numpy as np
import sys

#brute forcing is fine since number of monomers is small

def prime_factors(n):
    i = 2
    factors = []
    while i * i <= n:
        if n % i:
            i += 1
        else:
            n //= i
            factors.append(i)
    if n > 1:
        factors.append(n)
    return factors

print np.sum(prime_factors(int(sys.argv[1]))) + np.max(prime_factors(int(sys.argv[1]))) + 5