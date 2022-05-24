import numpy as np
from math import exp,sqrt
from scipy.stats import norm

bottom = -1.96
top = 1.96
ndraws=100000

print('built-in:    ',norm.cdf(top) - norm.cdf(bottom))

# draw from uniform and put it into an numpy vector of size ndraws
a = np.random.uniform(bottom,top,size=ndraws)

# compute the normal pdf of the elements of a
# to exponentiate the elements of a vector at once, I need np.exp
# (exp exponentiates only scalars)
b = 1/sqrt(2*np.pi) * np.exp(-1/2*a**2)

# note that montecarlo is less efficient than the built-in norm.cdf
# routine, with 100K draws it's still not great
print('monte-carlo: ',b.mean()*(top-bottom))