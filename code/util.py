"""  Defines the class that runs the Green Mobility (GM) project code."""
#%%

import numpy as np
from scipy import optimize

#%%

def check_grad(x0, f, jac, args_f=None, kwargs_jac=None):
    """ Approximates the gradient of f and compares it to the calculating gradients from jac."""
    if args_f is None:
        args_f = []
    if kwargs_jac is None:
        kwargs_jac = {}
    approx_fprime = optimize.approx_fprime(xk=x0, f=f, epsilon=1.4901161193847656e-08, *args_f)
    err = approx_fprime - jac(x0, **kwargs_jac)
    return np.sqrt(np.sum(np.square(err)))
