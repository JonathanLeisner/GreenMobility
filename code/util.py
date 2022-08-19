"""  Defines the class that runs the Green Mobility (GM) project code."""
#%%

import numpy as np
from scipy import optimize

#%%

def check_grad(x0, f, jac, args_f=None, kwargs_jac=None, epsilon=1.4901161193847656e-08):
    """ Approximates the gradient of f and compares it to the calculating gradients from jac."""
    if args_f is None:
        args_f = []
    if kwargs_jac is None:
        kwargs_jac = {}
    approx_fprime = optimize.approx_fprime(xk=x0, f=f, epsilon=epsilon, *args_f)
    err = approx_fprime - jac(x0, **kwargs_jac)
    return np.sqrt(np.sum(np.square(err)))

def c_m1(xi_in, xi_out):
    m1 = np.exp(np.sum(np.meshgrid(xi_in, xi_out), axis=0)) #add cost of going OUT of a sector with the cost of going IN to a sector.
    np.fill_diagonal(m1, 0) #Not switching sector is costless
    return m1

def c_m2(ages, kappa, scale_kappa):
    return np.exp(np.array(ages * kappa[0] / scale_kappa + np.square(ages) * kappa[1] / scale_kappa))