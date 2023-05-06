"""  Defines the class that runs the Green Mobility (GM) project code."""
#%%

import numpy as np
from scipy import optimize

#%%

def check_grad(x0, f, jac, args_f=None, kwargs_jac=None, epsilon=1.4901161193847656e-08, relative=False):
    """ Approximates the gradient of f and compares it to the calculating gradients from jac."""
    if args_f is None:
        args_f = []
    if kwargs_jac is None:
        kwargs_jac = {}
    approx_fprime = optimize.approx_fprime(xk=x0, f=f, epsilon=epsilon, *args_f)
    if not relative:
        err = approx_fprime - jac(x0, **kwargs_jac)
        return (err, np.sqrt(np.sum(np.square(err))))
    else: 
        err = (np.abs(approx_fprime) / np.abs(jac(x0, **kwargs_jac))) - 1
        return (err, np.max(np.abs(err)))
        
def update_attrs(ns, **kwargs):
    """ Helper function for updating the attributes of a namespace (ns) according to the named keyword arguments. """
    for k, v in kwargs.items():
        setattr(ns, k, v)

def get_key_from_value(dictionary, value):
    """ Helper function for finding the key of specific value in a dictionary. 
        The value must only exist under one key, otherwise an error is thrown.
    """
    keys = [k for k in dictionary.keys() if value in dictionary[k]]
    assert len(keys) == 1
    return keys[0]

def shuffle_vector(size):
    return (np.random.randint(20, 50, size) * np.random.choice([-1, 1], size)) / 100

def suffix_to_str(string, true_or_suffix):
    if true_or_suffix is None:
        return None
    elif true_or_suffix is True:
        return string
    elif isinstance(true_or_suffix, str):
        return string + "_" + true_or_suffix
