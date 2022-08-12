#"""  Defines the class that calculates gradients for the Green Mobility (GM) project code."""
#%%

# from operator import index
# from turtle import update
# from types import SimpleNamespace
# from collections import OrderedDict
import numpy as np
# import pandas as pd
# import matplotlib.pyplot as plt
# import seaborn as sns
# from scipy import optimize
# from numpy.random import default_rng
# from mpl_toolkits import mplot3d

# np.set_printoptions(precision=3)
# sns.set_theme()
# import warnings
# warnings.filterwarnings("error", category=RuntimeWarning)
# import pyfuncs as pyf

# def update_attrs(ns, **kwargs):
#     for k, v in kwargs.items():
#         setattr(ns, k, v)

#%%

class gradients:
    """ Class to calculate gradients for the model presented in Leisner (2023). 
        Generally, the notation c_[name] for methods means construct/create/calculate an object or attribute called [name].
    """
    # def __init__(self, partial_model, par, sol, est):
    def __init__(self, par, sol, est, partial_model):
        self.par = par
        self.sol = sol
        self.est = est
        self.partial_model = partial_model

        # Helper object for calculating gradients of excess labor demand.
        self.ts_indexes = np.transpose([np.tile(np.arange(0, par.T), par.S), np.repeat(np.arange(0, par.S), par.T)])

    def dP(self, dv, g):
        """ Calculates dP with respect to either r or parameters, controlled by the input variable g. 
            dv must also be provided.
        """
        #unpack
        P = self.sol.P
        par = self.par

        ndim_add = {"r":(np.newaxis, np.newaxis), "utility":(np.newaxis,)}.get(g)
        
        dP = P[(...,) + ndim_add]/par.sigma * (dv -  np.sum(P[(...,) + ndim_add] * dv, axis=3)[:, :, :, np.newaxis, ...])
        #this functional form might change with sigma!
        return dP


    def dHsup(self, dP, g):
        """ Calculate the partial derivative of human capital supply with respect to either skill prices or theta. 
            The input g controls this. 
            Expressions for D and dP must be supplied as well. Implicitly, this function assumes that D is constant and so D' can be ignored.
        """

        ndim_add = {"r":(np.newaxis, np.newaxis), "utility":(np.newaxis,)}.get(g)
        dHsup =  self.par.MASS[(slice(None), np.newaxis) + ndim_add] * \
                        np.sum(self.est.D[(..., np.newaxis) + ndim_add] * dP, axis=tuple(np.arange(0, len(self.sol.EV.shape) - 1)))
        return dHsup

        # def partial_ED(self, theta_or_r):
        #     """ Calculate the partial derivative of the excess labor demand functions (S x T of these) with respect to either theta or r.
        #         This partial derivative is used to find the equilibrium on the labor market but also to calculate the derivative dll/dtheta 
        #         when the model is not partial.
        #     """
        #     #unpack
        #     D = self.est.D
        #     par = self.par

        #     if theta_or_r == "theta":
        #         assert len(self.est.gs) == 1
        #         if "utility" in self.est.gs:
        #             du = self.du("utility")
        #             dv = self.dv_and_dEV(du, "utility")
        #             dP = self.dP(dv, "utility")
        #             dHsup = self.dHsup(D, dP, "utility")
        #             # dP =  P[..., np.newaxis]/par.sigma * (dv -  np.sum(P[..., np.newaxis] * dv, axis=3)[:, :, :, np.newaxis, ...])
        #             # dHsup =  par.MASS[(slice(None), np.newaxis) + (np.newaxis,)] * \
        #             #         np.sum(D[(..., np.newaxis) + (np.newaxis,)] * dP, axis=tuple(np.arange(0, len(self.sol.EV.shape) - 1)))
        #             dED = - dHsup #parameters cannot affect demand in the model so I simply leave it out

        #     elif theta_or_r == "r":
                
        #         ts_indexes = self.ts_indexes
                
        #         du = self.du("r")
        #         dv = self.dv_and_dEV(du, "r")
        #         dP = self.dP(dv, "r")
        #         dHsup = self.dHsup(D, dP, "r")

        #         # dP =  P[..., np.newaxis, np.newaxis]/par.sigma * (dv -  np.sum(P[..., np.newaxis, np.newaxis] * dv, axis=3)[:, :, :, np.newaxis, ...])
        #         # dHsup =  par.MASS[(slice(None), np.newaxis) + (np.newaxis,)*2] * \
        #         #          np.sum(D[(..., np.newaxis) + (np.newaxis,)*2] * dP, axis=tuple(np.arange(0, len(self.sol.EV.shape) - 1)))

        #         dHdem = np.zeros((par.T, par.S, par.T, par.S))
        #         dHdem[ts_indexes[:, 0], ts_indexes[:, 1], ts_indexes[:, 0], ts_indexes[:, 1]] = np.reshape(par.alpha1 * par.pY / (- np.square(par.r)), (par.T*par.S), order="F")

        #         dED = dHdem - dHsup
        #     return dED

    def du(self, g):
        """ Calculates du with respect to some group of parameters, e.g. utility parameters, sigma or r. This choice determines the functional form
            implicitly (or if you look in the code actually, explicitly) used to calculate du. """

        #I make it such that du always has the same dimension for a given parameter type (r or actual parameter).
        # assert g in ["utility", "sigma", "price"]. g refers to parameter type groups.

        #unpack
        par = self.par
        
        if g == "utility":
            #unpack
            pte = self.est.pars_to_estimate
            theta_idx = self.est.theta_idx
            w = self.sol.w

            du = np.zeros(self.sol.P.shape + (self.est.n_params,)) #wrong when sigma is also estimated I think.
            for para in [k for k in pte.keys() if k in par.groups.utility]: #could this be precomputed? Yes, in setup_estimation perhaps. 
                if para == "beta0":
                    du[..., pte[para], theta_idx[para]] = w[..., pte[para]] * par.ages[np.newaxis, :, np.newaxis, np.newaxis] / getattr(par.scale, para)
                else:
                    raise Exception()
        elif g == "r":
            #unpack
            H = self.sol.H
            ts_indexes = self.ts_indexes #Object for picking combinations where st == s't'.

            #Allocate
            du = np.zeros((par.S, par.N_ages, par.T, par.S, par.T, par.S)) #does not depend on slag when we diff wrt. r
                    
            # H has no t-dimension (given a), and is repeated for each t.
            du[:, :, ts_indexes[:, 0], ts_indexes[:, 1], ts_indexes[:, 0], ts_indexes[:, 1]] = np.repeat(H, par.T, axis=1)[np.newaxis, ...]  
        #note: dimensions of du vary depending on g
        return du

        # def dlnP(self, theta_or_r):
        #     """ Calculates dlnP with respect to either the parameter vector (theta) or skill prices (r). The choice of this determined by the argument
        #         theta_or_r. The function calculates du, and from this dv, and then finally dlnP from dv. 
        #     """
        #     assert theta_or_r in ["theta", "r"]
        #     if theta_or_r == "theta":
        #         assert len(self.est.gs) == 1
        #         if "utility" in self.est.gs:
        #             #Assumes only utility for now
        #             du = self.du("utility")
        #             dv = self.dv_and_dEV(du, "utility") #så længe det er utility-parametre, kan denne køre på samme måde med alle parametrene (ligning 44)
        #             dlnP = 1/self.par.sigma * (dv - np.sum(self.sol.P[..., np.newaxis] * dv, axis=-2)[:, :, :, np.newaxis, :])
                    
        #     if theta_or_r == "r":
        #         du = self.du("r")
        #         dv = self.dv_and_dEV(du, "r")
        #         dlnP = 1/self.par.sigma * (dv - np.sum(self.sol.P[..., np.newaxis, np.newaxis] * dv, axis=-3)[:, :, :, np.newaxis, :, :])
                
        #     return dlnP

    def dv_and_dEV(self, du, g):
        """ Calculates dv and dEV from du. g specifies whether the derivative is with respect to r or parameters, since this determines whether we have
            to add two dimensions or one (since r is defined over s and t while parameters are collected in one dimension). 
            The function only returns dv since dEV is never necessary to know on its own if we know dv.
            The calculation of derivatives follows the same backwards recursion structure as the solution to the worker's problem.
        """

        #unpack
        par = self.par
        P = self.sol.P

        
        if g is "r":
            ndim_add = (np.newaxis, np.newaxis) #when differentiating wrt. r, there are two additional dimensions (t' and s')
            shape_add = (par.T, par.S)
        else:
            ndim_add = (np.newaxis,) #all other derivatives are simply a parameter vector dimension
            shape_add = (self.est.n_params,)

        dEV = np.zeros(self.sol.EV.shape + shape_add) #her tager den alle parametrene, selvom det kun er dem i utility vi regner på nu. kommer fejl senere.
        dv = np.zeros(P.shape + shape_add) #derivatives of choice-specific value functions

        #Now we can fill in dv for the last age. There are no continuation values, so only the marginal wage effect matters
        a = par.a_max - par.a_min
        dv[:, a, ...] = du[:, a, ...]

        #Next, the last age, dEV
        #Equation 41: Multiply P and dv for all s, and sum over s afterwards. We get positive values for all combinations of (slag x t),
        #because the expected value of the terminal period goes up no matter which of the wages we change through beta0. 
        dEV[:, a, :, ...] = np.sum(P[(slice(None), a, slice(None), slice(None)) + ndim_add] * dv[:, a, ...], axis=2)

        # Derivative for younger people in the last period: use continuation values from age + 1 in the same period.
        # We have to iterate on age, because the continuation value of someone with age a uses the continuation value of someone
        # with age a + 1. 
        while a > 0:
            a -= 1
            dv[:, a, par.T - 1, ...] = du[:, a, par.T - 1, ...] + par.rho * dEV[:, a + 1, par.T - 1, ...][np.newaxis, ...]
            #Here we match s in P with slag in dEV because the choice today becomes the lagged choice tomorrow
            dEV[:, a, par.T - 1, ...] = np.sum(P[(slice(None), a, par.T - 1, slice(None)) + ndim_add] * dv[:, a, par.T - 1, ...], axis=1)

        #Now we can perform proper time iteration for the remaining ages.
        t = par.T - 1
        while t > 0:
            t -= 1
            #Choice specific value function derivatives. 
            dv[:, :-1, t, ...] = du[:, :-1, t, ...] + par.rho * dEV[:, 1:, t + 1, :].swapaxes(0, 1)[np.newaxis, ...]
            dEV[:, :-1, t, ...] = np.sum(P[(slice(None), slice(None, -1), t, slice(None)) + ndim_add] * dv[:, :-1, t, ...], axis=2)

        # This concludes the calculation of the expected value function derivatives (wrt. beta0)
        return dv

        # def partial_ll(self, dlnP):
        #     """ The partial derivative of the log likelihood can be calculated from dlnP/dtheta which is the input to this function. 
        #         To calculate this derivative, dlnP is evaluated in all the data points from the estimation sample. 
        #     """
        #     d = self.est.simulated_data
        #     return np.sum(dlnP[d["slag"], d["a"] - self.par.a_min, d["t"], d["s"], ...], axis=0) / self.est.loglik_scale

    def dED_dr(self):
        """ This function is used only when we are only interested in the derivative wrt. r, i.e. for the equilibrum determination. """
        
        #unpack
        ts_indexes = self.ts_indexes
        par = self.par

        du = self.du("r")
        dv = self.dv_and_dEV(du, "r")
        dP = self.dP(dv, "r")

        dHsup = self.dHsup(dP, "r")
        dHdem = np.zeros((par.T, par.S, par.T, par.S))
        dHdem[ts_indexes[:, 0], ts_indexes[:, 1], ts_indexes[:, 0], ts_indexes[:, 1]] = \
            np.reshape(par.alpha1 * par.pY / (- np.square(par.r)), (par.T*par.S), order="F")
        dED = dHdem - dHsup
        return dED

    def ll_gradients(self):
        """ Calculates analytic gradients of the log likelihood function. 
            When the model is solved in partial mode, this assumes dr/dtheta = 0 and so only the partial effect dll/dtheta remains.
        """
        
        ts_indexes = self.ts_indexes
        par = self.par
        

        #noget omkring "utiltiy" her skal ændres senere til "sigma" osv også.
        if self.partial_model:
            du = self.du("utility")
            dv = self.dv_and_dEV(du, "utility") #så længe det er utility-parametre, kan denne køre på samme måde med alle parametrene (ligning 44)
            dlnP = 1/par.sigma * (dv - np.sum(self.sol.P[..., np.newaxis] * dv, axis=-2)[:, :, :, np.newaxis, :])
            return dlnP
        else:
            T = par.T
            S = par.S
            
            #du and dv enter both dED and dlnP
            du = self.du("r")
            dv = self.dv_and_dEV(du, "r")
            
            #Intermediate stuff for dED
            dP = self.dP(dv, "r")
            dHsup = self.dHsup(dP, "r")

            #lav en ny method: dHdem()
            #lav en ny method: dlnP(g)

            dHdem = np.zeros((par.T, par.S, par.T, par.S))
            dHdem[ts_indexes[:, 0], ts_indexes[:, 1], ts_indexes[:, 0], ts_indexes[:, 1]] = \
                np.reshape(par.alpha1 * par.pY / (- np.square(par.r)), (par.T*par.S), order="F")
            dED = dHdem - dHsup
            dED_dr_inv = np.linalg.inv(dED.swapaxes(0, 1).swapaxes(2, 3).reshape((T*S, T*S), order="C"))

            #To be outputted
            dlnP_dr = 1/par.sigma * (dv - np.sum(self.sol.P[..., np.newaxis, np.newaxis] * dv, axis=-3)[:, :, :, np.newaxis, :, :])
            
            #Theta derivatives
            #du and dv enter both deD and dlnP
            du = self.du("utility")
            dv = self.dv_and_dEV(du, "utility")

            dP = self.dP(dv, "utility")
            dHsup = self.dHsup(dP, "utility")
            dED = - dHsup #parameters cannot affect demand in the model so I simply leave it out

            #the two outputs
            dlnP_dtheta = 1/par.sigma * (dv - np.sum(self.sol.P[..., np.newaxis] * dv, axis=-2)[:, :, :, np.newaxis, :])
            dr_dtheta = np.matmul(dED_dr_inv, - dED.reshape((T*S, self.est.n_params), order="F"))

            return (dr_dtheta, dlnP_dtheta, dlnP_dr)

                #equation with label dll_dtheta in the paper
                #done dED_dr_inv = np.linalg.inv(self.partial_ED("r").swapaxes(0, 1).swapaxes(2, 3).reshape((T*S, T*S), order="C")) 
                #done dr_dtheta = np.matmul(dED_dr_inv, - self.partial_ED("theta").reshape((T*S, self.est.n_params), order="F"))
                # dll_dtheta = self.partial_ll(self.dlnP("theta"))
                # done (uden ll) dll_dr = self.partial_ll(self.dlnP("r"))
                
                #1-dimensional array with shape (nparams,)
                # return dll_dtheta + np.matmul(dll_dr.reshape((T*S), order="F"), dr_dtheta)
            
            # #Unpack
            # # w = self.sol.w
            # P = self.sol.P
            # # pte = self.est.pars_to_estimate
            # # theta_idx = self.est.theta_idx

            # #Containers for results
            # dEV = np.zeros(self.sol.EV.shape + (self.est.n_params,))
            # # du = np.zeros(P.shape + (self.est.n_params,))

            # #Assumes that only utility parameters are to be estimated FOR NOW.
            
            # du = self.du("utility")
            # dv = np.zeros(P.shape + (self.est.n_params,)) #derivatives of choice-specific value functions

            # #todo: VI KAN KUN ESTIMERE 1 PARAMETER LIGE NU, OG DET ER BETA0. 
            # #beta0 og xi-parametre skal sættes sammen, da alt fra dv kan vektoriseres. Men sigma kan ikke det samme.
            # #til sidst skal det hele så samles i en enkelt dlnP. 
            # #JEG udvikler dette når jeg har skrevet overleaf omkring ligevægtsgradienter OG sendt dem til Bertel.
            # # parameter = "beta0"
            # # for parameter in pte:
            # #Derivative of utility wrt. beta0. Calculated for all state space points once-and-for-all

            # # du[..., pte[parameter], theta_idx[parameter]] = w[..., pte[parameter]] * self.par.ages[np.newaxis, :, np.newaxis, np.newaxis] / getattr(self.par.scale, parameter)

            # a = self.par.a_max - self.par.a_min

            # #Now we can fill in dv for the last age. There are no continuation values, so only the marginal wage effect matters
            # dv[:, a, ...] = du[:, a, ...]

            # #Next, the last age, dEV
            # #Equation 41: Multiply P and dv for all s, and sum over s afterwards. We get positive values for all combinations of (slag x t),
            # #because the expected value of the terminal period goes up no matter which of the wages we change through beta0. 
            # dEV[:, a, :, :] = np.sum(P[:, a, :, :, np.newaxis] * dv[:, a, ...], axis=2)

            # # Derivative for younger people in the last period: use continuation values from age + 1 in the same period.
            # # We have to iterate on age, because the continuation value of someone with age a uses the continuation value of someone
            # # with age a + 1. 
            # while a > 0:
            #     a -= 1
            #     dv[:, a, self.par.T - 1, :, :] = du[:, a, self.par.T - 1, :, :] + self.par.rho * dEV[:, a + 1, self.par.T - 1, :][np.newaxis, :, :]
            #     #Here we match s in P with slag in dEV because the choice today becomes the lagged choice tomorrow
            #     dEV[:, a, self.par.T - 1, :] = np.sum(P[:, a, self.par.T - 1, :, np.newaxis] * dv[:, a, self.par.T - 1, :, :], axis=1)

            # #Now we can perform proper time iteration for the remaining ages.
            # t = self.par.T - 1
            # while t > 0:
            #     t -= 1
            #     #Choice specific value function derivatives. 
            #     dv[:, :-1, t, :, :] = du[:, :-1, t, :, :] + self.par.rho * dEV[:, 1:, t + 1, :].swapaxes(0, 1)[np.newaxis, ...]
            #     dEV[:, :-1, t, :] = np.sum(P[:, :-1, t, :, np.newaxis] * dv[:, :-1, t, :, :], axis=2)

            # # This concludes the calculation of the expected value function derivatives (wrt. beta0)
            # # Next, calculate the derivatives of the choice probabilities. 
            # # This is made easier by the fact that we have already calculated dv.

            # dlnP = np.zeros((P.shape) + (self.est.n_params,)) #can be commented since it actually not used.
            # #this is the dv of the choice minus a log sum, where after summing over k, we create the s dimension again
            # dlnP = 1/self.par.sigma * (dv - np.sum(P[..., np.newaxis] * dv, axis=-2)[:, :, :, np.newaxis, :])

            # return dlnP
