#"""  Defines the class that calculates gradients for the Green Mobility (GM) project code."""
#%%

import numpy as np

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

        if g is "sigma":
            dP = (P * ((dv * self.par.sigma - self.sol.v) / np.square(self.par.sigma) - \
                  np.sum(self.sol.P * (dv * self.par.sigma - self.sol.v) / np.square(self.par.sigma), axis=-1)[..., np.newaxis])
                 )[..., np.newaxis]
        else:
            ndim_add = {"r":(np.newaxis, np.newaxis), "utility":(np.newaxis,)}.get(g)
            dP = P[(...,) + ndim_add]/par.sigma * (dv -  np.sum(P[(...,) + ndim_add] * dv, axis=3)[:, :, :, np.newaxis, ...])
        
        return dP


    def dHsup(self, dP, g):
        """ Calculate the partial derivative of human capital supply with respect to either skill prices or theta. 
            The input g controls this. 
            Expressions for D and dP must be supplied as well. Implicitly, this function assumes that D is constant and so D' can be ignored.
        """

        ndim_add = {"r":(np.newaxis, np.newaxis), "utility":(np.newaxis,), "sigma":(np.newaxis,)}.get(g)
        dHsup =  self.par.MASS[(slice(None), np.newaxis) + ndim_add] * \
                        np.sum(self.est.D[(..., np.newaxis) + ndim_add] * dP, axis=tuple(np.arange(0, len(self.sol.EV.shape) - 1)))
        return dHsup

    def du(self, g):
        """ Calculates du with respect to some group of parameters, e.g. utility parameters, sigma or r. This choice determines the functional form
            implicitly (or if you look in the code actually, explicitly) used to calculate du. 
            The dimensions of du is state space + 1 when g is not "r" and state space + 2 when g is "r".    
        """

        #unpack
        par = self.par
        
        if g == "utility":
            #unpack
            pte = self.est.pars_to_estimate
            theta_idx = self.est.theta_idx
            w = self.sol.w

            du = np.zeros(self.sol.P.shape + (self.est.n_params[g],))
            for para in [k for k in pte.keys() if k in par.groups[g]]: #could this be precomputed? Yes, in setup_estimation perhaps. 
                if para == "beta0":
                    du[..., pte[para], theta_idx[para][g]] = w[..., pte[para]] * par.ages[np.newaxis, :, np.newaxis, np.newaxis] / getattr(par.scale, para)
                elif para == "xi_in":
                    du[:, :, :, pte[para], theta_idx[para][g]] = - par.M[:, :, np.newaxis, pte[para]] / getattr(par.scale, para)
                elif para == "xi_out":
                    du[pte[para], :, :, :, theta_idx[para][g]] = - par.M[pte[para], :, np.newaxis, :] / getattr(par.scale, para)
                elif para == "kappa0":
                    du[:, :, :, :, theta_idx[para][g]] = - par.M[:, :, np.newaxis, :, np.newaxis] * \
                                                      par.ages[np.newaxis, :, np.newaxis, np.newaxis, np.newaxis] / getattr(par.scale, para)
                elif para == "kappa1":
                    du[:, :, :, :, theta_idx[para][g]] = - par.M[:, :, np.newaxis, :, np.newaxis] * \
                                                      np.square(par.ages[np.newaxis, :, np.newaxis, np.newaxis, np.newaxis]) / getattr(par.scale, para)
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

    def dv_and_dEV_sigma(self):
        #unpack
        par = self.par
        P = self.sol.P
        v = self.sol.v
        EV = self.sol.EV

        dEV = np.zeros(EV.shape) #with sigma, no need to add an extra dimension, since it has length 1 anyways and so is redundant.
        dv = np.zeros(v.shape)

        a = par.a_max - par.a_min
        dEV[:, a, :] = EV[:, a, :] / par.sigma + np.sum(P[:, a, ...] * ((dv[:, a, ...] * par.sigma - v[:, a, ...]) / par.sigma), axis=-1)

        #age iteration for the final period
        while a > 0:
            a -= 1
            dv[:, a, par.T - 1, :] = par.rho * dEV[np.newaxis, :, a, par.T - 1]
            dEV[:, a, par.T - 1] = EV[:, a, par.T - 1] / par.sigma + \
                                np.sum(P[:, a, par.T - 1, :] * ((dv[:, a, par.T - 1, :] * par.sigma - v[:, a, par.T - 1, :]) / par.sigma), axis=-1)

        #time iteration for all remaining ages
        t = par.T - 1
        while t > 0:
            t -= 1
            dv[:, :-1, t, :] = par.rho * dEV[:, 1:, t + 1].swapaxes(0, 1)[np.newaxis, ...]
            dEV[:, :-1, t] = EV[:, :-1, t] / par.sigma + np.sum(P[:, :-1, t, :] * ((dv[:, :-1, t, :] * par.sigma - v[:, :-1, t, :]) / par.sigma), axis=-1)

        return dv

    def dv_and_dEV(self, du, g):
        """ Calculates dv (alternative-specific value function derivatives) and dEV (expected value function derivatives) from du. 
            g specifies whether the derivative is with respect to r or parameters, since this determines whether we have
            to add two dimensions or one (since r is defined over s and t while parameters are collected in one dimension). 
            The function only returns dv since dEV is never necessary to know on its own if we know dv.
            The calculation of derivatives follows the same backwards recursion structure as the solution to the worker's problem.
            If g == "sigma", the functional forms change and I use the method dv_and_dEV_sigma to calculate the derivative dv.
        """

        #unpack
        par = self.par
        P = self.sol.P

        if g is "r":
            ndim_add = (np.newaxis, np.newaxis)
            shape_add = (par.T, par.S)
        elif g is "utility":
            ndim_add = (np.newaxis,)
            shape_add = (self.est.n_params[g],)

        dEV = np.zeros(self.sol.EV.shape + shape_add)
        dv = np.zeros(P.shape + shape_add)

        #Age 65
        a = par.a_max - par.a_min
        dv[:, a, ...] = du[:, a, ...]
        dEV[:, a, :, ...] = np.sum(P[(slice(None), a, slice(None), slice(None)) + ndim_add] * dv[:, a, ...], axis=2)

        # a < 65, t = terminal period. Use continuation values from age + 1 in the same period. We have to iterate on age
        # because the continuation value of someone with age a uses the continuation value of someone with age a + 1. 
        while a > 0:
            a -= 1
            dv[:, a, par.T - 1, ...] = du[:, a, par.T - 1, ...] + par.rho * dEV[:, a + 1, par.T - 1, ...][np.newaxis, ...]
            #Here we match s in P with slag in dEV because the choice today becomes the lagged choice tomorrow
            dEV[:, a, par.T - 1, ...] = np.sum(P[(slice(None), a, par.T - 1, slice(None)) + ndim_add] * dv[:, a, par.T - 1, ...], axis=1)

        #Time iteration for remaining ages
        t = par.T - 1
        while t > 0:
            t -= 1
            dv[:, :-1, t, ...] = du[:, :-1, t, ...] + par.rho * dEV[:, 1:, t + 1, :].swapaxes(0, 1)[np.newaxis, ...]
            dEV[:, :-1, t, ...] = np.sum(P[(slice(None), slice(None, -1), t, slice(None)) + ndim_add] * dv[:, :-1, t, ...], axis=2)
        return dv

    def dED_dr(self):
        """ Calculates the derivatives of the excess labor demands with respect to skill prices.
            The method is used only for calculating derivatives for the equilibrum conditions, and not log likelihood. """
        du = self.du("r")
        dv = self.dv_and_dEV(du, "r")
        dP = self.dP(dv, "r")

        dED = self.dHdem_dr() - self.dHsup(dP, "r")
        return dED

    def dHdem_dr(self):
        """ Calculate the derivative of labor demand with respect to skill prices r."""
        par = self.par
        ts_indexes = self.ts_indexes
        dHdem = np.zeros((par.T, par.S, par.T, par.S))
        dHdem[ts_indexes[:, 0], ts_indexes[:, 1], ts_indexes[:, 0], ts_indexes[:, 1]] = \
            np.reshape(par.alpha1 * par.pY / (- np.square(par.r)), (par.T*par.S), order="F")
        return dHdem

    def dlnP(self, dv, g):
        """ Calculate the derivatives of log choice probabilities, the only component of the log likelihood function. 
            They are calculated from the derivatives of alternative-specific value functions (dv). 
            The returned array has dimensions (state space + S + 1)
        """
        if g is "r":
            return 1/self.par.sigma * (dv - np.sum(self.sol.P[..., np.newaxis, np.newaxis] * dv, axis=-3)[:, :, :, np.newaxis, :, :])
        elif g is "utility":
            return 1/self.par.sigma * (dv - np.sum(self.sol.P[..., np.newaxis] * dv, axis=-2)[:, :, :, np.newaxis, :])
        elif g is "sigma":
            return ((dv * self.par.sigma - self.sol.v) / np.square(self.par.sigma) - \
                    np.sum(self.sol.P * (dv * self.par.sigma - self.sol.v) / np.square(self.par.sigma), axis=-1)[..., np.newaxis])[..., np.newaxis]

    def ll_gradients(self):
        """ Calculates analytic gradients of the log likelihood function. 
            When the model is solved in partial mode, this assumes dr/dtheta = 0 and so only the partial effect dll/dtheta remains.
        """
        #todo: når sigma kommer ind skal der laves om her.
        par = self.par

        #noget omkring "utility" her skal ændres senere til "sigma" osv også.
        if self.partial_model:
            dlnP = []
            for g in self.est.groups_to_estimate:
                if g is "sigma":
                    dv = self.dv_and_dEV_sigma()
                else:
                    du = self.du(g)
                    dv = self.dv_and_dEV(du, g)
                dlnP.append(self.dlnP(dv, g))
            return np.concatenate(dlnP, axis=-1)

            # du = self.du("utility")
            # dv = self.dv_and_dEV(du, "utility")
            # return self.dlnP(dv, "utility")
        else:
            T = par.T
            S = par.S
            
            du = self.du("r")
            dv = self.dv_and_dEV(du, "r")
            dP = self.dP(dv, "r")
            dHsup = self.dHsup(dP, "r")
            dED = self.dHdem_dr() - dHsup
            dED_dr_inv = np.linalg.inv(dED.swapaxes(0, 1).swapaxes(2, 3).reshape((T*S, T*S), order="C"))
            dlnP_dr = self.dlnP(dv, "r")
            
            #Theta derivatives
            dED = []
            dlnP = []
            for g in self.est.groups_to_estimate:
                if g is "sigma":
                    dv = self.dv_and_dEV_sigma()
                else:
                    du = self.du(g)
                    dv = self.dv_and_dEV(du, g)
                dlnP.append(self.dlnP(dv, g))
                dP = self.dP(dv, g)
                dED.append(- self.dHsup(dP, g)) #parameters cannot affect demand in the model so I simply leave it out (dED = 0 - dHsup)
            dED = np.concatenate(dED, axis=-1) 
            dlnP_dtheta = np.concatenate(dlnP, axis=-1)
            
            dr_dtheta = np.matmul(dED_dr_inv, - dED.reshape((T*S, self.est.n_params["full"]), order="F"))

            return (dr_dtheta, dlnP_dtheta, dlnP_dr)
