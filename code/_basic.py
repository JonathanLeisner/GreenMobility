"""  Defines the class that runs the Green Mobility (GM) project code."""
#%%

from collections import OrderedDict
import numpy as np

#User-made modules
import _gradients as gradients


#%%

class basic:
    """ This class is meant to be called only as a parent class to either the estimation or the
    counterfactual simulation classes. It defines all the functions that are common to these two. """
    # pylint: disable=no-member

    def setup_ages(self, version="standard"):
        par = self.par
        if version == "standard":
            par.a_max = 65
        elif version == "simple":
            par.a_max = 40
        par.a_min = 30 #youngest cohort
        par.ages = np.arange(par.a_min, par.a_max + 1)
        par.N_ages = par.a_max + 1 - par.a_min

    def setup_tenure(self, version="standard"):
        par = self.par
        if version == "standard":
            par.ten_max = 7
        elif version == "simple":
            par.ten_max = 3
        par.tenures = np.arange(1, par.ten_max + 1)

    def setup_statespace(self, version="standard"):
        """ This function defines the attributes related to the state space. It allows setting 'version'
            to a string to get a smaller statespace when testing new features.
            version must be either 'standard' or 'simple'.
            Also defines the maximum value that tenure (ten) can attain.
        """
        self.setup_sectors(version)
        self.setup_ages(version)
        self.setup_periods(version)
        self.setup_tenure(version)

    def set_parameter_values(self, **kwargs):
        """
        Set the default values for all parameters.
        #Meaning true parameters are x times smaller than the values in the code
        """
        par = self.par
        #Parameters that are calibrated
        par.rho = 0.96

        #### Parameters to be estimated later ####
        #Gumbel shock scale parameter
        par.sigma = np.array([2.0])

        #Human capital function parameters
        par.beta0 = np.array([0.2, 0.6, 1.2, 0.8, 1, 1])[:par.S - 1] #1-dimensional over sectors. 
        par.scale.beta0 = 100

        par.beta1 = np.array([-0.1, -0.2, -0.3, -0.4, -0.5])[:par.S - 1]
        par.scale.beta1 = 10000

        # par.beta2 = np.array([0.0, 0.0, 0.0, 0.0, 0.0])[:par.S - 1]
        par.beta2 = np.array([0.768, 0.768, 0.768, 0.768, 0.473, 1.0])[:par.S - 1]
        par.scale.beta2 = 10

        #Unemployment utility parameters
        par.phi0 = np.array([-2.0], dtype=float) #-2.0
        par.scale.phi0 = 100
        par.phi1 = np.array([0.7], dtype=float) #0.7
        par.scale.phi1 = 1000
        par.phi2 = np.array([-2.0])
        par.scale.phi2 = 1

        #Moving parameters
        par.xi_in = np.array([-10.0, 5.0, 6.0, 7.0, 4.0, 4.0])[:par.S - 1]
        par.xi_out = np.array([5.0, 0.0, 4.0, 7.0, 3.0, 3.0])[:par.S].transpose()
        par.scale.xi_in = 10
        par.scale.xi_out = 10

        #Individual characteristics shifter of moving costs
        par.kappa0 = np.array([2.0])
        par.scale.kappa0 = 100
        par.kappa1 = np.array([-2.0])
        par.scale.kappa1 = 10000

        #Define parameter groups. These are used in self.gradient().
        par.groups = OrderedDict([("sigma", ["sigma"]), 
                                  ("utility", ["beta0", "beta1", "xi_in", "xi_out", "kappa0", "kappa1", "phi0", "phi1", "beta2"])])
        par.all_parameters = [p for sublist in par.groups.values() for p in sublist]
        # #Replace default values by those given explicitly to the method
        for k, v in kwargs.items():
            setattr(par, k, v)

        #Post-setup stuff (e.g. calculating M based on xi-parameters)
        self.c_switching_cost_matrix()

    def c_m1(self):
        """ Calculate the m1 matrix of switching costs, a component of M."""
        par = self.par
        #add cost of going OUT of a sector with the cost of going IN to a sector:
        m1 = np.exp(np.sum(np.meshgrid(par.xi_in / par.scale.xi_in, par.xi_out / par.scale.xi_out), axis=0)) 
        m1 = np.append(np.zeros(par.S)[:, np.newaxis], m1, axis=1) #Switching into s=0 is costless.
        np.fill_diagonal(m1, 0) #Not switching sector is costless
        return m1

    def c_m2(self):
        """ Calculate the switching cost shifter m2 that depends on individual characteristics. m2 is a component of M."""
        par = self.par
        return np.exp(np.array(par.ages * par.kappa0 / par.scale.kappa0 + np.square(par.ages) * par.kappa1 / par.scale.kappa1))

    def c_switching_cost_matrix(self):
        """ Takes the xi-parameters and constructs the matrix of switching costs."""

        #Switching cost (M) has dimensions (s_{t-1} x age x s_{t}) i.e. (s_out x age x s_in) where age component comes from m2 in M = m1 * m2
        m1 = self.c_m1()
        m2 = self.c_m2()
        self.par.M = m1[:, np.newaxis, :] * m2[np.newaxis, :, np.newaxis]

    # def l_psi(self):
    #     """ Loads the psi parameter from data. There is one value for each year, corresponding to the base
    #         year used for calibration """ 
    #     return np.load("../data_prep/output/psi.npy")

    def l_e(self):
        """ Loads the emission intensity parameters (symbol e)."""
        return np.load("../data_prep/output/e.npy")

    def l_tau(self):
        """ Loads the tau parameter, i.e. the average cost of emissions in the data """
        return np.load("../data_prep/output/tau.npy")

    def l_rK(self):
        """ Loads the return to physical capital (symbol rK) from the data """
        return np.load("../data_prep/output/rK.npy")

    def l_Z(self):
        """ Loads emissions (symbol Z) from data.""" 
        return np.load("../data_prep/output/Z.npy")

    def allocate_worker(self):
        """ Allocate empty containers for solution objects. 
        P contains the CCPs and EV contains expected value functions"""
        #todo: Update docstring
        sol = self.sol
        par = self.par

        sol.P = np.zeros((par.ten_max, par.S, par.N_ages, par.T, par.S)) #4 dimensions: slag, a, t, s (choice)
        sol.EV = np.zeros((par.ten_max, par.S, par.N_ages, par.T)) #3 dimensions: slag, a, t
        sol.v = np.zeros((par.ten_max, par.S, par.N_ages, par.T, par.S))

    def precompute(self):
        """ Calculates the wages offered in each sector at each point in the state space. 
        The result is sol.w which we can then use directly, rather than using r and H. """
        #Precompute (for all ages) human capital #H has dimensions: (ten, slag, s, a)
        self.precompute_H()
        #Precompute wages (a x s x t)
        self.precompute_w()

    def precompute_H(self):
        """ Calculate all possible values of the human capital function. Dimensions: ten x slag x a x s"""
        # self.sol.H = np.exp(self.par.beta0[np.newaxis, :]/self.par.scale.beta0 * (self.par.ages[:, np.newaxis] - self.par.a_min))
        par = self.par
        self.sol.H = np.exp(
            (par.beta0[np.newaxis, :]/par.scale.beta0 * (par.ages[:, np.newaxis] - par.a_min))[np.newaxis, np.newaxis, ...] + \
            par.beta1[np.newaxis, :]/par.scale.beta1 * (np.square(par.ages[:, np.newaxis] - par.a_min))[np.newaxis, np.newaxis, ...] + \
            ((par.beta2[np.newaxis, :]/par.scale.beta2 * par.tenures[:, np.newaxis])[:, np.newaxis, :] * \
            np.append(np.zeros(par.S - 1)[np.newaxis, :], np.eye(par.S - 1), axis=0)[np.newaxis, ...])[..., np.newaxis, :]
            )

    def precompute_w(self):
        """ Calculate all possible values of wages. w = r*H. Dimensions: (ten, slag, a, t, s)"""
        # par = self.par
        # w = np.zeros((par.ten_max, par.S, par.N_ages, par.T, par.S))
        # #Non-employment
        # # w[:, :, 0] = (par.phi0/par.scale.phi0 * (par.ages - par.a_min) + par.phi1/par.scale.phi1 * 
        # #               np.square(par.ages - par.a_min) + par.phi2/par.scale.phi2)[:, np.newaxis]
        # w[:, :, 0] = (par.phi0/par.scale.phi0 * (par.ages - par.a_min) + par.phi1/par.scale.phi1 * 
        #               np.square(par.ages - par.a_min))[:, np.newaxis]
        # #Productive sectors
        # w[:, :, 1:] = self.sol.H[:, np.newaxis, :] * par.r[np.newaxis, :, :]
        par = self.par
        w = np.zeros((par.ten_max, par.S, par.N_ages, par.T, par.S))
        #Non-employment
        # w[:, :, 0] = (par.phi0/par.scale.phi0 * (par.ages - par.a_min) + par.phi1/par.scale.phi1 * 
        #               np.square(par.ages - par.a_min) + par.phi2/par.scale.phi2)[:, np.newaxis]
        w[..., 0] = (par.phi0/par.scale.phi0 * (par.ages - par.a_min) + par.phi1/par.scale.phi1 * 
                    np.square(par.ages - par.a_min))[..., np.newaxis]

        #Productive sectors
        w[..., 1:] = self.sol.H[..., np.newaxis, :] * par.r
        self.sol.w = w

        # self.sol.w = self.sol.H[:, np.newaxis, :] * self.par.r[np.newaxis, :, :]

    @staticmethod
    def closed_forms(sigma, arr):
        """Calculate the closed-form solutions for EV and CCP. """
        if arr.ndim == 1:
            constant = np.max(arr / sigma)
            sumexp = np.sum(np.exp(arr / sigma - constant))
            EV = sigma * (constant + np.log(sumexp))
            P = np.exp(arr / sigma - constant) / sumexp
        if arr.ndim > 1:
            constant = np.max(arr / sigma, axis=-1) #max for a given state space point, across choices. This becomes last axis when I move choice.
            sumexp = np.sum(np.exp(arr / sigma - constant[..., np.newaxis]), axis=-1) #sumexp is used in both formulas, so we only calculate it once.
            EV = sigma * (constant + np.log(sumexp))
            P = np.exp(arr / sigma - constant[..., np.newaxis]) / sumexp[..., np.newaxis]
        return (EV, P)

    def solve_worker(self):
        """ Solve the worker's problem by backwards induction for a given set of skill prices (r). 
            First we solve the last period using static expectations. Then we perform backwards 
            iteration from the remaining periods until period t = 0 using rational expectations. """

        #Unpack solution objects
        par = self.par
        sol = self.sol
        w = sol.w #offered wages.
        M = par.M #utility switching cost
        P = sol.P #conditional choice probabilities
        EV = sol.EV #expected value functions
        v = sol.v #alternative-specific value functions

        #Indicator/dummy array with 1 whenever a combination of slag and s are equal.
        stay = np.eye(par.S)

        #PART 0: Precompute EV at the retirement age for all periods (no continuation values)
        a = par.N_ages - 1
        v[..., a, :, :] = w[..., a, :, :] - M[..., a, np.newaxis, :]
        (EV[..., a, :], P[..., a, :, :]) = self.closed_forms(par.sigma, v[..., a, :, :])

        #PART I (age iteration (64 -> 30) within the terminal period, static expectation)
        t = par.T - 1
        for a in reversed(np.arange(par.N_ages - 1)):
            #Staying adds +1 to tenure and lagged sector is matched to previously chosen one
            EVstay = np.concatenate((EV[1:, :, a + 1, t], EV[-1, :, a + 1, t][np.newaxis, ...]))[:, np.newaxis, :]
            EVnonstay = EV[0, :, a + 1, t][np.newaxis, np.newaxis, ...] #nonstaying means tenure gets reset (index 0).
            # My choice of sector today enters EV tomorrow as slag. 
            # Therefore the last dimension of w (the choice) must match the first dimension of EV (slag) 
            v[..., a, t, :] = w[..., a, t, :] - M[..., a, :] + par.rho * (EVstay * stay + EVnonstay * (1 - stay))
            (EV[..., a, t], P[..., a, t, :]) = self.closed_forms(par.sigma, v[..., a, t, :])

        #PART II (time iteration from T - 2 to 0.)
        t = par.T - 2
        for t in reversed(np.arange(par.T - 1)):
            # The value of an alternative is wage minus mobility costs + next period's continuation value discounted. 
            # Again, the transpose() and newaxis on EV makes sure that the choice dimension of w and M lines up with the slag dimension of EV.
            EVstay = np.swapaxes(np.concatenate((EV[1:, :, 1:, t + 1], EV[-1, :, 1:, t + 1][np.newaxis, ...])), -1, -2)[..., np.newaxis, :, :]
            EVnonstay = np.swapaxes(EV[0, :, 1:, t + 1], -1, -2)[np.newaxis, np.newaxis, ...]
            v[..., :-1, t, :] = w[..., :-1, t, :] - M[..., :-1, :] + \
                par.rho * (EVstay * stay[:, np.newaxis, :] + EVnonstay * (1 - stay[:, np.newaxis, :]))

            # v[:, :-1, t, :] = w[np.newaxis, :-1, t, :] - M[:, :-1, :] + par.rho * EV[:, 1:, t + 1].transpose()[np.newaxis, :, :]
            (EV[..., :-1, t], P[..., :-1, t, :]) = self.closed_forms(par.sigma, v[..., :-1, t, :])

        assert np.all(np.isclose(np.sum(P, axis=-1), 1)) 

    def dED_dr(self, r=None):
        """ Calculates the derivative dED/dr. If a vector (1-D or 2-D) of skill prices is given, the model is solved and simulated
            for these skill prices before the derivative is calculated. If r is None, the derivative is calculated from the currently 
            stored solution and vector of prices.
        """
        if r is not None:
            if r.ndim == 1:
                self.par.r = r.reshape(self.par.r.shape, order="F")
            elif r.ndim == 2:
                self.par.r = r
            self.precompute_w() #only wages need to be updated when the only thing we change is r
            self.solve_worker()
            self.simulate()

        #Draws upon the gradient class in the gradient module:
        return gradients.gradients(self.par, self.sol, self.sim, self.est, self.partial_model).dED_dr()

    def ED_from_r(self, r_1d):
        """ Calculates the excess labor demands from a new vector (1-D) of skill prices. 
            In this vector, index t moves the fastest, implying for example 
            that the first entries vary the year and keep the sector fixed. 
        """
        self.par.r = r_1d.reshape(self.par.r.shape, order="F") #Update rental prices
        self.precompute_w()
        self.solve_worker()
        self.simulate()
        #Calculate excess labor demand
        ED = self.c_ED()
        return ED.reshape((self.par.T * (self.par.S - 1)), order="F")

#%% Test new feautures of GM before implementing them in the class