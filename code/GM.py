"""  Defines the class that runs the Green Mobility (GM) project code."""
#%%

# from operator import index
# from turtle import update
import importlib
from types import SimpleNamespace
from collections import OrderedDict
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import optimize
from numpy.random import default_rng
# from mpl_toolkits import mplot3d

np.set_printoptions(precision=3)
sns.set_theme()
import warnings
warnings.filterwarnings("error", category=RuntimeWarning)

#User-made modules
import util
import gradients
import output

def update_attrs(ns, **kwargs):
    for k, v in kwargs.items():
        setattr(ns, k, v)

def get_key_from_value(dict, value):
    keys = [k for k in dict.keys() if value in dict[k]]
    assert len(keys) == 1
    return keys[0]




#%%

class GM:
    """ Class to solve the model presented in Leisner (2023). 
        Generally, the notation c_[name] for methods means construct/create/calculate an object or attribute called [name].
    """
    def __init__(self, name = None, endo_D=True):
        self.sol = SimpleNamespace()
        self.par = SimpleNamespace()
        self.par.scale = SimpleNamespace() #collects scalar scales for all parameters.
        self.sim = SimpleNamespace()
        self.est = SimpleNamespace()
        self.results = SimpleNamespace()
        self.results.tables = SimpleNamespace()
        self.results.figures = SimpleNamespace()
        self.diag = SimpleNamespace()
        self.diag.tables = SimpleNamespace()

        if name is None:
            self.name = "standard"
        else:
            self.name = name


        self.version = "v3"
        self.rng = default_rng(123456) #seed

        self.partial_model = False
        self.agrad_ED0 = True
        self.agrad_quadED = True

        self.est.default_method = "BFGS"
        update_attrs(self.sim, endo_D=endo_D)
        # self.setup_sim(endo_D=endo_D)

    # def setup_sim(self, **kwargs):
    #     """ Sets up settings needed to simulate. 
    #         The keywords that can be used are:
    #         endo_D : {True/False}. Switching between endogenous and exogenous D (density).
    #         agrad_quadED : {True/False}. Switch between numeric jacobian and analytic jacobian for outer GE problem. """

    #     for k, v in kwargs.items():
    #         setattr(self.sim, k, v)


    def setup_statespace(self, simple_statespace):
        """ This function defines the attributes related to the state space. It allows setting 'simple_statespace'
            to True to get a smaller statespace when testing new features. 
        """
        
        #Unpack
        par = self.par
        
        if simple_statespace:        
            par.T = 4 #periods
            par.a_max = 40 #eldest cohort
        else:
            par.T = 8
            par.a_max = 65
        
        #Sectors
        par.sector_names = ["Unemployment", "Clean", "Dirty"]
        par.S = len(par.sector_names)

        #Age
        par.a_min = 30 #youngest cohort
        par.ages = np.arange(par.a_min, par.a_max + 1)
        par.N_ages = par.a_max + 1 - par.a_min

    def setup(self, simple_statespace=False, **kwargs):
        """ Initiates parameter values at their default values. Also populates meta settings. Specify named
            'keyword=value' to the function to replace default values by those specified in 'value'. 
        """
        
        par = self.par
        self.setup_statespace(simple_statespace)

        #Fixed quantities from data
        par.MASS = self.load_MASS()
        par.pY = self.load_nominal_output()

        #Parameters that are calibrated
        par.rho = 0.96

        #Should be loaded from data later
        par.alpha1 = np.repeat(np.array([0.2, 0.2, 0.3, 0.1])[:par.S][np.newaxis, :], par.T, axis=0)
        par.alpha2 = np.repeat(np.array([0.2, 0.2, 0.2, 0.1])[:par.S][np.newaxis, :], par.T, axis=0)
        par.alpha3 = 1 - par.alpha1 - par.alpha2

        #Parameters to be estimated later

        #Gumbel shock scale parameter
        par.sigma = 0.8

        #Human capital function parameters
        par.beta0 = np.array([0.2, 0.6, 1.2, 0.8])[:par.S] #1-dimensional over sectors. 
        par.scale.beta0 = 100 #Meaning true parameters are x times smaller than the values in the code

        #Switching costs (direct utility costs)
        if self.name == "high_switchingcost":
            par.xi_in = 20 * np.array([0.04, 0.05, 0.06, 0.07])[:par.S]
            par.xi_out = 20 * np.array([0.03, 0.06, 0.09, 0.07])[:par.S].transpose()    
        else:
            par.xi_in = np.array([0.04, 0.05, 0.06, 0.07])[:par.S]
            par.xi_out = np.array([0.03, 0.06, 0.09, 0.07])[:par.S].transpose()

        #Individual characteristics scaler of moving costs
        par.kappa0 = 1
        par.scale.kappa0 = 100
        par.kappa1 = -0.001
        par.scale.kappa1 = 1000

        #Define parameter groups. These are used in self.gradient().
        par.groups = {"sigma":["sigma"], "utility":["beta0", "xi_in", "xi_out", "kappa0", "kappa1"], "r":["r"]}

        # #Replace default values by those given explicitly to the method
        for k, v in kwargs.items():
            setattr(self.par, k, v)

        #Post-setup stuff (e.g. calculating M based on xi-parameters)
        self.c_switching_cost_matrix()

        #Initial values for skill prices (r)
        par.r = np.linspace(np.arange(1, par.S+1), np.arange(1, par.S+1)[::-1], par.T, axis=0)

        # self.diag.tables.skillprice_converge = \
        #     pd.DataFrame(np.nan, 
        #                  index=pd.MultiIndex.from_product([self.par.sector_names, range(self.par.T)], 
        #                                                   names=["Sector", "Year"]), 
        #                  columns=pd.MultiIndex.from_product([["r0", "r1"], range(0, self.sol.maxiter)], 
        #                                                     names=["Pre/post", "Iteration"]))

    def c_m1(self):
        par = self.par
        m1 = np.exp(np.sum(np.meshgrid(par.xi_in, par.xi_out), axis=0)) #add cost of going OUT of a sector with the cost of going IN to a sector.
        np.fill_diagonal(m1, 0) #Not switching sector is costless
        return m1

    def c_m2(self):
        par = self.par
        return np.exp(np.array(par.ages * par.kappa0 / par.scale.kappa0 + np.square(par.ages) * par.kappa1 / par.scale.kappa1))

    def c_switching_cost_matrix(self):
        """ Takes the xi-parameters and constructs the matrix of switching costs."""

        #Switching cost (M) has dimensions (s_{t-1} x age x s_{t}) i.e. (s_out x age x s_in) where age component comes from m2 in M = m1 * m2
        m1 = self.c_m1()
        m2 = self.c_m2()
        self.par.M = m1[:, np.newaxis, :] * m2[np.newaxis, :, np.newaxis]

    def load_MASS(self, version="constant"):
        """ Creates a time series of MASS."""
        #for testing, this makes the population double halfway through the sample
        if version == "constant":
            return np.ones(self.par.T)
        else:
            return np.concatenate([np.ones(int(round(self.par.T/2, 0))), np.ones(self.par.T - int(round(self.par.T/2, 0))) + 1])

    def load_nominal_output(self):
        """ Creates a time series of nominal output, i.e. p^Y times Y in the paper."""
        return np.repeat((np.ones(self.par.S))[np.newaxis, :], self.par.T, axis=0)

    def allocate(self):
        """ Allocate empty containers for solution objects. 
        P contains the CCPs and EV contains expected value functions"""
        sol = self.sol
        par = self.par

        sol.P = np.zeros((par.S, par.N_ages, par.T, par.S)) - 1 #4 dimensions: slag, a, t, s (choice)
        sol.EV = np.zeros((par.S, par.N_ages, par.T)) - 99 #3 dimensions: slag, a, t

    def precompute(self):
        """ Calculates the wages offered in each sector at each point in the state space. 
        The result is sol.w which we can then use directly, rather than using r and H. """
        #Precompute (for all ages) human capital #H has dimensions: s x a
        self.precompute_H()
        #Precompute wages (a x s x t) = (36, 3, 4)
        self.precompute_w()

    def precompute_H(self):
        """ Calculate all possible values of the human capital function. Dimensions: a x s"""
        self.sol.H = np.exp(self.par.beta0[np.newaxis, :]/self.par.scale.beta0 * self.par.ages[:, np.newaxis])

    def precompute_w(self):
        """ Calculate all possible values of wages. w = r*H. Dimensions: a x t x s"""
        self.sol.w = self.sol.H[:, np.newaxis, :] * self.par.r[np.newaxis, :, :]

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

        #PART 0: Precompute EV at the retirement age for all periods (no continuation values)
        a = par.N_ages - 1        
        (EV[:, a, :], P[:, a, :, :]) = self.closed_forms(par.sigma, w[np.newaxis, a, :, :] - M[:, a, np.newaxis, :])

        #PART I (age iteration (64 -> 30) within the terminal period, static expectation)
        t = par.T - 1
        for a in reversed(np.arange(par.N_ages - 1)):
            V_alternatives = w[np.newaxis, a, t, :] - M[:, a, :] + par.rho * EV[np.newaxis, :, a + 1, t]
            # My choice of sector today enters EV tomorrow as slag. 
            # Therefore the last dimension of w (the choice) must match the first dimension of EV (slag) 
            (EV[:, a, t], P[:, a, t, :]) = self.closed_forms(par.sigma, V_alternatives)

        #PART II (time iteration from T - 2 to 0.)
        t = par.T - 2
        for t in reversed(np.arange(par.T - 1)):
            # The value of an alternative is wage minus mobility costs + next period's continuation value discounted. 
            # Again, the transpose() and newaxis on EV makes sure that the choice dimension of w and M lines up with the slag dimension of EV.
            V_alternatives =  w[np.newaxis, :-1, t, :] - M[:, :-1, :] + par.rho * EV[:, 1:, t + 1].transpose()[np.newaxis, :, :]
            (EV[:, :-1, t], P[:, :-1, t, :]) = self.closed_forms(par.sigma, V_alternatives)

    def c_D_from_data(self, simulated_data=True):
        """ Calculates the variable D in the paper. The variable is defined over all state space points (Omega), 
            and measures how large a fraction of individuals are characterized by that state space point. 
            The variable sums to one within each time period by construction.
            When using this method with simulated data, c_simulated_data must be run first.
        """
        if simulated_data:
            d = gm.est.simulated_data
        else:
            raise Exception()
        mi = pd.MultiIndex.from_product([[s for s in np.arange(0, self.par.S)], list(self.par.ages), [t for t in np.arange(0, self.par.T)]], 
                                        names=["slag", "a", "t"])
        D = (pd.Series(0, index=mi).add((d.groupby(["slag", "a", "t"]).size() / d.groupby("t").size()), fill_value=0)
             .to_numpy()
             .reshape((self.par.S, self.par.N_ages, self.par.T)))
        self.est.D = D

    def simulate(self):
        """ Method to simulate. This is a wrapper method that chooses between various ways of simulating from the model."""
        if self.sim.endo_D:
            self.simulate_endoD()
        elif not self.sim.endo_D:
            self.simulate_exoD()

    def simulate_exoD(self):
        """ Simulate the model forward for a given set of skill prices and a given quantity of D (density) in each period. 
            Strictly speaking, this implies that the model contradicts itself, because it does not incorporate the number of 
            individuals choosing a specific action causing the true D to shifts. Instead, here we keep D fixed to its data target.
            Presumes that the method c_D_from_data() has been run beforehand. """
        self.sim.density = self.est.D[..., np.newaxis] * self.sol.P


    def simulate_endoD(self):
        """ Simulate the model forward for a given set of skill prices. D is endogenous, hence it gets updated according to actual choices in model. 
        We initialize the model in a point where all points of the state space are equally 
        populated. Then, we iterate forward in time and calculate the share of total employment accounted for by each point 
        in the state space. We do not have to draw gumbel shocks since we implicitly invoke a 'law of large numbers' and simply
        use the conditional choice probabilities as "transition" probabilities. """

        P = self.sol.P
        par = self.par

        #Measure how large a fraction of the people are at each point in the {state space x sector}. Symbol D in the paper + a choice dimension.
        density = np.zeros((par.S, par.N_ages, par.T, par.S))

        #This should be based on data later. Here we assume that all points in the state space are equally populated.
        #
        # init_share = 1 / (par.S * par.N_ages * par.S)
        density[..., 0, :] = self.sim.initDist[..., np.newaxis] * P[..., 0, :] #initial choice

        #Do something different than equally populated for the sake of illustration and debugging
        # density[:, :, 0, :] += np.array([- 1/2 * init_share, 0, 1/2 * init_share])[np.newaxis, np.newaxis, :]
        assert np.around(np.sum(density[:, :, 0, :]), 8) == 1

        #Specify the distribution across the state space for entering cohorts.
        #Eventually this could be based on the distribution of the cohort that entered the last year we have data (e.g. 2016)
        #Right now the only part of the state space that is relevant to specify is the lagged sector choice, since a = 30 by construction.
        #Just insert some non-uniform values
        # enter_share = 1 / par.S 
        # EnteringStateSpaceDist = np.zeros(shape=(par.S)) + enter_share + np.array([- 1/2 * enter_share, 0, 1/2 * enter_share])
        EnterDist = self.sim.EnterDist

        # assert np.sum(EnterDist, axis=0)

        #loop over t starting from t = 0 going until the second to last period (which inserts values into the last). We only replace values in density in years 1, 2, ... 
        #Since the starting year is data.
        # t = 1
        for t in np.arange(1, par.T):
            ## How many retired at the end the previous year?
            # The non-standard part is making sure people come into the model with age 30. 
            # np.sum(density[:, -1, :, t - 1]) sums over all parts of the state except age and time. Because we simply need ALL retiring people.
            retiring = np.sum(density[:, -1, t - 1, :])
            ## Entering cohort of age 30:
            density[:, 0, t, :] = retiring * EnterDist[:, t, np.newaxis] * P[:, 0, t, :]

            ## What will people of ages 31-65 do in this period.
            # The transpose moves the sector dimension in front of the age dimension. This is intuitive, since the current sector
            # becomes next period's lagged sector. And it needs to line up with the policy function which has dimensions
            # (s_lag, age, s_curr, time) 
            # We sum over period t-2 sectors (previous period's lagged sector) since they do not matter
            # for current choices. #We move previous' period's choice from the last index to the place of slag so it matches to P
            density[:, 1:, t, :] = np.sum(density[:, :-1, t-1, :], axis=0).swapaxes(0, -1)[:, :, np.newaxis] * P[:, 1:, t, :]

        assert all(np.around(np.sum(density, axis=(0, 1, 3)), 7) == 1), "Summing over state space (excluding time) and choices does not yield density == 1"
        self.sim.density = density

    def l_cohorts(self, simulated_data=True):
        """ Loads data on cohorts from the simulated data."""
        if simulated_data:
            d = self.est.simulated_data
        else:
            raise Exception("Can only do simulated data at the moment.")
        self.sim.EnterDist = d[d.a == self.par.a_min].groupby("t")["slag"].value_counts(normalize=True).to_numpy().reshape((self.par.T, self.par.S)).swapaxes(0, 1)

    def l_init_distribution(self, simulated_data=True):
        """ Loads data on initial (1996) distribution of the state space. This is used to initialize the model when simulating from 1996. """
        if simulated_data:
            d = self.est.simulated_data
        else:
            raise Exception("Can only do simulated data at the moment.")
        self.sim.initDist = d[d.t == 0][["slag", "a"]].value_counts(normalize=True).sort_index().to_numpy().reshape((self.par.S, self.par.N_ages))

    # def post_estimation(self):
    #     """This method calculates post-estimation quantities such as TFP (A), 
    #         physical capital (K), clean energy (E) and dirty energy (O)"""
    #     #This could change if we switch to the abatement form (JBE).

    def c_pars_to_estimate(self, parameters=None, indexes=None):
        """ Constructs an ordered dictionary where keys are parameter names
            and each value is an array of indexes. This indexes the parameter itself in self.par. 
            These two are controlled using the optional arguments 'parameters' and 'indexes'. 
            'indexes' should be either (1): an array (if all parameters to be estimated belong to the same
            attribute in par, or (2): if parameters from multiple attributes are to be estimated, 'indexes' 
            should be a list. The list can contain either arrays or None values. The value None implies that all
            parameters in the attribute will be estimated. """
        pars_to_estimate = OrderedDict()

        #default values
        if parameters is None:
            parameters = ["beta0"]
        if indexes is None:
            indexes = [None] * len(parameters)

        if isinstance(indexes, np.ndarray):
            #if only an array is given in 'indexes' it must refer to the index of only one parameter to be estimated
            assert len(parameters) == 1 
            pars_to_estimate[parameters[0]] = indexes
        else:
            for i, par in enumerate(parameters):
                if indexes[i] is None:
                    pars_to_estimate[par] = self.c_index_array(getattr(self.par, par))
                else:
                    assert isinstance(indexes[i], np.ndarray)
                    pars_to_estimate[par] = indexes[i]
        self.est.pars_to_estimate = pars_to_estimate
        #Follow-up helper objects
        self.c_n_params()
        self.c_gs()

    def c_n_params(self):
        n_params = {**{"total":0}, **{k:0 for k in self.par.groups.keys() if k is not "r"}}
        for p, index in self.est.pars_to_estimate.items():
            if index is None:
                n_params["total"] += 1
                n_params[get_key_from_value(self.par.groups, p)] += 1
            else:
                n_params["total"] += len(index)
                n_params[get_key_from_value(self.par.groups, p)] += len(index)
        self.est.n_params = n_params

    def c_gs(self):
        """Define the parameter groups (shorthand: gs) to be estimated. It is inferred directly from self.est.pars_to_estimate."""
        self.est.gs = set([g for g in self.par.groups.keys() for p in self.est.pars_to_estimate if p in self.par.groups[g]])

    def c_theta_idx(self):
        pte = self.est.pars_to_estimate

        n = 0
        theta_idx = OrderedDict()
        for parameter, index in pte.items():
            if index is not None:
                theta_idx[parameter] = np.arange(n, n + len(index))
                n += len(index)
            else:
                #scalar parameter vector
                theta_idx[parameter] = np.array([n])
                n += 1
        self.est.theta_idx = theta_idx

    def c_theta0(self):
        """ Constructs a single, 1-dimensional array containing all the parameters to be estimated.
            Requires self.pars_to_estimate to be defined (aliased pte here). I use the name theta0 as 
            the starting values of the parameters to be estimated. 
            The actual parameter values used are those stored in par when run. """
        pte = self.est.pars_to_estimate
        theta_idx = self.est.theta_idx
        #Initialize theta0 with the correct length.
        # theta0 = np.zeros(np.sum([len(pte[k]) if pte[k] is not None else 1 for k in pte]))
        # theta0 = np.zeros(max(theta_idx[next(reversed(theta_idx))]) + 1) #last element of theta_idx stores the
        theta0 = np.zeros(self.est.n_params["total"])
        # n = 0 #counter for index in theta0 itself

        for par, index in pte.items():
            if index is None:
                theta0[theta_idx[par]] = getattr(self.par, par)
                # n += 1
            else:
                theta0[theta_idx[par]] = getattr(self.par, par)[index]
                # n += len(index)
        self.est.theta0 = theta0

    @staticmethod
    def c_index_array(parameter):
        """ Constructs an index array from a parameter's values in self.par. If the variable is a scalar, 
            this method returns None, since a scalar cannot be indexed. 
            Examples: 
                np.array([0.5, 1.5, 2.5])   =>  np.array([0, 1, 2)]
                10                          =>  None
        """
        if isinstance(parameter, np.ndarray):
            return np.arange(parameter.size)
        elif isinstance(parameter, (float, int)):
            return None
        else:
            raise Exception("Parameter must be array or scalar")

    @staticmethod
    def update_par(par, theta, pars_to_estimate):
        """ For each of the values in stated in pars_to_estimate (and their indexes), the 
        corresponding parameter vectors in par are updated. The new values are stored in the 
        one-dimensional vector theta. """
        n = 0
        for key, index in pars_to_estimate.items():
            if index is not None:
                #array parameter
                vals = getattr(par, key) #previous value is retrieved
                vals[pars_to_estimate[key]] = theta[n : n + len(index)] #update
                setattr(par, key, vals) #new value is set
                n += len(index)
            else:
                #scalar parameter
                setattr(par, key, theta[n]) #new value is set
                n += 1

    def setup_estimation(self, parameters=None, indexes=None, method=None, agrad_loglik=True):
        """ This method controls a few different options for estimation.
            'partial' controls whether skill prices are kept fixed or endogenously determined in equilibrium.
            'method' controls the optimization algorithm. Common choice is 'BFGS'. 
            'agrad_loglik' controls whether to use numeric approximations for gradients or 
            the analytic gradients. """
        self.c_pars_to_estimate(parameters=parameters, indexes=indexes)
        self.c_theta_idx()
        self.c_theta0()
        if method is None:
            method = self.est.default_method
        self.c_est_options(method)
        self.est.agrad_loglik = agrad_loglik
        self.est.loglik_scale = len(self.est.simulated_data)

    def c_est_options(self, method):
        """Constructs a dictionary given to scipy.optimize.minimize when estimating. """
        options = {"disp":False, "maxiter":1000}
        if method == "BFGS":
            options["gtol"] = 1e-7
        self.est.options = options

    @staticmethod
    def est_ftol(method):
        """ Determines the tolerance for the objective function for estimation. This method allows it to
            vary depending on the chosen algorithm/method. """
        if method == "Nelder-Mead":
            return 1e-10
        else:
            return 1e-9

    def estimate(self, method=None):
        """ Estimate the parameters stored in self.est.pars_to_estimate. 
            Uses the starting values stored in self.est.theta0. 
        """
        if method is None:
            method = self.est.default_method
        #Prep estimation (assumes est.pars_to_estimate is already defined)
        theta0 = self.est.theta0
        pte = self.est.pars_to_estimate
        self.update_par(self.par, theta0, pte)
        self.c_switching_cost_matrix()

        print(f'objective function at starting values: {self.ll_objfunc(theta0)}')

        res = optimize.minimize(self.ll_objfunc, 
                                theta0, 
                                method=method, 
                                jac=self.est.agrad_loglik, 
                                tol=self.est_ftol(method), 
                                options=self.est.options)
        return res

    def solve_and_simulate(self):
        """ Runs the methods necessary to solve and simulate the model given parameters. """
        # self.precompute() #H changes with parameters and so needs to be precomputed again.
        self.solve_worker()
        if not self.partial_model:
            # self.solve_humancap_equilibrium()
            self.find_humcap_equilibrium()
        self.simulate()

    def ll_objfunc(self, theta):
        """ Update parameters, calculate equilibrium and evaluate loglik. If partial = True, 
            the skill prices (r) are kept fixed, and hence the human capital equilibrium conditions 
            are ignored. """
        # print(theta)
        self.update_par(self.par, theta, self.est.pars_to_estimate)
        self.c_switching_cost_matrix()
        self.precompute()
        self.solve_and_simulate()
        
        #data and CCPs 
        d = self.est.simulated_data
        P = self.sol.P
        ll = self.loglik(d, P)
        if self.est.agrad_loglik:
            return (ll, self.score())
        else:
            return ll

    def loglik(self, d, P):
        """Calculate the (negative of the) log likelihood for a sample d with choice probabilities c"""
        return - np.sum(np.log(P[d["slag"], d["a"] - self.par.a_min, d["t"], d["s"]])) / self.est.loglik_scale

    def score(self):
        """Calculate the score of the (negative of the) log likelihood. """
        #unpack
        # d = self.est.simulated_data #make something that chooses between simulated and actual datasets somehow.
        dll =  - self.gradients()
        return dll
        # return - np.sum(dlnP[d["slag"], d["a"] - self.par.a_min, d["t"], d["s"], :], axis=0) / len(d)

    def score_from_theta(self, theta):
        """ Helper function for self.score. This makes score a function of theta so that it can be used 
            when running functions such as scipy.optimize.check_grad(f, g, theta). 
        """
        #unpack
        self.update_par(self.par, theta, self.est.pars_to_estimate)
        self.c_switching_cost_matrix()
        self.precompute()
        self.solve_and_simulate()
        return self.score()

    def c_simulated_data(self, n_individuals=100_000):
        """ Simulate data from the model currently stored. We simulate N = 'n_individuals' individuals. 
            The total number of observations (NxT) is random because we simulate individuals of all ages.
            This means that some people only enter the sample 1 year because that is the year they turn 30.
            Other individuals only enter the sample in 3 years because they start at the age of 63. 
            To run this method, the worker's problem must already have been solved. It simulates data given
            the current parameter and skill price values stored implicitly in the choice probabilities (self.sol.P). 
        """
        self.rng = default_rng(123456)
        data = pd.DataFrame(-1, 
                            index=pd.MultiIndex.from_product([range(0, n_individuals), range(self.par.T)], names=["i", "t"]), 
                            columns=["s", "a", "slag"]).reset_index()[["i", "slag", "a", "t", "s"]] 
                            #i for individual, slag for lagged sector, a for age, t for year, s for sector (chosen) 

        #Randomly generate starting ages that mean some are younger than 30 in the first period, so they enter model later
        data.loc[data.t == 0, "a"] = self.rng.integers(self.par.a_min - self.par.T + 1, self.par.a_max + 1, n_individuals)
        #Update the ages through time deterministically
        for t in range(1, self.par.T):
            data.loc[(data.t == t), "a"] = (data.loc[data.t == t - 1, "a"] + 1).to_numpy()
        #Delete rows where workers are not active in the labor market anymore
        #For example, this deletes period 0 observations for those that are 29 or younger in that year
        data = data.loc[((data.a <= self.par.a_max) & (data.a >= self.par.a_min))]

        #Randomly (uniformly) generate period -1 sector choices (period 0 sector lagged)
        data.loc[data.groupby("i").apply(lambda x: min(x.index)).to_numpy(), "slag"] = self.rng.integers(0, self.par.S, n_individuals)

        for t in np.arange(0, self.par.T):
            #Update from previous period if we're not in the first one
            if t > 0:
                #Transformation function (For those that do not enter or exit, add lagged sector from previous choice)
                data.loc[((data.t == t) & (data.a > self.par.a_min)), "slag"] = data.loc[((data.t == t - 1) & (data.a < self.par.a_max)), "s"].to_numpy()
            for a in np.arange(self.par.a_min, self.par.a_max + 1):
                for slag in np.arange(0, self.par.S):
                    size = len(data.loc[((data.t == t) & (data.a == a) & (data.slag == slag)), "s"])
                    #If there are any individuals characterized by this state space point:
                    if size > 0:
                        #Randomly draw a chosen sector (based on choice probabilities i.e. indirectly the gumbel shocks)
                        data.loc[((data.t == t) & (data.a == a) & (data.slag == slag)), "s"] = \
                            self.rng.choice(self.par.S, size=size, p=self.sol.P[slag, a - self.par.a_min, t, :])
        self.est.simulated_data = data.reset_index(drop=True)

    def s_data(self, filename=None, simulated_data=True):
        """ Save the dataset currently stored in the model. This method makes it convenient to save and 
            later load a dataset rather than resimulating with each new instance of GM. 
        """
        if simulated_data:
            data = self.est.simulated_data
            if filename is None:
                filename = "simulated_data.pkl"
            elif not filename.endswith(".pkl"):
                filename += ".pkl"
        #Save
        data.to_pickle("../data/" + filename)

    def l_data(self, filename, simulated_data=True):
        """ Load a dataset stored as a .pkl file. """
        if not filename.endswith(".pkl"):
            filename += ".pkl"
        if simulated_data:
            self.est.simulated_data = pd.read_pickle("../data/" + filename)
        else:
            raise Exception()

    def plot_ll_3d(self, x_values, y_values):
        """ Plot the loglikelihood as a function of two parameters."""
        assert self.est.n_params["total"] == 2, "Exactly 2 parameters must be chosen for plotting."

        z = np.zeros((len(x_values)*(len(y_values))))

        n = 0
        for xval in x_values:
            for yval in y_values:
                z[n] = self.ll_objfunc(theta = np.array([xval, yval]))
                n += 1

        X, Y = np.meshgrid(x_values, y_values)
        Z = z.reshape(X.shape)

        fig = plt.figure()
        ax = plt.axes(projection="3d")
        ax.plot_surface(X, Y, Z, cmap ='viridis', edgecolor ='green')
        return (fig, ax)

    def partial_ll(self, dlnP):
        """ The partial derivative of the log likelihood can be calculated from dlnP/dtheta which is the input to this function. 
            To calculate this derivative, dlnP is evaluated in all the data points from the estimation sample. 
        """
        d = self.est.simulated_data
        return np.sum(dlnP[d["slag"], d["a"] - self.par.a_min, d["t"], d["s"], ...], axis=0) / self.est.loglik_scale

    def gradients(self):
        if self.partial_model:
            dlnP_dtheta = gradients.gradients(self.par, self.sol, self.est, self.partial_model).ll_gradients()
            return self.partial_ll(dlnP_dtheta)
        else:
            (dr_dtheta, dlnP_dtheta, dlnP_dr) = gradients.gradients(self.par, self.sol, self.est, self.partial_model).ll_gradients()
            dll_dtheta = self.partial_ll(dlnP_dtheta)
            dll_dr = self.partial_ll(dlnP_dr)
            return dll_dtheta + np.matmul(dll_dr.reshape((self.par.T*self.par.S), order="F"), dr_dtheta)


    def minimize_quadED(self):
        """ This function minimizes the objective function of the labor market equilibrium to find equilibrium skill prices. 
            Initial values for the skill prices are the ones currently stored in par.r. 
            Uses analytic gradients when self.agrad_quadED is True.
        """
        res = optimize.minimize(self.quadED, 
                                self.par.r.reshape((self.par.T * self.par.S), order="F"), 
                                method="BFGS", 
                                jac=self.agrad_quadED,
                                tol=1e-7,
                                options={"disp":False, "maxiter":1000, "gtol":1e-5})
        return res

    def c_ED(self):
        """ Calculates excess demand from the currently stored simulation density. """
        return self.par.alpha1 * self.par.pY / self.par.r - \
               np.sum(self.sim.density, axis=tuple(i for i in range(self.sim.density.ndim - 2))) * self.par.MASS[:, np.newaxis]

    def dED_dr(self, r=None):
        """ Calculates the derivative dED/dr. If a vector (1-D or 2-D) of skill prices is given, the model is solved and simulated
            for these skill prices before the derivative is calculated. If r is None, the derivative is calculated from the currently 
            stored solution and vector of prices.
        """
        if (r is not None):
            if r.ndim == 1:
                self.par.r = r.reshape(self.par.r.shape, order="F")
            elif r.ndim == 2:
                self.par.r = r
            self.precompute_w() #only wages need to be updated when the only thing we change is r
            self.solve_worker()
            self.simulate()
            
        return gradients.gradients(self.par, self.sol, self.est, self.partial_model).dED_dr() #Draws upon the gradient class in the gradient module

    def c_jac_quadED(self, r=None, ED=None):
        """ r1_id is the collection of skill prices as a 1-D array. 
            ED can be given directly so I avoid recomputing it."""
        # if r_1d is not None:
        #     self.par.r = r_1d.reshape(self.par.r.shape, order="F")
        #     self.solve_and_simulate()
        dED = self.dED_dr(r)
        if ED is None:
            ED = self.c_ED()
        # dED = self.partial_ED("r")
        #The 'ED / abs(ED)' enters to flip the sign of derivatives where ED < 0. What if ED is zero here? Could cause issues.
        # return np.sum(((ED < 0) * (-1) + (ED > 0) * 1)[..., np.newaxis, np.newaxis] * dED, axis=(0, 1)).reshape((self.par.T * self.par.S), order="F")
        # return np.sum((ED/np.abs(ED))[..., np.newaxis, np.newaxis] * dED, axis=(0, 1)).reshape((self.par.T * self.par.S), order="F")
        return np.sum(2 * ED[..., np.newaxis, np.newaxis] * dED, axis=(0, 1)).reshape((self.par.T * self.par.S), order="F")


    def ED_from_r(self, r_1d):
        """ Evaluates the objective function for finding the equilibrium on the labor market. 
            The function takes a 1-dimensional vector of skill prices. In this vector, index t moves the fastest, implying for example 
            that the first entries vary the year and keep the sector fixed. 
            The function resolves the model for this new vector of skill prices and evaluates whether the sum of absolute excess labor demands are zero.
            The function also return analytic gradients when self.agrad_quadED is True.
            REWRITE THIS DOCSTRING.
        """
        self.par.r = r_1d.reshape(self.par.r.shape, order="F") #Update rental prices
        self.precompute_w()
        self.solve_worker()
        self.simulate()
        #Calculate excess labor demand
        ED = self.c_ED()
        return ED.reshape((self.par.T * self.par.S), order="F")

    def solve_ED0(self):
        """ Solve the equation system ED = 0. The initial value for r is that stored in .par currently."""
        initvals = self.par.r.reshape((self.par.T * self.par.S), order="F").copy()
        if self.agrad_ED0:
            fprime = lambda x: self.dED_dr(x).swapaxes(0, 1).swapaxes(2, 3).reshape((self.par.T*self.par.S, 
                                                                                     self.par.T*self.par.S), order="C")
        else:
            fprime = None
        (x, infodict, ier, mesg) = optimize.fsolve(self.ED_from_r, initvals, fprime=fprime, full_output=True)
        assert ier == 1, "Solving for equilibrium not successful."
        #Replace r with the solution.
        self.par.r = x.reshape(self.par.r.shape, order="F")
        return (x, infodict, ier, mesg)

    def find_humcap_equilibrium(self):
        _ = self.minimize_quadED() #Solve using quadratic deviances.
        _ = self.solve_ED0() #Solve by settings ED=0 in equation system.

    def quadED(self, r_1d):
        """ Evaluates the objective function for finding the equilibrium on the labor market. 
            The function takes a 1-dimensional vector of skill prices. In this vector, index t moves the fastest, implying for example 
            that the first entries vary the year and keep the sector fixed. 
            The function resolves the model for this new vector of skill prices and evaluates whether the sum of absolute excess labor demands are zero.
            The function also return analytic gradients when self.sim.analytic_grad_ED is True.
            """
        self.par.r = r_1d.reshape(self.par.r.shape, order="F") #Update rental prices
        self.precompute_w() #With new skill prices, update wages
        self.solve_worker()
        self.simulate()
        # self.solve_worker() #solve the worker's problem
        # self.simulate()
        #Calculate excess labor demand
        ED = self.c_ED()
        # objfunc = np.sum(np.abs(ED)) #sum of absolute excess labor demands (scalar function)
        quadED = np.sum(np.square(ED))
        if self.agrad_quadED:
            #Calculate jacobian of objective function (sum(ED^2)) and reshape into 1-D
            return (quadED, self.c_jac_quadED(r=None, ED=ED))
        else:
            return quadED

    def avg_yearly_transition_rates(self, data=True):
        print(output.output(self.par, self.sol, self.sim, self.est, self.name, self.version).avg_yearly_transition_rates(data=data))

    def age_profile_switching(self, save=False):
        output.output(self.par, self.sol, self.sim, self.est, self.name, self.version).age_profile_switching(save=save)

    def time_profile_switching(self, save=False):
        output.output(self.par, self.sol, self.sim, self.est, self.name, self.version).time_profile_switching(save=save)

    def fig_employment_shares(self, save=False):
        output.output(self.par, self.sol, self.sim, self.est, self.name, self.version).fig_employment_shares(save=save)

    def fig_avg_wages(self, save=False):
        output.output(self.par, self.sol, self.sim, self.est, self.name, self.version).fig_avg_wages(save=save)

#%% Initiate model using the saved data from the testing cell below

gm = GM(endo_D=False)
gm.l_data("data_equilibrium_default_values.pkl")
gm.setup()
gm.allocate()
gm.precompute()
gm.l_cohorts()
gm.l_init_distribution()
gm.solve_worker()
gm.c_D_from_data()
gm.simulate()


# importlib.reload(gradients)
#%%

gm.estimate()
gm.c_simulated_data()


#%% Now, to develop new stuff! 


#%% Code for testing the general equilibrium gradients (quadED, ED0 and loglik) as well as simulate data after estimation.
#Should eventually be moved to a different script?

gm = GM(endo_D=True)
update_attrs(gm, agrad_quadED=False)
gm.setup(simple_statespace=False)
gm.allocate()
gm.precompute()
gm.solve_worker()
gm.c_simulated_data()
gm.l_init_distribution()
gm.l_cohorts()
#gm.c_D_from_data()
gm.simulate()
update_attrs(gm, partial_model=False)
gm.setup_estimation(method="BFGS", agrad_loglik=True)
# update_attrs(gm.sim, agrad_quadED=True)
res = gm.minimize_quadED() #find equilibrium without using analytic gradients (cannot since it needs D, which needs data.)
#Simulate new data using the current skill prices, which are in equilibrium.
gm.c_simulated_data()
update_attrs(gm.sim, endo_D=False)
gm.c_D_from_data()
update_attrs(gm, agrad_quadED=True)
update_attrs(gm.sim, endo_D=False)
#Find equilibrium while using gradients to check that it works.
res = gm.minimize_quadED()
# Find equilibrium skill prices using fsolve:
(x, infodict, ier, mesg) = gm.solve_ED0()
update_attrs(gm, agrad_quadED=False)
r_1d = gm.par.r.reshape((gm.par.T * gm.par.S), order="F").copy()

#Outer equilibrium loop gradient check
g = lambda x: gm.c_jac_quadED(x).reshape((gm.par.T * gm.par.S), order="F")
assert util.check_grad(r_1d, gm.quadED, g) < 1e-7

#Inner equilibrium loop gradient check
n = 0
for s in np.arange(gm.par.S):
    for t in np.arange(gm.par.T):
        f = lambda x: gm.ED_from_r(x)[n]
        g = lambda x: gm.dED_dr(x)[t, s, :, :].reshape((gm.par.T * gm.par.S), order="F")
        assert util.check_grad(x0=r_1d, f=f, jac=g) < 1e-7, "Inner solver gradient wrong"
        n += 1

#Loglik gradient checks (partial model and in general equilibrium)
gm.setup_estimation(parameters=["kappa1", "beta0", "kappa0", "xi_in", "xi_out"], indexes=None, agrad_loglik=False)
gm.partial_model = True

# loglik gradient in the partial model:
initvals = gm.est.theta0.copy()
assert util.check_grad(x0=initvals, f=gm.ll_objfunc, jac=gm.score_from_theta) < 1e-7

# loglik gradient in general equlibrium model
gm.partial_model = False
theta0 = gm.est.theta0.copy() 
res = optimize.approx_fprime(theta0, gm.ll_objfunc, 1.4901161193847656e-08)
gm.score_from_theta(theta0)
assert util.check_grad(theta0, gm.ll_objfunc, gm.score_from_theta) < 1e-7, "Loglikelihood does not work"

#Save data before moving parameters
gm.find_humcap_equilibrium()
gm.c_simulated_data()
gm.s_data("data_equilibrium_default_values")

update_attrs(gm.est, agrad_loglik=True)
gm.estimate()

#Tjek om vi kan estimere os tilbage til (tilpas tæt på) de sande parametre.
# gm.par.beta0[gm.est.pars_to_estimate["beta0"]] *= np.array([1.25, 0.75])
# gm.c_theta0()
# gm.estimate()



#%%

