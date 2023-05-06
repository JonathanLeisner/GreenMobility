"""  Defines the class that runs the Green Mobility (GM) project code."""
#%%
import time
import pickle
from types import SimpleNamespace
from collections import OrderedDict
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import optimize

#User-made modules
import util
import _gradients as gradients
# import output
import _basic
import _output

# import importlib
# importlib.reload(_output)

#%%
class GM_estimate(_basic.basic, _output.output_estimate):
    """ Class to solve the model presented in Leisner (2023). 
        Generally, the notation c_[name] for methods means construct/create/calculate an object or attribute called [name].
        The class relies on three helper modules:
        util    : Contains helper/utility functions for various simple functionality.
        gradient: Calculates gradients of the model.
        output  : Constructs figures and tables from the model.
    """
    def __init__(self, endo_D=None, name=None, from_pickle=False):
        if isinstance(from_pickle, str):
            with open("../model_instances/" + from_pickle + ".pickle", "rb") as handle:
                d = pickle.load(handle)
                for attr in ["sol", "par", "sim", "est", "partial_model"]:
                    setattr(self, attr, d[attr])
                # if "data_name" in d.keys():
                #     self.l_data(d["data_name"])
        else:
            self.sol = SimpleNamespace()
            self.par = SimpleNamespace()
            self.par.scale = SimpleNamespace()
            self.par.default = SimpleNamespace()
            self.sim = SimpleNamespace()
            self.est = SimpleNamespace()
            self.partial_model = False
            
        if name is None:
            if isinstance(from_pickle, str):
                self.name = from_pickle
            else:    
                self.name = "standard"
        else:
            self.name = name

        self.results = SimpleNamespace()
        self.results.tables = []
        self.results.figures = []
        # self.diag = SimpleNamespace()
        # self.diag.tables = SimpleNamespace()


        # self.version = "v4"
        # self.rng = default_rng(123456) #seed #need numpy version 1.17.00 or above to work
        
        if endo_D is None and not hasattr(self.sim, "endo_D"):        
            self.sim.endo_D = False #default value
        elif endo_D is not None:
            self.sim.endo_D = endo_D

        self.agrad_ED0 = True
        self.agrad_quadED = True

    def c_default_pars(self):
        for p in self.par.all_parameters:
            setattr(self.par.default, p, getattr(self.par, p).copy())

    def setup(self, version="standard", from_data=False, **kwargs):
        """ Initiates parameter values at their default values. Also populates meta settings. Named arguments 
            'keyword=value' of the function replace default values. 
        """
        self.setup_statespace(version)
        self.c_output_settings()
        self.l_calibrated(from_data)
        self.set_parameter_values(**kwargs)
        self.c_default_pars()
        self.set_parameter_restrictions()
        self.reset_skillprices()

    def setup_sectors(self, version="standard"):
        par = self.par
        if version == "standard":
            par.sector_names = ["Non-employment", "Food", "Manufacturing", "Energy", "Transport", "Services"]
        elif version == "simple":
            par.sector_names = ["Non-employment", "Food", "Manufacturing"]
        par.S = len(par.sector_names)

    def setup_periods(self, version="standard"):
        if version == "standard":
            self.par.T = 2016 - 2002 + 1
        elif version == "simple":
            self.par.T = 8

    def reset_skillprices(self):
        """
        Reset the skill prices (r) to default values.
        """
        par = self.par
        #Initial values for skill prices (r)
        par.r = np.linspace(np.arange(1, par.S), np.arange(1, par.S)[::-1], par.T, axis=0)

    def set_parameter_restrictions(self):
        """ 
        Set the parameter restrictions, i.e. the elements of various parameter vectors that will never be estimated.
        """
        self.par.restrictions = {"xi_out":np.array([1])}

    def l_calibrated(self, from_data):
        """
        Loads the data and parameters that are calibrated before estimation.
        """
        par = self.par
        #Fixed quantities from data
        par.MASS = self.l_MASS(from_data) #Number of workers. Scaled by average number of workers.
        par.pY1z = self.l_pY1z(from_data) #Real Value added (in million, per (average) workers)
        par.gYtau = self.l_gYtau(from_data) #Emission taxes paid
        par.alpha1 = self.l_alpha1(from_data) #Human capital share
        par.e = self.l_e()[:self.par.T, :self.par.S - 1]
        par.Z = self.l_Z()[:self.par.T, :self.par.S - 1]
        par.tau = self.l_tau()[:self.par.T, :self.par.S - 1]
        # par.psi = self.l_psi()[:self.par.T, :self.par.S - 1]
        par.rK = self.l_rK()[:self.par.T, :self.par.S - 1]
        par.z = self.l_z()[:self.par.S - 1]

    def l_z(self):
        return np.load("../data_prep/output/init_z.npy")

    def reset_par(self):
        """ Resets the parameters that can be estimated to their default values"""
        for p in self.par.all_parameters:
            setattr(self.par, p, getattr(self.par.default, p).copy())

    def l_MASS(self, from_data):
        """ Creates a time series of MASS."""
        if not from_data:
            MASS = np.ones(self.par.T)
            self.par.scale.MASS = 1
        else:
            MASS = np.load("../data_prep/output/MASS.npy")[:self.par.T]
            self.par.scale.MASS = int(MASS.mean())
        return MASS/self.par.scale.MASS

    def l_pY1z(self, from_data):
        """ Creates a time series of nominal output, i.e. p^Y times Y times (1-z) in the paper."""
        assert hasattr(self.par.scale, "MASS"), "Mass must be loaded first, so we can compute per-worker measures"
        if from_data:
            pY1z = np.load("../data_prep/output/pY1z.npy")[:self.par.T, :self.par.S - 1]
        else:
            pY1z = np.repeat((np.ones(self.par.S - 1))[np.newaxis, :], self.par.T, axis=0)
        return pY1z/self.par.scale.MASS

    def l_gYtau(self, from_data):
        """ Creates a time series of nominal output, i.e. p^Y times Y times (1-z) in the paper."""
        assert hasattr(self.par.scale, "MASS"), "Mass must be loaded first, so we can compute per-worker measures"
        if from_data:
            gYtau = np.load("../data_prep/output/gYtau.npy")[:self.par.T, :self.par.S - 1]
        else:
            gYtau = 0.01 * np.repeat((np.ones(self.par.S - 1))[np.newaxis, :], self.par.T, axis=0)
        return gYtau/self.par.scale.MASS

    def l_alpha1(self, from_data):
        """ Loads the alpha1 parameter (value added share of human capital) from the data """
        if from_data:
            return np.load("../data_prep/output/alpha1.npy")[:self.par.T, :self.par.S - 1]
        else:
            return np.repeat(np.array([0.5, 0.2, 0.6, 0.6, 0.6, 0.6])[:self.par.S - 1][np.newaxis, :], self.par.T, axis=0)

    def l_D_from_data(self):
        """ Calculates the variable D in the paper. The variable is defined over all state space points (Omega), 
            and measures how large a fraction of individuals are characterized by that state space point. 
            The variable sums to one within each time period by construction.
            When using this method with simulated data, c_simulated_data must be run first.
        """
        d = self.est.data

        mi = pd.MultiIndex.from_product([list(self.par.tenures), list(np.arange(0, self.par.S)), list(self.par.ages), list(np.arange(0, self.par.T))], 
                                        names=["ten", "slag", "a", "t"])
        D = (pd.Series(0, index=mi, dtype=float).add((d.groupby(["ten", "slag", "a", "t"]).size() / d.groupby("t").size()), fill_value=0)
                .to_numpy()
                .reshape((self.par.ten_max, self.par.S, self.par.N_ages, self.par.T)))
        self.est.D = D

    def simulate(self):
        """ Method to simulate. This is a wrapper method that chooses between various ways of simulating from the model."""
        if self.sim.endo_D:
            self.simulate_endoD()
        elif not self.sim.endo_D:
            self.simulate_exoD()

    def simulate_exoD(self):
        """ Simulate the model forward for a given set of skill prices and a given quantity of D (density) in each period. 
            Strictly speaking, this implies that the simulation contradicts itself, because it does not incorporate the number of 
            individuals choosing a specific action causing the true D to shifts. Instead, here we keep D fixed to its data target.
            Presumes that the method l_D_from_data() has been run beforehand. """
        self.sim.density = self.est.D[..., np.newaxis] * self.sol.P

    def simulate_endoD(self):
        """ Simulate the model forward for a given set of skill prices. D is endogenous, hence it gets updated according to actual choices in model. 
        We initialize the model in the initial state space distribution calculated from the data, using l_init_distribution(). 
        Then, we iterate forward in time and calculate the share of total employment accounted for by each point 
        in the state space. We do not have to draw gumbel shocks since we implicitly invoke a 'law of large numbers' and simply
        use the conditional choice probabilities as "transition" probabilities. 
        Note: This method assumes that the number of retiring people always equals next year's entering people, which is not literally true in the data.
        We do not use this method for actual estimation, so it does not affect estimation results. """

        P = self.sol.P
        par = self.par

        #Measure how large a fraction of the people are at each point in the {state space x sector}. Symbol D in the paper + a choice dimension.
        density = np.zeros_like(P)

        density[..., 0, :] = self.sim.initDist[..., np.newaxis] * P[..., 0, :] #initial choice
        assert np.isclose(np.sum(density[..., 0, :]), 1)

        #Specify the distribution across the state space for entering cohorts.
        enterDist = self.sim.enterDist

        #loop over t starting from t = 0 going until the second to last period (which inserts values into the last). 
        # We only replace values in density in years 1, 2, ... since the starting year is data.
        for t in np.arange(1, par.T):
            retiring = np.sum(density[..., -1, t - 1, :])
            ## Entering cohort of age 30:
            density[..., 0, t, :] = retiring * enterDist[..., t, np.newaxis] * P[..., 0, t, :]

            ## What will people of ages 31-65 do in this period.
            # The transpose moves the sector dimension in front of the age dimension. This is intuitive, since the current sector
            # becomes next period's lagged sector. And it needs to line up with the policy function which has dimensions
            # (s_lag, age, s_curr, time) 
            # We sum over period t-2 sectors (previous period's lagged sector) since they do not matter
            # for current choices. #We move previous' period's choice from the last index to the place of slag so it matches to P.
            # Tenure gets reset if the worker switches.
            #Switchers
            density[0, :, 1:, t, :] = np.sum(density[..., :-1, t-1, :] * (1 - np.eye(par.S)[..., np.newaxis, :]), axis=(-3, -4)).swapaxes(-1, -2)[..., np.newaxis] * \
                P[0, :, 1:, t, :]
            #Stayers
            stayers_today = np.sum(density[..., :-1, t-1, :] * np.eye(par.S)[..., np.newaxis, :], axis=-3)
            stayers_tomor = stayers_today[:-1, ...].copy() #the three groups with tenure below max yesterday
            stayers_tomor[-1, ...] += stayers_today[-1, ...] #those with maximum tenure yesterday keep maximum tenure
            density[1:, :, 1:, t, :] = stayers_tomor.swapaxes(-1, -2)[..., np.newaxis] * P[1:, :, 1:, t, :]

        self.sim.density = density
        assert all(np.isclose(np.sum(density, axis=tuple(s for s in np.arange(0, density.ndim - 2)) + (-1,)), 1)), "Summing over state space (excluding time) and choices does not yield density == 1"
        

    def l_enterDist(self):
        """ Loads data on cohorts from the currently stored data."""
        d = self.est.data
        enter = pd.Series(0, 
                          index=pd.MultiIndex.from_product([range(1, self.par.ten_max + 1), range(self.par.S), range(self.par.T)], names=["ten", "slag", "t"]))
        enter = enter.add(d[d.a == self.par.a_min].groupby(["ten", "slag", "t"]).size() / d[d.a == self.par.a_min].groupby("t").size(), fill_value=0)
        self.sim.enterDist = enter.to_numpy().reshape((self.par.ten_max, self.par.S, self.par.T))

    def l_init_distribution(self):
        """ Loads data on initial (t=0) distribution of the state space. This is used to initialize the model when simulating from the initial year. """
        d = self.est.data
        init = pd.Series(0, index=pd.MultiIndex.from_product([range(1, self.par.ten_max + 1), range(self.par.S), self.par.ages], names=["ten", "slag", "a"]))
        init = init.add(d[d.t == 0][["ten", "slag", "a"]].groupby(["ten", "slag", "a"]).size() / len(d[d.t == 0]), fill_value=0)
        self.sim.initDist = init.to_numpy().reshape((self.par.ten_max, self.par.S, self.par.N_ages))

    def post_estimation(self):
        """
        Calculate post-estimation quantities such as TFP (A), 
        physical capital (K)
        """
        self.est.terminal_K = (1 - self.par.alpha1[-1, :]) * (self.par.pY1z[-1, :] - self.par.gYtau[-1, :]) / self.par.rK[-1, :]
        self.est.terminal_A = self.par.pY1z[-1, :] / ((self.c_Hdem()[-1, :] ** self.par.alpha1[-1, :]) * 
                                                      (self.est.terminal_K ** (1 - self.par.alpha1[-1, :])) *
                                                      (1 - self.par.z)
                                                     )

    def c_indexes(self, parameters, indexes):
        """ 
        Constructs index objects for each parameter. An index object for a given
        parameter is an array which tells the program which elements in the parameter are supposed to be estimated.
        So for example, if we choose to estimate elements 2 and 3 of beta0, the function will return np.array([1, 2])
        for that specific parameter. The index arrays are organized in a list, with one array for each parameter in
        'parameters'. This function makes sure that if a None is given in indexes, it is replaced by an array which
        refers to all elements of the given parameter, with one exception: If the element is part of the list of
        parameter restrictions, i.e. is contained in par.restrictions, the element will not be referred to (and hence
        it is kept fixed).
        Example:
        Inputs: parameters = ["beta0"] and indexes = [None]. 
        Outputs: [np.array([0, 1, ... , len(beta0) - 1])].
        """
        new_indexes = []
        for p in parameters:
            old_index = indexes[parameters.index(p)]
            if old_index is None:
                old_index = self.c_index_array(getattr(self.par, p))
            if p in self.par.restrictions:
                new_index = old_index[old_index != self.par.restrictions[p]]
                assert len(new_index) > 0, f"You chose only non-estimatable parameters in {p}."
            else:
                new_index = old_index
            new_indexes.append(new_index)
        return new_indexes

    def c_pars_to_estimate(self, parameters=None, indexes=None):
        """ Constructs an ordered dictionary where keys are parameter names
            and each value is an array of indexes. This indexes the parameter itself in self.par. 
            These two are controlled using the optional arguments 'parameters' and 'indexes'. 
            'indexes' should be either (1): an array (if all parameters to be estimated belong to the same
            attribute in par, or (2): a list (if parameters from multiple attributes are to be estimated). 
            The list can contain either arrays or None values. The value None implies that all
            parameters in the attribute will be estimated. """
        
        #default values
        if parameters is None:
            parameters = self.par.all_parameters
        elif isinstance(parameters, str):
            parameters = [parameters]
        if indexes is None:
            indexes = [None] * len(parameters)
        elif isinstance(indexes, np.ndarray):
            assert len(parameters) == 1
            indexes = [indexes]

        self.est.groups_to_estimate = OrderedDict()
        for g, p_list in self.par.groups.items():
            pars = [p for p in p_list if p in parameters]
            if len(pars) > 0:
                self.est.groups_to_estimate[g] = pars

        pars_to_estimate = OrderedDict()            
        #looping through groups insures the ordering of groups is intact. Will never be util0, sigma, util1 etc.
        new_indexes = self.c_indexes(parameters, indexes)
        for _, pars in self.est.groups_to_estimate.items(): 
            for par in pars:
                pars_to_estimate[par] = new_indexes[parameters.index(par)]
        self.est.pars_to_estimate = pars_to_estimate

    def c_n_params(self):
        """ Calculates n_params, which is a dictionary with parameter groups as keys and the number of parameters 
            to be estimated within each group as values.
        """
        n_params = {**{"full":0}, **{k:0 for k in self.par.groups.keys() if k is not "r"}}
        for p, index in self.est.pars_to_estimate.items():
            if index is None:
                n_params["full"] += 1
                n_params[util.get_key_from_value(self.par.groups, p)] += 1
            else:
                n_params["full"] += len(index)
                n_params[util.get_key_from_value(self.par.groups, p)] += len(index)
        self.est.n_params = n_params

    def c_gs(self):
        """Define the parameter groups (shorthand: gs) to be estimated. It is inferred directly from self.est.pars_to_estimate."""
        self.est.gs = set([g for g in self.par.groups.keys() for p in self.est.pars_to_estimate if p in self.par.groups[g]])

    def c_theta_idx(self):
        """ Construct a two-layered dictionary. The upper level is a dictionary with parameters to be estimated as keys.
            The lower level of keys defines parameter groups g ('full', 'utility' or 'sigma') depending on the parameter. 
            The values are arrays that indicate these parameters' positions in the theta-vector (the 1-D vector of 
            parameters to be estimated) and the vector parameters in the corresponding subgroup. 
        """
        pte = self.est.pars_to_estimate

        #helper object for counting
        n = {**{"full":0}, **{k:0 for k in self.par.groups.keys()}}
        theta_idx = OrderedDict()
        for parameter, index in pte.items():
            g = util.get_key_from_value(self.par.groups, parameter)
            if index is not None:
                theta_idx[parameter] = {"full":np.arange(n["full"], n["full"] + len(index)), g:np.arange(n[g], n[g] + len(index))}
                n["full"] += len(index)
                n[g] += len(index)
            else:
                #scalar parameter vector
                theta_idx[parameter] = {"full":np.array([n["full"]]), g:np.array([n[g]])}
                n["full"] += 1
                n[g] += 1

        self.est.theta_idx = theta_idx

    def c_theta0(self, inplace=True, default=False):
        """ Constructs a single, 1-dimensional array containing all the initial values of the parameters to be estimated.
            Requires self.pars_to_estimate to be defined (aliased pte here). 
            The parameter values used are those stored in self.par when this method is run. 
        """
        pte = self.est.pars_to_estimate
        theta_idx = self.est.theta_idx
        theta0 = np.zeros(self.est.n_params["full"])
        if default:
            par = self.par.default
        else:
            par = self.par

        for p, index in pte.items():
            if index is None:
                theta0[theta_idx[p]["full"]] = getattr(par, p)
            else:
                theta0[theta_idx[p]["full"]] = getattr(par, p)[index]
        if inplace:
            self.est.theta0 = theta0
        else:
            return theta0

    @staticmethod
    def c_index_array(parameter):
        """ 
        Constructs an array of indexes from a parameter's values in self.par.
        Examples: 
            np.array([0.5, 1.5, 2.5])   =>  np.array([0, 1, 2)]
        """
        if isinstance(parameter, np.ndarray):
            return np.arange(parameter.size)
        else:
            raise Exception("Parameter must be array")

    @staticmethod
    def update_par(par, theta, pars_to_estimate):
        """ For each of the values in stated in pars_to_estimate (and their indexes), the 
        corresponding parameter vectors in par are updated. The new values are stored in the 
        one-dimensional vector theta which is given as an argument. """
        n = 0
        for key, index in pars_to_estimate.items():
            vals = getattr(par, key) #previous value is retrieved
            vals[pars_to_estimate[key]] = theta[n : n + len(index)] #update
            setattr(par, key, vals) #new value is set
            n += len(index)

    def setup_estimation(self, parameters=None, indexes=None, agrad_loglik=True):
        """ This method controls a few different options for estimation.
            'agrad_loglik' controls whether to use numeric approximations for gradients or 
            the analytic gradients. """
        self.c_pars_to_estimate(parameters=parameters, indexes=indexes)
        self.c_n_params()
        self.c_gs()
        self.c_theta_idx()
        self.c_theta0(inplace=True)
        self.c_est_options()
        self.est.agrad_loglik = agrad_loglik
        self.est.n_obs = len(self.est.data)

    def c_est_options(self):
        """Constructs a dictionary given to scipy.optimize.minimize when estimating. """
        self.est.options = {"disp":True, "maxiter":1000, "gtol":1e-8}
        self.est.ftol = 1e-9

    def estimate(self):
        """ Estimates the parameters stored in self.est.pars_to_estimate. 
            Uses the starting values stored in self.est.theta0. 
        """
        theta0 = self.est.theta0
        pte = self.est.pars_to_estimate
        self.update_par(self.par, theta0, pte)
        self.c_switching_cost_matrix()

        init_objfunc = self.ll_objfunc(theta0)
        if isinstance(init_objfunc, tuple):
            init_objfunc = init_objfunc[0]
        print(f'\nObjective function at starting values: {round(init_objfunc, 5)}')
        start_time = time.time()
        res = optimize.minimize(self.ll_objfunc, 
                                theta0, 
                                method="BFGS", 
                                jac=self.est.agrad_loglik, 
                                tol=self.est.ftol, 
                                options=self.est.options)
        self.est.theta = res.x
        self.est.SE = np.sqrt(np.diagonal(res.hess_inv/self.est.n_obs))

        print(f"""Estimation of {self.est.n_params['full']} parameter(s)""" + \
              f""" took a total of {round((time.time()-start_time)/60, 2)} minutes.""")
        self.post_estimation()
        return res

    def compare_init_to_estimated(self):
        """ 
        Constructs a dataframe that summarizes estimation results. The dataframe contains estimates,
        initial values, standard errors and t-statistics. It also contains the default values for parameters.
        This is especially interesting for simulated data, since in this case, these are also the true values.
        """
        tuples = []
        for p, indexes in self.est.pars_to_estimate.items():
            if indexes is None:
                tuples.append((p, " "))
            else:
                for idx in indexes:
                    tuples.append((p, idx))

        df = pd.DataFrame(index=pd.MultiIndex.from_tuples(tuples, names=("Parameter", "Index")), 
                          columns=["Default/True", "Initial", "Estimated", "SE", 
                                   "t-stat (p=0)", "t-stat (p=default)"])

        df.loc[:, "Default/True"] = self.c_theta0(inplace=False, default=True)
        df.loc[:, "Initial"] = self.est.theta0
        df.loc[:, "Estimated"] = self.c_theta0(inplace=False)
        if hasattr(self.est, "SE"):
            df.loc[:, "SE"] = self.est.SE
        else:
            df.loc[:, "SE"] = np.nan
        df.loc[:, "t-stat (p=0)"] = df.loc[:, "Estimated"]/df.loc[:, "SE"]
        df.loc[:, "t-stat (p=default)"] = (df.loc[:, "Estimated"] - df.loc[:, "Default/True"])/df.loc[:, "SE"]
        return df

    def compare_gradients_loglik(self):
        """ 
        Constructs a dataframe which compares gradients of the loglikelihood function. It compares
        gradients calculated by numerical approximation to those written analytically using the gradients module
        """
        tuples = []
        for p, indexes in self.est.pars_to_estimate.items():
            if indexes is None:
                tuples.append((p, " "))
            else:
                for idx in indexes:
                    tuples.append((p, idx))

        df = pd.DataFrame(index=pd.MultiIndex.from_tuples(tuples, names=("Parameter", "Index")), 
                          columns=["Numeric", "Analytic", "Difference", "Relative difference"])

        df["Analytic"] = self.score_from_theta(self.est.theta0.copy())
        self.reset_par()
        df["Numeric"] = optimize.approx_fprime(xk=self.est.theta0.copy(), f=self.ll_objfunc, epsilon=1.4901161193847656e-08)
        df["Difference"] = df["Numeric"] - df["Analytic"]
        df["Relative difference"] = (np.abs(df["Numeric"]) / np.abs(df["Analytic"])) - 1
        return df


    def solve_and_simulate(self):
        """ Runs the methods necessary to solve and simulate the model given parameters. 
            If self.partial_model is False, r is adjusted when this method is run too.
        """
        self.solve_worker()
        if not self.partial_model:
            self.find_humcap_equilibrium()
        self.simulate()

    def ll_objfunc(self, theta):
        """ Update parameters, calculate equilibrium and evaluate loglik. If partial = True, 
            the skill prices (r) are kept fixed, and hence the human capital equilibrium conditions 
            are ignored. 
        """
        self.update_par(self.par, theta, self.est.pars_to_estimate)
        # print(f"theta vector in model:{theta}")
        self.c_switching_cost_matrix()
        self.precompute()
        self.solve_and_simulate()
        # print(f"rental prices: {self.par.r}")
        # print(f"Lowest choice probability: {np.min(self.sol.P)}")
        
        #data and CCPs 
        ll = self.loglik(self.est.data, self.sol.P)
        if self.est.agrad_loglik:
            return (ll, self.score())
        else:
            return ll

    def loglik(self, d, P):
        """Calculate the (negative of the) log likelihood for a sample d with choice probabilities P"""
        return - np.sum(np.log(P[d["ten"] - 1, d["slag"], d["a"] - self.par.a_min, d["t"], d["s"]])) / self.est.n_obs

    def score(self):
        """Calculate the score of the (negative of the) log likelihood. """
        #unpack
        dll = - self.gradients()
        return dll
        # return - np.sum(dlnP[d["slag"], d["a"] - self.par.a_min, d["t"], d["s"], :], axis=0) / len(d)

    def score_from_theta(self, theta):
        """ Helper function for self.score. This makes score a function of theta so that it can be used 
            when running functions such as util.check_grad(f, g, theta). 
        """
        #unpack
        self.update_par(self.par, theta, self.est.pars_to_estimate)
        self.c_switching_cost_matrix()
        self.precompute()
        self.solve_and_simulate()
        return self.score()

    def setup_simulated_data(self, n_individuals=100_000):
        """ 
        Runs the methods that simulate data from the model and loads the data moments needed for simulation.
        This includes the initial state space distribution, the entering cohorts as well as the state space
        density in all years, symbol D. This is kept fixed throughout estimation. 
        """
        self.c_simulated_data(n_individuals=n_individuals)
        self.l_init_distribution()
        self.l_enterDist()
        # self.l_enterShare(from_data=False)
        self.l_D_from_data()

    def c_simulated_data(self, n_individuals=100_000):
        """ Simulate data from the model currently stored. We simulate N = 'n_individuals' individuals. 
            The total number of observations (NxT) is random because we simulate individuals of all ages.
            This means that some people only enter the sample 1 year because that is the year they turn 30.
            Other individuals only enter the sample in 3 years because they start at the age of 63. 
            To run this method, the worker's problem must already have been solved. It simulates data given
            the current parameter and skill price values stored implicitly in the choice probabilities (self.sol.P). 
        """
        # self.rng = default_rng(123456)
        start_time = time.time()
        self.rng = np.random
        data = pd.DataFrame(-1, 
                            index=pd.MultiIndex.from_product([range(0, n_individuals), range(self.par.T)], names=["i", "t"]), 
                            columns=["s", "a", "slag", "ten"]).reset_index()[["i", "ten", "slag", "a", "t", "s"]] 
                            #i for individual, slag for lagged sector, a for age, t for year, s for sector (chosen) 

        #Randomly generate starting ages that means some are younger than 30 in the first period, so they enter model later
        # data.loc[data.t == 0, "a"] = self.rng.integers(self.par.a_min - self.par.T + 1, self.par.a_max + 1, n_individuals)
        data.loc[data.t == 0, "a"] = self.rng.randint(self.par.a_min - self.par.T + 1, self.par.a_max + 1, n_individuals)
        #Update the ages through time deterministically
        for t in range(1, self.par.T):
            data.loc[(data.t == t), "a"] = (data.loc[data.t == t - 1, "a"] + 1).to_numpy()
        #Delete rows where workers are not active in the labor market anymore
        #For example, this deletes period 0 observations for those that are 29 or younger in that year
        data = data.loc[((data.a <= self.par.a_max) & (data.a >= self.par.a_min))].reset_index(drop=True)

        #Randomly (uniformly) generate period -1 sector choices (period 0 sector lagged)
        data.loc[data.groupby("i").apply(lambda x: min(x.index)).to_numpy(), "slag"] = self.rng.randint(0, self.par.S, n_individuals)

        #Randomly (uniformly) generate period 0 tenure
        data.loc[data.groupby("i").apply(lambda x: min(x.index)).to_numpy(), "ten"] = self.rng.randint(1, self.par.ten_max + 1, n_individuals)

        for t in np.arange(0, self.par.T):
            print(f"Simulating individuals. Currently year t={t}, the terminal year is {self.par.T - 1}.")
            #Update from previous period if we're not in the first one
            if t > 0:
                #Transformation function (For those that do not enter or exit, add lagged sector from previous choice)
                data.loc[((data.t == t) & (data.a > self.par.a_min)), "slag"] = data.loc[((data.t == t - 1) & (data.a < self.par.a_max)), "s"].to_numpy()
                #Update tenure (Those that stay get +1 tenure. Those that switch get tenure =1.)
                idx = ((data.s == data.slag) & (data.t == t - 1) & (data.a < self.par.a_max))
                data.loc[idx[idx].index + 1, "ten"] = (data.loc[idx, "ten"] + 1).clip(upper=self.par.ten_max).to_numpy()
                idx = ((data.s != data.slag) & (data.t == t - 1) & (data.a < self.par.a_max))
                data.loc[idx[idx].index + 1, "ten"] = 1

            for a in np.arange(self.par.a_min, self.par.a_max + 1):
                for slag in np.arange(0, self.par.S):
                    for ten in np.arange(1, self.par.ten_max + 1):
                        size = len(data.loc[((data.t == t) & (data.a == a) & (data.slag == slag) & (data.ten == ten)), "s"])
                        #If there are any individuals characterized by this state space point:
                        if size > 0:
                            #Randomly draw a chosen sector (based on choice probabilities i.e. indirectly the gumbel shocks)
                            data.loc[((data.t == t) & (data.a == a) & (data.slag == slag) & (data.ten == ten)), "s"] = self.rng.choice(self.par.S, size=size, p=self.sol.P[ten - 1, slag, a - self.par.a_min, t, :])
        self.est.data = data.reset_index(drop=True)
        print(f"""Simulation of {n_individuals} individuals """ + \
              f"""took a total of {round((time.time()-start_time)/60, 2)} minutes.""")


    def s_data(self, filename):
        """ Save the dataset currently stored in the model. This method makes it convenient to save and 
            later load a dataset rather than resimulating with each new instance of GM. 
        """
        data = self.est.data
        if not filename.endswith(".pkl"):
            filename += ".pkl"
        #Save
        data.to_pickle("../data/" + filename)

    def l_data(self, filename):
        """ Load a dataset stored as a .pkl file. """
        if not filename.endswith(".pkl"):
            filename += ".pkl"
        self.est.data = pd.read_pickle("../data/" + filename)
        self.l_init_distribution()
        self.l_enterDist()
        # self.l_enterShare(from_data=True)
        self.l_D_from_data()

    def plot_ll_3d(self, x_values, y_values):
        """ Plot the loglikelihood as a function of two parameters."""
        assert self.est.n_params["full"] == 2, "Exactly 2 parameters must be chosen for plotting."

        z = np.zeros((len(x_values)*(len(y_values))))

        n = 0
        for xval in x_values:
            for yval in y_values:
                z[n] = self.ll_objfunc(theta=np.array([xval, yval]))
                n += 1

        X, Y = np.meshgrid(x_values, y_values)
        Z = z.reshape(X.shape)

        fig = plt.figure()
        ax = plt.axes(projection="3d")
        ax.plot_surface(X, Y, Z, cmap='viridis', edgecolor='green')
        return (fig, ax)

    def partial_ll(self, dlnP):
        """ The partial derivative of the log likelihood can be calculated from dlnP/dtheta which is the input to this function. 
            To calculate this derivative, dlnP is evaluated in all the data points from the estimation sample. 
        """
        d = self.est.data
        return np.sum(dlnP[d["ten"] - 1, d["slag"], d["a"] - self.par.a_min, d["t"], d["s"], ...], axis=0) / self.est.n_obs

    def gradients(self):
        """ Calculates the gradients of the loglikelihood function, using the helper module 'gradients'."""
        if self.partial_model:
            dlnP_dtheta = gradients.gradients(self.par, self.sol, self.sim, self.est, self.partial_model).ll_gradients()
            return self.partial_ll(dlnP_dtheta)
        else:
            (dr_dtheta, dlnP_dtheta, dlnP_dr) = gradients.gradients(self.par, self.sol, self.sim, self.est, self.partial_model).ll_gradients()
            dll_dtheta = self.partial_ll(dlnP_dtheta)
            dll_dr = self.partial_ll(dlnP_dr)
            return dll_dtheta + np.matmul(dll_dr.reshape((self.par.T*(self.par.S - 1)), order="F"), dr_dtheta)

    def c_Hdem(self):
        """ Calculate human capital demand as given from the firm side. """
        return self.par.alpha1 / self.par.r * (self.par.pY1z - self.par.gYtau)

    def c_Hsup(self):
        """ Calculate human capital supply as given from the worker side."""
        return np.sum(self.sim.density[..., 1:] * self.sol.H[..., np.newaxis, :], 
                      axis=tuple(i for i in range(self.sim.density.ndim - 2))) * self.par.MASS[:, np.newaxis]


    def c_ED(self):
        """ Calculates excess demand from the currently stored simulation density. """
        return self.c_Hdem() - self.c_Hsup()

    def solve_ED0(self):
        """ Solve the equation system ED = 0. The initial value for r is that stored in .par currently."""
        initvals = self.par.r.reshape((self.par.T * (self.par.S - 1)), order="F").copy()
        if self.agrad_ED0:
            fprime = lambda x: self.dED_dr(x).swapaxes(0, 1).swapaxes(2, 3).reshape((self.par.T*(self.par.S - 1), 
                                                                                     self.par.T*(self.par.S - 1)), order="C")
        else:
            fprime = None
        (x, infodict, ier, mesg) = optimize.fsolve(self.ED_from_r, initvals, fprime=fprime, full_output=True)
        if ier != 1: 
            print("Solving for equilibrium not successful.")
        #Replace r with the solution.
        self.par.r = x.reshape(self.par.r.shape, order="F")
        return (x, infodict, ier, mesg)

    def successive_approximations(self, min_iterations=10, print_out=False):
        """ 
        Update the skill prices by successive approximation of the labor market equilibrium condition. 
        """
        if not hasattr(self.sim, "density"):
            self.precompute_w()
            self.solve_worker()
            self.simulate()
        max_iterations = 1000
        doubled_r_at_least_once = False
        for i in range(max_iterations):
            Hsup = self.c_Hsup()
            zeroHsup = np.isclose(Hsup, 0)
            Hsup[zeroHsup] = 1.0 #arbitrary number above 0 to avoid dividing by zero before replacing
            r1 = self.c_Hdem() * self.par.r / Hsup
            if np.sum(zeroHsup) > 0:
                #If supply is zero somewhere, double its price.
                r1[zeroHsup] = self.par.r[zeroHsup] * 2
                doubled_r_at_least_once = True
            err = np.sum(np.abs(r1 - self.par.r))
            if print_out:
                print(f"Sum of absolute deviances in skill prices is {err:.7f} after {i + 1} iterations.")
            if err > 100:
                stepsize = 0.05
            else:
                stepsize = 0.25     
            r2 = self.par.r * (1-stepsize) + stepsize * r1
            self.par.r = r2
            self.precompute_w()
            self.solve_worker()
            self.simulate()
            if (i >= min_iterations - 1) and (err < 0.005 * self.par.T * (self.par.S - 1)):
                if i > min_iterations - 1:
                    if print_out:
                        print(f"Sum of absolute deviances in skill prices was {err:.7f} after {i + 1} iterations.")
                        if doubled_r_at_least_once:
                            print("A human capital price was doubled at least once due to Hsup = 0.")
                break
        if i == max_iterations - 1:
            print(f"Successive approximations could not reach tolerance after {max_iterations} iterations.")
            print(f"Error ended at {err:.7f}.")
            if doubled_r_at_least_once:
                print("A human capital price was doubled at least once due to Hsup = 0.")

    def find_humcap_equilibrium(self):
        """ Find the prices r that set all excess labor demands equal to zero. 
            The equation system ED = 0 is solved."""
        self.successive_approximations()
        # _ = self.minimize_quadED()
        #assert res.success
        (_, _, ier, _) = self.solve_ED0()
        assert ier == 1, "Solving for equilibrium not succesful."

    def to_pickle(self):
        """
        Save the model object to a pickle. If save_data is True, it also saves the dataset currently stored in 
        the attribute self.est.data. Saving a model as a pickle allows new class instances to be instantiated from it.
        """
        # assert isinstance(save_data, bool), "Argument save_data must be a boolean."
        d = {"par":self.par, "sol":self.sol, "sim":self.sim, "est":self.est, "partial_model":self.partial_model}
        # if save_data:
        #     self.s_data(filename=self.name + "_data")
        #     d["data_name"] = self.name + "_data"
        with open("../model_instances/" + self.name + ".pickle", "wb") as handle:
            pickle.dump(d, handle, protocol=pickle.HIGHEST_PROTOCOL)

#%% Test new feautures of GM before implementing them in the class