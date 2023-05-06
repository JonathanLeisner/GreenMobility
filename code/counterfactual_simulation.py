"""  Defines the class that runs the Green Mobility (GM) project code."""
#%%
import os
os.chdir("I:/Workdata/706651/users/Jonathan/Projects/GreenMobility/code")
import importlib
import pickle
from types import SimpleNamespace
import numpy as np
from scipy import optimize
np.set_printoptions(precision=3)

import warnings
warnings.filterwarnings("error", category=RuntimeWarning)
warnings.filterwarnings("ignore", message="numpy.ufunc size has changed")
warnings.filterwarnings("ignore", message="numpy.ufunc size changed, may indicate binary incompatibility. Expected 192 from C header, got 216 from PyObject")
import matplotlib.pyplot as plt
import copy

#User-made modules
import _basic
import _output

#%% Counterfactually simulate

class counterfactual_simulation(_basic.basic, _output.output_counterfactual_simulation):
    """ 
    Class to run counterfactual simulation from 2016 onward. The class relies on methods from 
    the modules _basic and _output. 
    """
    def __init__(self, from_estimation_pickle=False):
        """ Initializes the containers sol, par, sim, est etc.
        if from_estimation_pickle is a string, the parameters and initial values are taken from the corresponding pickle.
        """
        self.sol = SimpleNamespace()
        self.par = SimpleNamespace()
        self.par.scale = SimpleNamespace()
        # self.par.default = SimpleNamespace()
        self.sim = SimpleNamespace()
        self.est = SimpleNamespace()
        self.results = SimpleNamespace()
        self.results.tables = []
        self.results.figures = []
        self.agrad_ED0 = False
        if isinstance(from_estimation_pickle, str):
            self.l_from_estimation(from_estimation_pickle)
            self.default_from_pickle = True
        else:
            self.default_from_pickle = False
        self.tax_scenarios = ["fixed tax", "rising tax food", "rising tax manufacturing"]
        self.price_scenarios = ["fixed price", "manufacturing price down"]
        self.current_tax_scenario = None
        self.current_price_scenario = None
        self.current_parameter_values_scenario = None
        self.comparison_scenario = SimpleNamespace() #container for results of the scenario for comparison (usually the fixed tax one)
        self.tau_scenario_description = ""


    def update_comparison_scenario(self):
        for attribute in ["par", "sol", "sim"]:
            setattr(self.comparison_scenario, attribute, copy.deepcopy(getattr(self, attribute)))

    def l_from_estimation(self, from_estimation_pickle):
        """
        Loads parameters and initial values from the pickle called 'from_estimation_pickle' and 
        saves them under the attribute 'self.from_estimation'.
        """
        with open("../model_instances/" + from_estimation_pickle + ".pickle", "rb") as handle:
            d = pickle.load(handle)
        self.from_estimation = SimpleNamespace()
        fe = self.from_estimation
        fe.par = SimpleNamespace()
        fe.par.scale = SimpleNamespace()
        fe.par.all_parameters = d["par"].all_parameters + ["rho"]
        for p in fe.par.all_parameters:
            setattr(fe.par, p, getattr(d["par"], p))
            if hasattr(d["par"].scale, p):
                setattr(fe.par.scale, p, getattr(d["par"].scale, p))
        fe.terminal_r = d["par"].r[-1, :] #time by s
        #change to enterDist instead of EnterDist
        # fe.enterDist = d["sim"].EnterDist[:, -1] #slag by time
        # data = d["est"].data
        # fe.initDist = ((data[data.t == data.t.max()][["slag", "a"]].groupby(["slag", "a"]).size() / len(data[data.t == data.t.max()]))
        #                .to_numpy().reshape((data["slag"].nunique(), data["a"].nunique())))
        fe.terminal_K = d["est"].terminal_K
        fe.terminal_A = d["est"].terminal_A
        fe.terminal_z = d["par"].z
        # fe.terminal_MASS = d["par"].MASS[-1]
        # fe.data = data

    def l_gamma(self):
        """ Loads gamma parameters for each sector. """ 
        self.par.gamma = np.load("../data_prep/output/gamma.npy")

    def c_g(self):
        """ Calculates the value of g given the currently stored value of z"""
        gamma = self.par.gamma
        z = self.sol.z
        self.sol.g = gamma[0, :] * (gamma[1, :] + gamma[2, :] * (1 - z) + gamma[3, :] * (1 - z) ** (1 + gamma[4, :]))

    def l_parameter_values(self):
        """ Load parameter values from self.from_estimation and assign them to self.par"""
        assert self.default_from_pickle
        for p in self.from_estimation.par.all_parameters + ["rho"]:
            setattr(self.par, p, getattr(self.from_estimation.par, p))
            if hasattr(self.from_estimation.par.scale, p):
                setattr(self.par.scale, p, getattr(self.from_estimation.par.scale, p))

    def c_parameter_values_scenario(self, scenario):
        """ Make an exogenous change in parameters. """
        if scenario == "low switching costs":
            self.par.xi_in *= 0.8
            self.c_switching_cost_matrix()

    def setup_sectors(self, version="standard"):
        """ Defines sector names as well as which sectors are tradeable and non-tradeable. """
        par = self.par
        if version == "standard":
            par.trade_names = ["Food", "Manufacturing"]
            par.nontrade_names = ["Energy", "Transport", "Services"]
        elif version == "simple":
            par.trade_names = ["Food"]
            par.nontrade_names = ["Services"]
        par.sector_names = ["Non-employment"] + par.trade_names + par.nontrade_names
        par.trade = np.array([par.sector_names.index(s) for s in par.trade_names])
        par.nontrade = np.array([par.sector_names.index(s) for s in par.nontrade_names])
        par.S = len(par.sector_names)

    def setup_periods(self, version="standard"):
        """ Defines the periods of simulation. Default is to simulate until 2060."""
        if version == "standard":
            self.par.T = 2060 - 2016 + 1
        elif version == "simple":
            self.par.T = 8

    def l_calibrated(self, from_data):
        """
        Loads the data and parameters that are calibrated before estimation.
        These do not rely on the estimated model but are values that are calibrated even before estimation.
        """
        par = self.par
        #Fixed quantities from data
        # par.MASS = self.l_MASS(from_data) #Number of workers. Scaled by average number of workers.
        par.alpha1 = self.l_alpha1(from_data) #Human capital share
        par.e = self.l_e()[-1, :]
        par.tau_init = self.l_tau()[-1, :]
        # par.psi = self.l_psi()[-1, :]
        par.rK = self.l_rK()[-1, :]
        par.mu = self.l_mu(from_data)

    def l_init_MASS(self, from_data=True):
        """ Creates a time series of MASS."""
        if not from_data:
            MASS = np.ones(1)
            self.par.scale.MASS = 1
        else:
            MASS = np.ones(1) * np.load("../data_prep/output/MASS.npy")[-1]
            self.par.scale.MASS = MASS.mean()
        return MASS/self.par.scale.MASS

    def c_future_tau(self, scenario, tau_factor_end=None):
        """ Constructs the time series of taxes. Different tax series developments are 
            chosen using the argument 'scenario'. 
        """
        assert scenario in self.tax_scenarios
        t_phase_in, t_phase_in_end = None, None
        self.current_tax_scenario = scenario
        if scenario == "fixed tax":
            self.par.tau = np.repeat(self.par.tau_init[np.newaxis, :], self.par.T, axis=0)
        elif scenario == "rising tax food":
            if tau_factor_end is None:
                tau_factor_end = 6
            self.par.tau = np.repeat(self.par.tau_init[np.newaxis, :], self.par.T, axis=0)
            self.par.tau[:, self.par.sector_names[1:].index("Food")] *= np.linspace(1, tau_factor_end, self.par.T)
        elif scenario == "rising tax manufacturing":
            if tau_factor_end is None:
                tau_factor_end = 6
            self.par.tau = np.repeat(self.par.tau_init[np.newaxis, :], self.par.T, axis=0)
            t_phase_in, t_phase_in_end = 2025, 2040
            self.par.tau[:, self.par.sector_names[1:].index("Manufacturing")] *= \
                np.concatenate((np.ones(t_phase_in - 2016 - 1), 
                                np.linspace(1, tau_factor_end, t_phase_in_end - t_phase_in + 2), 
                                tau_factor_end * np.ones(self.par.T - (t_phase_in_end - 2016) - 1))
                              )
        self.tau_scenario_description = self.tax_scenario_descriptions(scenario, tau_factor_end, t_phase_in, t_phase_in_end)

    def tax_scenario_descriptions(self, scenario, tau_factor_end=None, t_phase_in=None, t_phase_in_end=None):
        if scenario == "fixed tax":
            string = "The simulation assumes a constant tax rate in each sector in all years and represents the baseline scenario."
        elif scenario == "rising tax food":
            string = "The policy scenario assumes a constant tax rate in all sectors but the food sector." +\
               " The tax rate in the food sector rises linearly throughout the full time period" +\
               f" to reach a level {tau_factor_end} times the initial one."
        elif scenario == "rising tax manufacturing":
            string = "The policy scenario assumes a constant tax rate in all sectors but manufacturing." + \
                f" The tax rate in manufacturing rises linearly between {t_phase_in} and {t_phase_in_end}" +\
                f" to reach a level {tau_factor_end} times the initial one." 
        return string

    def setup_simulation(self, from_data=True):
        """ 
        Prepares a simulation by running other methods. 
        Loads the initial distribution of the state space, the distribution of entering cohorts, 
        the size of MASS over time and the share of the state space that entering cohorts constitute.
        """
        self.l_initDist(from_data)
        self.sim.enterDist = self.sim.initDist[..., 0] / np.sum(self.sim.initDist[..., 0])
        self.c_future_MASS(from_data) #also calculates enterShare

    def l_initDist(self, from_data=True):
        """ Loads the initial distribution of the state space from data."""
        if from_data:
            self.sim.initDist = np.load("../data_prep/output/terminal_Dist.npy")
        else:
            self.sim.initDist = 1 / (self.par.S * self.par.N_ages) * np.ones(shape=(self.par.S, self.par.N_ages))

    def c_future_MASS(self, from_data=True):
        """ Calculates the time series of MASS in the future under the assumption that entering generations are identical to those 
        that entered in 2016. Also calculated the share of the MASS that entering generations constitute in the future.
        Stores self.par.MASS and self.sim.enterShare. """
        MASS = np.zeros(self.par.T)
        MASS[0] = self.l_init_MASS(from_data)
        cohort_sizes = np.sum(self.sim.initDist, axis=tuple(i for i in range(self.sim.initDist.ndim - 1))) * MASS[0]
        for t in np.arange(1, self.par.N_ages):
            MASS[t] = MASS[t-1] + cohort_sizes[0] - cohort_sizes[self.par.N_ages - 1 - (t - 1)]
        for t in np.arange(self.par.N_ages, self.par.T):
            MASS[t] = MASS[t-1]
        self.sim.enterShare = cohort_sizes[0] / MASS
        self.par.MASS = MASS

    def l_alpha1(self, from_data):
        """ Loads the alpha1 parameter (value added share of human capital) from the data """
        if from_data:
            return np.load("../data_prep/output/alpha1.npy")[-1, :self.par.S - 1]
        else:
            return np.array([0.5, 0.2, 0.6, 0.6, 0.6, 0.6])[:self.par.S - 1]

    def l_mu(self, from_data):
        """ Loads the calibrated consumption share parameters, mu """
        if from_data:
            return np.load("../data_prep/output/mu.npy")
        else:
            return (1/(self.par.S - 1)) * np.ones(self.par.S - 1)

    def reset_skillprices(self):
        """
        Reset the skill prices (r) to default values. 
        If the model is initialized from a pickle, skill prices are reset to the values from that.
        """
        if self.default_from_pickle:
            self.par.r = np.repeat(self.from_estimation.terminal_r[np.newaxis, :], self.par.T, axis=0)        
        else:
            self.par.r = 0.05 * np.repeat(np.ones(self.par.S - 1)[np.newaxis, :], self.par.T, axis=0)

    def reset_trade_prices(self):
        """ Reset the prices of tradeables """
        self.sol.p_trade = np.repeat(np.ones(len(self.par.trade))[np.newaxis, :], self.par.T, axis=0)

    def c_future_trade_prices(self, scenario):
        assert scenario in self.price_scenarios
        self.current_price_scenario = scenario
        if scenario == "manufacturing price down":
            self.sol.p_trade = np.repeat(np.ones(len(self.par.trade))[np.newaxis, :], self.par.T, axis=0)
            self.sol.p_trade[:, self.par.sector_names[1:].index("Manufacturing")] = \
                np.concatenate([np.ones(5), 0.9 * np.ones(self.par.T - 5)])

    def reset_nontrade_prices(self):
        """ Resets the output prices of non-tradeables. """
        self.sol.p_nontrade = np.repeat(np.ones(len(self.par.nontrade))[np.newaxis, :], self.par.T, axis=0)
        self.update_p()

    def update_p(self):
        """ Update the output prices, self.sol.p by using the values of tradeable prices and non-tradeables prices
        currently stored. """
        self.sol.p = np.concatenate((self.sol.p_trade, self.sol.p_nontrade), axis=1)

    def setup(self, version="standard", from_data=True, **kwargs):
        """ 
        Runs methods that prepare statespace, parameter values, initial values etc.
        Initiates parameter values at their default values. Also populates meta settings. 
        """
        self.setup_statespace(version)
        self.c_output_settings()
        self.l_calibrated(from_data)
        if not self.default_from_pickle:
            self.set_parameter_values(**kwargs)
        else:
            self.l_parameter_values()
            self.c_switching_cost_matrix()
        self.reset_skillprices()
        self.reset_trade_prices()
        self.reset_nontrade_prices()
        self.update_p()
        self.reset_K()
        self.l_A()
        self.l_gamma()

    def reset_K(self):
        """ Reset the values of physical capital in self.sol.K, to their default values. """
        if self.default_from_pickle:
            self.sol.K = np.repeat(self.from_estimation.terminal_K[np.newaxis, :], self.par.T, axis=0)
        else:
            self.sol.K = np.repeat(np.array([0.215, 0.285, 0.316, 0.614, 3.006])[np.newaxis, :], self.par.T, axis=0)

    def l_A(self):
        """ Sets the (constant) values of TFP. If a pickle was loaded initially, it uses that values from that. """
        if self.default_from_pickle:
            self.par.A = self.from_estimation.terminal_A
        else:
            self.par.A = np.array([1.533, 2.214, 0.082, 2.929, 1.774])    

    def reset_z(self):
        """ Resets the values of the abatement shares z to their default values. 
        If a pickle was loaded initially, the values are reset to those."""
        if self.default_from_pickle:
            self.sol.z = np.repeat(self.from_estimation.terminal_z[np.newaxis, :], self.par.T, axis=0)
        else:
            self.sol.z = 0.01 * np.ones(shape=(self.par.T, self.par.S - 1))

    def simulate_endoD(self):
        """ Simulate the model forward for a given set of skill prices. D is endogenous, hence it gets updated according to actual choices in model. 
        We initialize the model in the initial state space distribution calculated from the data, using l_init_distribution(). 
        Then, we iterate forward in time and calculate the share of total employment accounted for by each point 
        in the state space. We do not have to draw gumbel shocks since we implicitly invoke a 'law of large numbers' and simply
        use the conditional choice probabilities as "transition" probabilities. """

        P = self.sol.P
        par = self.par

        #Measure how large a fraction of the people are at each point in the {state space x sector}. Symbol D in the paper + a choice dimension.
        density = np.zeros_like(P)

        density[..., 0, :] = self.sim.initDist[..., np.newaxis] * P[..., 0, :] #initial choice
        assert np.isclose(np.sum(density[..., 0, :]), 1)

        enterDist = self.sim.enterDist
        enterShare = self.sim.enterShare

        for t in np.arange(1, par.T):
            ## Entering cohort of age 30:
            density[..., 0, t, :] = enterShare[t] * enterDist[..., np.newaxis] * P[..., 0, t, :]

            # Tenure gets reset if the worker switches.
            #Switchers
            density[0, :, 1:, t, :] = np.sum(density[..., :-1, t-1, :] * (1 - np.eye(par.S)[..., np.newaxis, :]), axis=(-3, -4)).swapaxes(-1, -2)[..., np.newaxis] * \
                P[0, :, 1:, t, :] * (par.MASS[t - 1] / par.MASS[t])
            #Stayers
            stayers_today = np.sum(density[..., :-1, t-1, :] * np.eye(par.S)[..., np.newaxis, :], axis=-3)
            stayers_tomor = stayers_today[:-1, ...].copy() #the three groups with tenure below max yesterday
            stayers_tomor[-1, ...] += stayers_today[-1, ...] #those with maximum tenure yesterday keep maximum tenure
            density[1:, :, 1:, t, :] = stayers_tomor.swapaxes(-1, -2)[..., np.newaxis] * P[1:, :, 1:, t, :] * (par.MASS[t - 1] / par.MASS[t])

        self.sim.density = density
        assert all(np.isclose(np.sum(density, axis=tuple(s for s in np.arange(0, density.ndim - 2)) + (-1,)), 1)), "Summing over state space (excluding time) and choices does not yield density == 1"

    def allocate_firm(self):
        """  Initialize the values of Y1z (meaning Y * (1-z)) and z itself."""
        sol = self.sol
        par = self.par
        sol.Y1z = np.zeros((par.T, par.S - 1))
        self.reset_z()

    def c_Hdem_FOC_H(self):
        """ Calculate the values of human capital demand consistent with the first order condition for human capital."""
        self.sol.Hdem_FOC_H = ((1 / self.par.r) * self.par.A[np.newaxis, :] * self.par.alpha1[np.newaxis, :] *
                (self.sol.K ** (1 - self.par.alpha1[np.newaxis, :])) * \
                    (self.sol.p * (1 - self.sol.z) - self.sol.g * self.par.tau)) ** \
                (1 / (1 - self.par.alpha1[np.newaxis, :]))

    def c_Hdem_FOC_K(self):
        """ Calculate the values of human capital consistent with the first order condition for physical capital."""
        self.sol.Hdem_FOC_K = self.par.alpha1[np.newaxis, :] / (1 - self.par.alpha1[np.newaxis, :]) * \
                        self.par.rK[np.newaxis, :] / self.par.r * self.sol.K

    def c_FOC_z(self):
        """ Calculate the optimal value of z from the condition p/tau = -inv(g'(z))""" 
        gamma = self.par.gamma
        tau = self.par.tau
        p = self.sol.p
        self.sol.z = 1 - (((p / (tau * gamma[0, :][np.newaxis, :])) - gamma[2, :][np.newaxis, :]) * \
                          (1 / (gamma[3, :][np.newaxis, :] * (1 + gamma[4, :][np.newaxis, :])))) \
                         ** (1 / gamma[4, :][np.newaxis, :])

    def c_Hsup(self):
        """ Calculate human capital supply as given from the worker side."""
        self.sol.Hsup = np.sum(self.sim.density[..., 1:] * self.sol.H[..., np.newaxis, :], 
                               axis=tuple(i for i in range(self.sim.density.ndim - 2))) * self.par.MASS[:, np.newaxis]

    def c_ELD_FOC_H(self):
        """ Calculate the excess labor demand condition as given from the first order condition for human capital"""
        return self.sol.Hdem_FOC_H - self.sol.Hsup

    def c_ELD_FOC_K(self):
        """ Calculate the excess labor demand condition as given from the first order condition for physical capital"""
        return self.sol.Hdem_FOC_K - self.sol.Hsup

    def c_Y1z(self):
        """ Calculate the values of Y * (1-z) from the stored values of H, K and z """
        self.sol.Y1z = (1 - self.sol.z) * self.par.A[np.newaxis, :] * (self.sol.Hdem_FOC_K ** self.par.alpha1[np.newaxis, :]) * \
                 (self.sol.K ** (1 - self.par.alpha1[np.newaxis, :]))

    def c_EOD(self):
        """ Calculate excess output demand given the currently stored values."""
        return self.par.mu[np.newaxis, self.par.nontrade - 1] * np.sum(self.sol.p * self.sol.Y1z, axis=1)[:, np.newaxis] - \
               self.sol.p_nontrade * self.sol.Y1z[:, self.par.nontrade - 1]

    def c_ELD_and_EOD(self, endo_1d=None):
        """ Calculate and return all excess demand functions. Runs all the methods necessary to resolve and simulate the model."""
        if endo_1d is not None:
            self.update_endo(endo_1d)
        if np.any(np.isnan(self.par.r)) or np.any(self.par.r < 0) or np.any(np.isclose(self.par.r, 0)):
            # print(self.par.r)
            return np.ones(shape=(self.par.T * (self.par.S - 1) * 2 + self.par.T * len(self.par.nontrade))) * np.inf
        self.precompute_w()
        self.solve_worker()
        self.simulate_endoD()
        self.c_Hsup()
        self.c_FOC_z()
        self.c_g()
        self.c_Hdem_FOC_H()
        self.c_Hdem_FOC_K()
        self.c_Y1z()
        ELD_H = self.c_ELD_FOC_H()
        ELD_K = self.c_ELD_FOC_K()
        EOD = self.c_EOD()
        out = np.concatenate((ELD_H.reshape((self.par.T * (self.par.S - 1)), order="F"), 
                              EOD.reshape((self.par.T * len(self.par.nontrade)), order="F"),
                              ELD_K.reshape((self.par.T * (self.par.S - 1)), order="F")),
                             axis=0)
        if any(np.isnan(out)):
            raise Exception()
        return out 

    def update_endo(self, endo_1d):
        """ Update all variables being solved for by replacing their values with the ones in the 'endo_1d'."""
        self.par.r = endo_1d[:self.par.r.size].reshape(self.par.r.shape, order="F")
        self.sol.p_nontrade = endo_1d[self.par.r.size:self.par.r.size + self.sol.p_nontrade.size].reshape(self.sol.p_nontrade.shape, order="F")
        self.sol.K = endo_1d[self.par.r.size + self.sol.p_nontrade.size:].reshape(self.sol.K.shape, order="F")
        self.update_p()

    def endo_to_1d(self):
        """ Initialize the initial values of the variables being solved for.""" 
        init_r = self.par.r.reshape((self.par.T * (self.par.S - 1)), order="F")
        init_p = self.sol.p_nontrade.reshape((self.par.T * len(self.par.nontrade)), order="F")
        init_K = self.sol.K.reshape((self.par.T * (self.par.S - 1)), order="F")
        endo_1d = np.concatenate((init_r, init_p, init_K), axis=0)
        return endo_1d

    def solve_ED0(self):
        """ Solve the three-part equation system ELD = 0 (from FOC H), ELD = 0 (from FOC K) and EOD = 0. 
        The initial values for r, K and p^NT are those stored in .par and .sol currently."""
        initvals = self.endo_to_1d()
        (x, infodict, ier, mesg) = optimize.fsolve(self.c_ELD_and_EOD, initvals, fprime=None, full_output=True)
        if ier != 1: 
            print("Solving for equilibrium not successful.")
        #Replace r, K and p^NT with the solution.
        self.update_endo(x)
        return (x, infodict, ier, mesg)

#%%