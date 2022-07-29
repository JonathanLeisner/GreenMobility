"""  Defines the class that runs the Green Mobility (GM) project code."""
#%%

from types import SimpleNamespace
from collections import OrderedDict
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import optimize
from numpy.random import default_rng
from mpl_toolkits import mplot3d

np.set_printoptions(precision=3)
sns.set_theme()
import warnings
warnings.filterwarnings("error", category=RuntimeWarning)



#%%

class GM:
    """ Class to solve the model presented in Leisner (2023). 
        Generally, the notation c_[name] for methods means construct/create/calculate an object or attribute called [name].
    """
    def __init__(self, name = None, endo_D=True, analytic_grad_ED=True):
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

        self.resultspath = "../results/"
        self.version = "v3"
        self.rng = default_rng(123456) #seed

        self.est.default_method = "BFGS"
        self.est.partial = True
        self.setup_sim(endo_D=endo_D, analytic_grad_ED=analytic_grad_ED)

    def setup_sim(self, **kwargs):
        """ Sets up settings needed to simulate. 
            The keywords that can be used are:
            endo_D : {True/False}. Switching between endogenous and exogenous D (density).
            analytic_grad_ED : {True/False}. Switch between numeric jacobian and analytic jacobian for GE problem. """

        for k, v in kwargs.items():
            setattr(self.sim, k, v)


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

        #Helper object for calculating gradients of excess labor demand. Used e.g. in partial_ED 
        self.sim.ts_indexes = np.transpose([np.tile(np.arange(0, par.T), par.S), np.repeat(np.arange(0, par.S), par.T)])

    def setup(self, simple_statespace=False, **kwargs):
        """ Initiates parameter values at their default values. Also populates meta settings. Specify named
            'keyword=value' to the function to replace default values by those specified in 'value'. 
        """
        
        par = self.par
        par.groups = SimpleNamespace()
        self.setup_statespace(simple_statespace)
        
        #Meta-settings
        sol = self.sol
        sol.step_fraction = 0.05
        sol.maxiter = 50000
        sol.tolerance_r = 1e-10

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
        # par.beta0 = np.array([0.002, 0.006, 0.012, 0.008])[:par.S] #1-dimensional over sectors. 
        par.beta0 = np.array([0.2, 0.6, 1.2, 0.8])[:par.S] #1-dimensional over sectors. 
        par.scale.beta0 = 100 #Meaning true parameters are x times smaller than the values in the code

        #Switching costs (direct utility costs)
        if self.name == "high_switchingcost":
            par.xi_in = 20 * np.array([0.04, 0.05, 0.06, 0.07])[:par.S]
            par.xi_out = 20 * np.array([0.03, 0.06, 0.09, 0.07])[:par.S].transpose()    
        else:
            par.xi_in = np.array([0.04, 0.05, 0.06, 0.07])[:par.S]
            par.xi_out = np.array([0.03, 0.06, 0.09, 0.07])[:par.S].transpose()

        #Define parameter groups. These are used in self.gradient().
        par.groups.sigma = ["sigma"]
        par.groups.utility = ["beta0", "xi_in", "xi_out"]
        par.groups.prices = ["r"]

        # #Replace default values by those given explicitly to the method
        for k, v in kwargs.items():
            setattr(self.par, k, v)

        #Post-setup stuff (e.g. calculating M based on xi-parameters)
        self.c_switching_cost_matrix()

        self.diag.tables.skillprice_converge = \
            pd.DataFrame(np.nan, 
                         index=pd.MultiIndex.from_product([self.par.sector_names, range(self.par.T)], 
                                                          names=["Sector", "Year"]), 
                         columns=pd.MultiIndex.from_product([["r0", "r1"], range(0, self.sol.maxiter)], 
                                                            names=["Pre/post", "Iteration"]))

    def c_switching_cost_matrix(self):
        """ Takes the xi-parameters and constructs the matrix of switching costs."""
        par = self.par

        #Switching cost (M) has dimensions (s_{t-1} x s_{t}) i.e. (s_out x s_in)
        par.M = sum(np.meshgrid(par.xi_in, par.xi_out)) #add cost of going OUT of a sector with the cost of going IN to a sector.
        par.M = np.exp(par.M)
        np.fill_diagonal(par.M, 0) #Not switching sector is costless

    def c_human_capital_unit_prices(self):
        """ Human capital unit prices, the symbol r in the paper. This is just an initial value. 
            This variable is endogenously determined in equilibrium (when self.est.partial == False). 
        """
        par = self.par
        # SxT
        par.r = np.linspace(np.arange(1, par.S+1), np.arange(1, par.S+1)[::-1], par.T, axis=0)

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
        (EV[:, a, :], P[:, a, :, :]) = self.closed_forms(par.sigma, w[np.newaxis, a, :, :] - M[:, np.newaxis, :])        

        #PART I (age iteration (64 -> 30) within the terminal period, static expectation)
        t = par.T - 1
        for a in reversed(np.arange(par.N_ages - 1)):
            V_alternatives = w[np.newaxis, a, t, :] - M[:, :] + par.rho * EV[np.newaxis, :, a + 1, t]
            # My choice of sector today enters EV tomorrow as slag. 
            # Therefore the last dimension of w (the choice) must match the first dimension of EV (slag) 
            (EV[:, a, t], P[:, a, t, :]) = self.closed_forms(par.sigma, V_alternatives)

        #PART II (time iteration from T - 2 to 0.)
        for t in reversed(np.arange(par.T - 1)):
            # The value of an alternative is wage minus mobility costs + next period's continuation value discounted. 
            # Again, the transpose() and newaxis on EV makes sure that the choice dimension of w and M lines up with the slag dimension of EV.
            V_alternatives =  w[np.newaxis, :-1, t, :] - M[:, np.newaxis, :] + par.rho * EV[:, 1:, t + 1].transpose()[np.newaxis, :, :]
            (EV[:, :-1, t], P[:, :-1, t, :]) = self.closed_forms(par.sigma, V_alternatives)

    def c_D_from_data(self, simulated_data=True):
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

# np.sum(density, axis=(0, 1, 3))[all(np.around(np.sum(density, axis=(0, 1, 3)), 7) != 1)]


    def solve_humancap_equilibrium(self, print_out=True):
        """ Calculate the skill price consistent with human capital demand equalizing human capital supply.
            This presumes that the worker's problem has already been solved once and simulated for some value of wages/skill prices.
            Currently, the equilibrium prices are found by successive approximations. """
        idx = pd.IndexSlice
        df = self.diag.tables.skillprice_converge

        #Proposed wages from the first iteration
        r1 = self.par.alpha1 * self.par.pY / (np.sum(self.sim.density, axis=(0,1)) * self.par.MASS[:, np.newaxis])

        #Insert r0 and r1 for the current iteration
        iteration = 0
        df.loc[idx[:, :], idx[:, iteration]] = np.array([self.par.r.reshape(self.par.T * self.par.S), r1.reshape(self.par.T * self.par.S)]).T

        err_r = np.sum(np.abs(r1 - self.par.r))
        for iteration in range(1, self.sol.maxiter):
            if err_r < self.sol.tolerance_r:
                if print_out:
                    print(f"Number of iterations when converged was: {iteration}")
                    break
            else:
                if iteration == self.sol.maxiter - 1:
                    raise Exception(f"Human capital equilibrium could not be found after {self.sol.maxiter} iterations.")
                #print(f"Current error of skill prices is: {err_r:.7f}")
                #Make another iteration
                #Update the skill prices (and wages) then resolve and simulate 
                self.par.r = r1 * self.sol.step_fraction + self.par.r * (1 - self.sol.step_fraction)
                self.precompute_w() #H need not be recomputed within this inner loop, since it is the unaltered given some vector of parameters.
                self.solve_worker()
                self.simulate()
                #Proposed skill prices (S x T)
                r1 = self.par.alpha1 * self.par.pY / (np.sum(self.sim.density, axis=(0,1)) * self.par.MASS[:, np.newaxis])
                #Save the initial and proposed skill prices (for diagnostics plotting)
                df.loc[idx[:, :], idx[:, iteration]] = np.array([self.par.r.reshape(self.par.S * self.par.T), r1.reshape(self.par.S * self.par.T)]).T
                #Calculate deviation
                err_r = np.sum(np.abs(r1 - self.par.r))

    def l_cohorts(self, simulated_data=True):
        if simulated_data:
            d = self.est.simulated_data
        else:
            raise Exception("Can only do simulated data at the moment.")
        self.sim.EnterDist = d[d.a == self.par.a_min].groupby("t")["slag"].value_counts(normalize=True).to_numpy().reshape((self.par.T, self.par.S)).swapaxes(0, 1)

    def l_init_distribution(self, simulated_data=True):
        if simulated_data:
            d = self.est.simulated_data
        else:
            raise Exception("Can only do simulated data at the moment.")
        self.sim.initDist = d[d.t == 0][["slag", "a"]].value_counts(normalize=True).sort_index().to_numpy().reshape((self.par.S, self.par.N_ages))

    # def post_estimation(self):
    #     """This method calculates post-estimation quantities such as TFP (A), 
    #         physical capital (K), clean energy (E) and dirty energy (O)"""
    #     #This could change if we switch to the abatement form (JBE).
        


    def fig_skillprice_converge(self, save=False):
        """Figure showing the time series of unit skill prices (r) at various iterations from the human capital equilibrium function. """
        idx = pd.IndexSlice
        df = self.diag.tables.skillprice_converge
        #Remove nans (iterations that were not reached before convergence)
        df = df.loc[:, ~df.isna().all(axis=0)]

        fig = plt.figure(figsize=(5, 3.5), dpi=100)
        ax = fig.add_subplot(1, 1, 1)
        r = "r0"
        ax.plot(df.loc[idx["Unemployment", :], idx[r, 0]].unstack(level="Sector"), marker="o", color="green")
        ax.plot(df.loc[idx["Unemployment", :], idx[r, 5::20]].unstack(level="Sector"), color="orange", alpha=0.7)
        ax.plot(df.loc[idx["Unemployment", :], idx[r, df.columns.get_level_values(level="Iteration")[-1]]].unstack(level="Sector"), marker="x", color="red")
        #todo: skriv grøn: start etc. i en legend
        #todo: latex-skrift for r og t
        # ax.legend(loc='upper center', frameon=True, bbox_to_anchor=(0.5, -0.15), ncol=self.par.S)
        ax.set_xlabel("Time (t)")
        ax.set_title("Convergence of skill prices")
        ax.set_ylabel("Skill price (r)")

        if save:
            self.savefig(fig, "skillprice_converge")

    def savefig(self, fig, figname):
        """Save the figure object 'fig' as a .pdf with name 'figname'. """
        fig.savefig(self.resultspath + figname +  "_" + self.name + "_" + self.version + ".pdf", bbox_inches='tight')

    def fig_employment_shares(self, save=False):
        """ Creates a figure of employment shares in each sector, over time"""
        #Unpack
        emp = np.sum(self.sim.density, axis=tuple(i for i in range(self.sim.density.ndim - 2)))
        sector_names = self.par.sector_names

        fig = plt.figure(figsize=(5, 3.5), dpi=100)
        ax = fig.add_subplot(1, 1, 1)

        for s, name in enumerate(sector_names):
            ax.plot(emp[:, s], label=name)

        ax.legend(loc='upper center', frameon=True, bbox_to_anchor=(0.5, -0.15), ncol=self.par.S)
        ax.set_xlabel("Time")
        ax.set_ylabel("Employment share")
        ax.set_title("Employment shares over time")

        if save:
            self.savefig(fig, "employment_shares")

    def fig_avg_wages(self, save=False):
        """ Creates a figure of average wages across ages, over time"""
        #Unpack
        w = np.mean(self.sol.w, 0)
        sector_names = self.par.sector_names
        
        fig = plt.figure(figsize=(5, 3.5), dpi=100)
        ax = fig.add_subplot(1, 1, 1)

        for s, name in enumerate(sector_names):
            ax.plot(w[:, s], label=name)

        ax.legend(loc='upper center', frameon=True, bbox_to_anchor=(0.5, -0.15), ncol=self.par.S)
        ax.set_xlabel("Time")
        ax.set_title("Average wages")
        ax.set_ylabel("Average wage across ages")

        if save:
            self.savefig(fig, "average_wages")

    def uncond_switching_probs(self):
        """ Construct table of unconditional switching probabilities. This calculates only net-flows. Should be recomputed for gross-flows.
            e.g. such that with just two sectors, if workers are 50/50 dsitributed  and they all switch sectors, the probabilities do not become 0. """
        d = self.sim.density
        #todo unconditional switching probabilities.
        # #Do not include period 0, since the transition calculated from that means period -1 to period 0, which lies outside the model
        # probs = np.mean(np.sum(d[:, :, 1:, :], axis=1) / np.sum(d[:, :, 1:, :], axis=(1, 3))[:, np.newaxis, :], axis=2)
        # self.results.tables.uncond_switching_probs = pd.DataFrame(probs, index=self.par.sector_names, columns=self.par.sector_names)

    def c_pars_to_estimate(self, parameters=None, indexes=None):
        """ Constructs an ordered dictionary where keys are parameter names
            and each value is an array of indexes. This indexes the parameter itself in self.par. 
            These two are controlled using the optional arguments 'parameters' and 'indexes'. 
            'indexes' should be either an array (if all parameters to be estimated belong to the same
            attribute in par. If parameters from multiple attributes are to be estimated, 'indexes' 
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
        n_params = 0
        for _, index in self.est.pars_to_estimate.items():
            if index is None:
                n_params += 1
            else:
                n_params += len(index)
        self.est.n_params = n_params

    def c_gs(self):
        """Define the parameter groups (shorthand: gs) to be estimated. It is in inferred directly from self.est.pars_to_estimate."""
        self.est.gs = set([g for g in vars(self.par.groups).keys() for p in self.est.pars_to_estimate if p in vars(self.par.groups)[g]])

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
        theta0 = np.zeros(self.est.n_params)
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

    def setup_estimation(self, partial=False, method=None, analytic_gradients=None):
        """ This method controls a few different options for estimation.
            'partial' controls whether skill prices are kept fixed or endogenously determined in equilibrium.
            'method' controls the optimization algorithm. Common choice is 'BFGS'. 
            'analytic_gradients' controls whether to use numeric approximations for gradients or 
            the analytic gradients. """
        if method is None:
            method = self.est.default_method
        self.c_est_options(method)
        #Fix skill prices arbitrarily, or by equilibrium conditions
        if partial:
            self.est.partial = True
        else:
            self.est.partial = False
        self.est.analytic_gradients = analytic_gradients

    def c_est_options(self, method):
        """Constructs a dictionary given to scipy.optimize.minimize when estimating. """
        options = {"disp":True, "maxiter":1000}
        if method == "BFGS":
            options["gtol"] = 1e-10
        self.est.options = options

    @staticmethod
    def est_ftol(method):
        """ Determines the tolerance for the objective function for estimation. This method allows it to
            vary depending on the chosen algorithm/method. """
        if method == "Nelder-Mead":
            return 1e-10
        else:
            return 1e-10

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

        print(f'objective function at starting values: {self.obj_func(theta0)}')

        res = optimize.minimize(self.obj_func, 
                                theta0, 
                                method=method, 
                                jac=self.est.analytic_gradients, 
                                tol=self.est_ftol(method), 
                                options=self.est.options)
        return res

    def solve_and_simulate(self):
        """ Runs the methods necessary to solve and simulate the model given parameters. """
        self.precompute() #H changes with parameters and so needs to be precomputed again.
        self.solve_worker()
        if not self.est.partial:
            # self.solve_humancap_equilibrium()
            self.GE_humcap()
        self.simulate()

    def obj_func(self, theta):
        """ Update parameters, calculate equilibrium and evaluate loglik. If partial = True, 
            the skill prices (r) are kept fixed, and hence the human capital equilibrium conditions 
            are ignored. """
        print(theta)
        self.update_par(self.par, theta, self.est.pars_to_estimate)
        self.solve_and_simulate()
        
        #data and CCPs 
        d = self.est.simulated_data
        P = self.sol.P
        ll = self.loglik(d, P)
        if self.est.analytic_gradients:
            return (ll, self.score())
        else:
            return ll

    def loglik(self, d, P):
        """Calculate the (negative of the) log likelihood for a sample d with choice probabilities c"""
        return - np.sum(np.log(P[d["slag"], d["a"] - self.par.a_min, d["t"], d["s"]])) / len(d)

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
        self.solve_and_simulate()
        return self.score()

    def c_simulated_data(self, n_individuals=100_000):
        """ Simulate data from the model currently stored. We simulate N = 'n_individuals' individuals. 
            The total number of observations (NxT) is random because we simulate individuals of all ages.
            This means that some people only enter the sample 1 year because that is the year they turn 30.
            Other individuals only enter the sample in 3 years because they start at the age of 63. """
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
        if simulated_data:
            data = self.est.simulated_data
            if filename is None:
                filename = "simulated_data.pkl"
            elif not filename.endswith(".pkl"):
                filename += ".pkl"
        #Save
        data.to_pickle(filename)

    def l_data(self, filename, simulated_data=True):
        if not filename.endswith(".pkl"):
            filename += ".pkl"
        if simulated_data:
            self.est.simulated_data = pd.read_pickle(filename)

    def plot_ll_3d(self, x_values, y_values):
        """ Plot the loglikelihood as a function of two parameters."""
        assert self.est.n_params == 2, "Exactly 2 parameters must be chosen for plotting."

        z = np.zeros((len(x_values)*(len(y_values))))

        n = 0
        for xval in x_values:
            for yval in y_values:
                z[n] = self.obj_func(theta = np.array([xval, yval]))
                n += 1

        X, Y = np.meshgrid(x_values, y_values)
        Z = z.reshape(X.shape)

        fig = plt.figure()
        ax = plt.axes(projection="3d")
        ax.plot_surface(X, Y, Z, cmap ='viridis', edgecolor ='green')
        return (fig, ax)

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

    def dHsup(self, D, dP, g):
        """ Calculate the partial derivative of human capital suppyly with respect to either skill prices or theta. 
            The input g controls this. 
            Expressions for D and dP must be supplied as well. Implicitly, this function assumes that D is constant and so D' can be ignored.
        """
        ndim_add = {"r":(np.newaxis, np.newaxis), "utility":(np.newaxis,)}.get(g)
        dHsup =  self.par.MASS[(slice(None), np.newaxis) + ndim_add] * \
                        np.sum(D[(..., np.newaxis) + ndim_add] * dP, axis=tuple(np.arange(0, len(self.sol.EV.shape) - 1)))
        return dHsup

    def partial_ED(self, theta_or_r):
        """ Calculate the partial derivative of the excess labor demand functions (S x T of these) with respect to either theta or r.
            This partial derivative is used to find the equilibrium on the labor market but also to calculate the derivative dll/dtheta 
            when the model is not partial.
        """

        #unpack
        D = self.est.D
        par = self.par

        if theta_or_r == "theta":
            assert len(self.est.gs) == 1
            if "utility" in self.est.gs:
                du = self.du("utility")
                dv = self.dv_and_dEV(du, "utility")
                dP = self.dP(dv, "utility")
                dHsup = self.dHsup(D, dP, "utility")
                # dP =  P[..., np.newaxis]/par.sigma * (dv -  np.sum(P[..., np.newaxis] * dv, axis=3)[:, :, :, np.newaxis, ...])
                # dHsup =  par.MASS[(slice(None), np.newaxis) + (np.newaxis,)] * \
                #         np.sum(D[(..., np.newaxis) + (np.newaxis,)] * dP, axis=tuple(np.arange(0, len(self.sol.EV.shape) - 1)))
                dH = - dHsup #parameters cannot affect demand in the model so I simply leave it out

        elif theta_or_r == "r":
            
            ts_indexes = self.sim.ts_indexes
            
            du = self.du("r")
            dv = self.dv_and_dEV(du, "r")
            dP = self.dP(dv, "r")
            dHsup = self.dHsup(D, dP, "r")

            # dP =  P[..., np.newaxis, np.newaxis]/par.sigma * (dv -  np.sum(P[..., np.newaxis, np.newaxis] * dv, axis=3)[:, :, :, np.newaxis, ...])
            # dHsup =  par.MASS[(slice(None), np.newaxis) + (np.newaxis,)*2] * \
            #          np.sum(D[(..., np.newaxis) + (np.newaxis,)*2] * dP, axis=tuple(np.arange(0, len(self.sol.EV.shape) - 1)))

            dHdem = np.zeros((par.T, par.S, par.T, par.S))
            dHdem[ts_indexes[:, 0], ts_indexes[:, 1], ts_indexes[:, 0], ts_indexes[:, 1]] = np.reshape(par.alpha1 * par.pY / (- np.square(par.r)), (par.T*par.S), order="F")

            dH = dHdem - dHsup
        return dH

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
        if g == "r":
            #unpack
            H = self.sol.H
            ts_indexes = self.sim.ts_indexes

            #Allocate
            du = np.zeros((par.S, par.N_ages, par.T, par.S, par.T, par.S)) #does not depend on slag when we diff wrt. r
            
            #for testing, delete later:
            # du = np.arange(0, par.N_ages*par.T*par.S*par.T*par.S).reshape((par.N_ages, par.T, par.S, par.T, par.S), order="F")
             
            # ts_indexes = np.transpose([np.tile(np.arange(0, par.T), par.S), np.repeat(np.arange(0, par.S), par.T)]) 
            # #Object for picking combinations where st == s't'.
            #H has no t-dimension (given a), and is repeated for each t.
            du[:, :, ts_indexes[:, 0], ts_indexes[:, 1], ts_indexes[:, 0], ts_indexes[:, 1]] = np.repeat(H, par.T, axis=1)[np.newaxis, ...]  
        #note: dimensions of du vary depending on g
        return du

    def dlnP(self, theta_or_r):
        """ Calculates dlnP with respect to either the parameter vector (theta) or skill prices (r). The choice of this determined by the argument
            theta_or_r. The function calculates du, and from this dv, and then finally dlnP from dv. 
        """
        assert theta_or_r in ["theta", "r"]
        if theta_or_r == "theta":
            assert len(self.est.gs) == 1
            if "utility" in self.est.gs:
                #Assumes only utility for now
                du = self.du("utility")
                dv = self.dv_and_dEV(du, "utility") #så længe det er utility-parametre, kan denne køre på samme måde med alle parametrene (ligning 44)
                dlnP = 1/self.par.sigma * (dv - np.sum(self.sol.P[..., np.newaxis] * dv, axis=-2)[:, :, :, np.newaxis, :])
                
        if theta_or_r == "r":
            du = self.du("r")
            dv = self.dv_and_dEV(du, "r")
            dlnP = 1/self.par.sigma * (dv - np.sum(self.sol.P[..., np.newaxis, np.newaxis] * dv, axis=-3)[:, :, :, np.newaxis, :, :])
            
        return dlnP

    def dv_and_dEV(self, du, g):
        """ Calculates dv and dEV from du. g specifies whether the derivative is with respect to r or parameters, since this determines whether we have
            to add two dimensions or one (since r is defined over s and t while parameters are collected in one dimension). 
            The function only returns dv since dEV is never necessary to know on its own if we know dv.
            The calculation of derivatives follows the same backwards recursion structure as the solution to the worker's problem.
        """

        #unpack
        P = self.sol.P
        par = self.par
        
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

    def partial_ll(self, dlnP):
        """ The partial derivative of the log likelihood can be calculated from dlnP/dtheta which is the input to this function. 
            To calculate this derivative, dlnP is evaluated in all the data points from the estimation sample. 
        """
        d = self.est.simulated_data
        return np.sum(dlnP[d["slag"], d["a"] - self.par.a_min, d["t"], d["s"], ...], axis=0) / len(d)

    def gradients(self):
        """ Calculates analytic gradients of the log likelihood function. 
            When the model is solves in partial mode, this assumes dr/dtheta = 0 and so only the partial effect dll/dtheta remains.
        """
        if self.est.partial:
            return self.partial_ll(self.dlnP("theta"))
        else:
            T = self.par.T
            S = self.par.S
            #equation with label dll_dtheta in the paper
            dED_dr_inv = np.linalg.inv(self.partial_ED("r").swapaxes(0, 1).swapaxes(2, 3).reshape((T*S, T*S), order="C")) 
            dr_dtheta = np.matmul(dED_dr_inv, - self.partial_ED("theta").reshape((T*S, self.est.n_params), order="F"))
            dll_dtheta = self.partial_ll(self.dlnP("theta"))
            dll_dr = self.partial_ll(self.dlnP("r"))
            
            #1-dimensional array with shape (nparams,)
            return dll_dtheta + np.matmul(dll_dr.reshape((T*S), order="F"), dr_dtheta)
        
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

    def GE_humcap(self):
        """ This function minimizes the objective function of the labor market equilibrium to find equilibrium skill prices. 
            Initial values for the skill prices are the ones currently stored in par.r. 
            Uses analytic gradients when self.sim.analytic_grad_ED is True.
        """
        res = optimize.minimize(self.objfunc_ED, 
                                self.par.r.reshape((self.par.T * self.par.S), order="F"), 
                                method="BFGS", 
                                jac=self.sim.analytic_grad_ED,
                                tol=1e-9,
                                options={"disp":True, "maxiter":1000})
        return res

    def c_ED(self):
        """ Calculates excess demand from the currently stored simulation density. """
        return self.par.alpha1 * self.par.pY / self.par.r - \
               np.sum(self.sim.density, axis=tuple(i for i in range(self.sim.density.ndim - 2))) * self.par.MASS[:, np.newaxis]

    def c_jacob_objfunc_ED(self, r_1d=None, ED=None):
        """ r1_id is the collection of skill prices as a 1-D array. 
            ED can be given directly so I avoid recomputing it."""
        if r_1d is not None:
            self.par.r = r_1d.reshape(self.par.r.shape, order="F")
            self.solve_and_simulate()
        if ED is None:
            ED = self.c_ED()
        dH = self.partial_ED("r")
        #The 'times 2' enters because it is a quadratic loss function being minimized. 
        return np.sum(2 * ED[..., np.newaxis, np.newaxis] * dH, axis=(0, 1)).reshape((self.par.T * self.par.S), order="F")

    def objfunc_ED(self, r_1d):
        """ Evaluates the objective function for finding the equilibrium on the labor market. 
            The function takes a 1-dimensional vector of skill prices. In this vector, index t moves the fastest, implying for example 
            that the first entries vary the year and keep the sector fixed. 
            The function resolves the model for this new vector of skill prices and evaluates whether the sum of absolute excess labor demands are zero.
            The function also return analytic gradients when self.sim.analytic_grad_ED is True.
            """
        self.par.r = r_1d.reshape(self.par.r.shape, order="F") #Update rental prices
        self.precompute_w() #With new skill prices, update wages
        self.solve_worker() #solve the worker's problem
        self.simulate()
        #Calculate excess labor demand
        ED = self.c_ED()
        objfunc = np.sum(np.square(ED)) #sum of squared excess labor demands (scalar function)
        if self.sim.analytic_grad_ED:
            #Calculate jacobian of objective function (sum(ED^2)) and reshape into 1-D
            return (objfunc, self.c_jacob_objfunc_ED(r_1d=None, ED=ED))
        else:
            return objfunc

#%%
self = gm
self.par.r.reshape((self.par.T * self.par.S), order="F")



r_1d.reshape(self.par.r.shape, order="F")

gm = GM(endo_D=True, analytic_grad_ED=True)
gm.setup(simple_statespace=False)
gm.c_human_capital_unit_prices()
gm.allocate()
gm.precompute()
gm.solve_worker()

# gm.s_data(filename="simdata", simulated_data=True)
# gm.l_data("simdata", simulated_data=True)


# Få simuleret data i en ligevægt og se om vi kan matche gradienter i ligevægten. Lige nu passer de slet ikke. 
gm.l_cohorts()
gm.l_init_distribution()


gm.simulate()

gm.GE_humcap()
gm.c_simulated_data()
gm.c_D_from_data()

#%%

gm.c_pars_to_estimate(parameters=["beta0"], indexes=np.array([2, 1]))
gm.c_theta0()
gm.c_theta_idx()
gm.setup_estimation(method="BFGS", partial=False, analytic_gradients=True)

#%% test Loglik med GE

optimize.check_grad(gm.obj_func, gm.score_from_theta, gm.par.beta0[gm.est.pars_to_estimate["beta0"]])

g = gm.gradients()

gm.est.analytic_gradients = False

gm.score_from_theta(gm.par.beta0[gm.est.pars_to_estimate["beta0"]])

# assert 

optimize.check_grad(gm.obj_func, gm.score_from_theta, gm.par.beta0[gm.est.pars_to_estimate["beta0"]])


#  < 1e-7, \
#        "Differences between analytic and numeric gradients are too large"


#%%
dED_theta = gm.partial_ED("theta")
dED_r = gm.partial_ED("r")

gm.par.r.reshape((gm.par.T * gm.par.S), order="F")


#order F means the first index, time, switches faster. So the first values loop through t' keeping s' fixed before switching to the next s'
# ED.reshape((gm.par.T * gm.par.S), order="F") 

#A mimicks dED_r
# A = np.arange(0, dED_r.size).reshape((S, T, S, T), order="C").swapaxes(0, 1).swapaxes(2, 3)
# B = A.swapaxes(0, 1).swapaxes(2, 3).reshape((T*S, T*S), order="C")

# dED_dr_inv = np.linalg.inv(self.partial_ED("r").swapaxes(0, 1).swapaxes(2, 3).reshape((T*S, T*S), order="C")) 
# dED_dr_inv = np.linalg.inv(dED_r.swapaxes(0, 1).swapaxes(2, 3).reshape((T*S, T*S), order="C"))

# C = np.arange(dED_theta.size).reshape((T*S, self.est.n_params), order="F")

# b = np.arange(36).reshape((6, 6), order="F")
# c = np.arange(12).reshape((6, 2), order="F")
# A = np.arange(0, dED_dr.size).reshape(dED_dr.shape, order="C")

# #t, s, t', s'
# A[0, 0, 0, 0]

# dED_dr.shape 

dED_theta_2d = - dED_theta.reshape((T*S, self.est.n_params), order="F")

dlnP_dr = gm.dlnP("r")

dr_dtheta = np.matmul(dED_dr_inv, dED_theta_2d)

dr_dtheta.shape

dlnP_dr.shape

d = self.est.simulated_data
# dlnP = self.gradients()
dL_dr = - np.sum(dlnP_dr[d["slag"], d["a"] - self.par.a_min, d["t"], d["s"], ...], axis=0) / len(d)

led2 = np.matmul(dL_dr.reshape((T*S), order="F"), dr_dtheta)
dL_dtheta = np.sum(gm.dlnP("theta")[d["slag"], d["a"] - self.par.a_min, d["t"], d["s"], ...], axis=0) / len(d)

dL_dtheta + led2



# led2.shape

# dED_dr_inv.shape
# dED_theta_2d.shape

# dED_theta[0, 0, 0]
# dED_theta[1, 0, 0]
# dED_theta_2d[0, 0]
# dED_theta_2d[1, 0]

# dED_dr_inv.shape

# x = gm.partial_ED("r")
# y = np.load("grad_r_old.npy")


# gm.c_theta0()
# gm.score_from_theta(self.par.beta0[gm.est.pars_to_estimate["beta0"]])

# res = gm.GE_humcap()

#
# 
# new = gm.c_jacob_objfunc_ED(r_1d).copy()

# r_1d *= 4
# optimize.check_grad(gm.objfunc_ED, gm.c_jacob_objfunc_ED, r_1d)
#%%
gm.setup_sim(analytic_grad_ED=False)
r_1d = gm.par.r.reshape((gm.par.T * gm.par.S), order="F")
assert optimize.check_grad(gm.objfunc_ED, gm.c_jacob_objfunc_ED, r_1d) < 1e-7, \
       "Differences between analytic and numeric gradients are too large"

#%%



# %% Run

gm = GM()
gm.setup(simple_statespace=False)
gm.c_human_capital_unit_prices()
gm.allocate()
gm.precompute()
gm.solve_worker()
gm.simulate()
gm.solve_humancap_equilibrium()
r_old = gm.par.r.copy()
# gm.par.r = gm.par.r * 1.5
res = gm.GE_humcap()

res.x.reshape(r_old.shape) - r_old
np.isclose(res.x.reshape(r_old.shape), r_old)

#%%
#Perturb and try again
gm.par.beta0[np.array([1, 2])] *= np.array([0.5, 1.5])
res = gm.GE_humcap()

res

#%%

#%%
# gm.c_simulated_data(n_individuals=50_000)
gm.setup_estimation(method="BFGS", partial=True, analytic_gradients=True)
gm.c_pars_to_estimate(parameters=["beta0"], indexes=np.array([1, 2]))
gm.c_theta_idx()
gm.c_theta0()

#%%


#%%

# np.save("P_old", gm.sol.P)
# np.save("density_old", gm.sim.density)
# np.save("score_old", gm.score())

#%%


#%% Compare analytic and numeric gradients

gm.est.analytic_gradients = False


assert optimize.check_grad(gm.obj_func, gm.score_from_theta, gm.par.beta0[gm.est.pars_to_estimate["beta0"]]) < 1e-7, \
       "Differences between analytic and numeric gradients are too large"
