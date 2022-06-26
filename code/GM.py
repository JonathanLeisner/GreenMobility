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
    def __init__(self, name = None):
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

        self.est.default_method = "Nelder-Mead"

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
        par.sigma = 0.8

        #Should be loaded from data later
        par.alpha1 = np.repeat(np.array([0.2, 0.2, 0.3, 0.1])[:par.S][:, np.newaxis], par.T, axis=1)
        par.alpha2 = np.repeat(np.array([0.2, 0.2, 0.2, 0.1])[:par.S][:, np.newaxis], par.T, axis=1)
        par.alpha3 = 1 - par.alpha1 - par.alpha2

        #Parameters to estimate later
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
        par.r = np.linspace(np.arange(1, par.S+1), np.arange(1, par.S+1)[::-1], par.T, axis=1)

    def load_MASS(self, version="constant"):
        """ Creates a time series of MASS."""
        #for testing, this makes the population double halfway through the sample
        if version == "constant":
            return np.ones(self.par.T)
        else:
            return np.concatenate([np.ones(int(round(self.par.T/2, 0))), np.ones(self.par.T - int(round(self.par.T/2, 0))) + 1])

    def load_nominal_output(self):
        """ Creates a time series of nominal output, i.e. p^Y times Y in the paper."""
        return np.repeat((np.ones(self.par.S))[:, np.newaxis], self.par.T, axis=1)

    def allocate(self):
        """ Allocate empty containers for solution objects. 
        c contains the CCPs and v contains expected value functions"""
        sol = self.sol
        par = self.par

        sol.c = np.zeros((par.S, par.N_ages, par.S, par.T)) - 1 #4 dimensions: slag, a, s (choice), t
        sol.v = np.zeros((par.S, par.N_ages, par.T)) - 99 #3 dimensions: slag, a, t

    def precompute(self):
        """ Calculates the wages offered in each sector at each point in the state space. 
        The result is sol.w which we can then use directly, rather than using r and H. """
        #Precompute (for all ages) human capital #H has dimensions: s x a
        self.precompute_H()
        #Precompute wages (a x s x t) = (36, 3, 4)
        self.precompute_w()


    def precompute_H(self):
        """ Calculate all possible values of the human capital function. """
        self.sol.H = np.exp(self.par.beta0[:, np.newaxis]/self.par.scale.beta0 * self.par.ages)

    def precompute_w(self):
        """ Calculate all possible values of wages. w = r*H"""
        self.sol.w = np.transpose(self.sol.H)[:, :, np.newaxis] * self.par.r[np.newaxis, :, :]

    @staticmethod
    def EV_closedform(sigma, arr, axis=0):
        """ Calculate the expected value using the closed form logsum formula from the logit structure. 
            Numerically stabilized by subtracting a max() value."""
        if arr.ndim == 1:
            return sigma * (np.max(arr / sigma) + np.log(np.sum(np.exp(arr / sigma - np.max(arr / sigma)), axis)))
        elif arr.ndim == 2:
            return sigma * (np.max(arr / sigma, axis=axis) + \
                np.log(np.sum(np.exp(arr / sigma - np.max(arr / sigma, axis=axis)[:, np.newaxis]), axis)))

    #måske slå de her to sammen? så CCP og EV regnes samtidig, det er måske hurtigere? Nogle af de samme elementer indgår nemlig. (sum exp())
    @staticmethod
    def CCP_closedform(sigma, arr, axis=0):
        """ Calculate the conditional choice probabilities (CCPs) using the closed form formula from the logit structure"""
        if arr.ndim == 1:
            return np.exp(arr / sigma - np.max(arr / sigma)) / np.sum(np.exp(arr / sigma - np.max(arr / sigma)))
        elif arr.ndim == 2:
            return np.divide(np.exp(arr / sigma - np.max(arr / sigma, axis=axis)[:, np.newaxis]), 
                             np.sum(np.exp(arr / sigma - np.max(arr / sigma, axis=axis)[:, np.newaxis]), axis=axis)[:, np.newaxis])

    def solve_worker(self):
        """ Solve the worker's problem by backwards induction for a given set of skill prices (r). 
            First we solve the last period using static expectations. Then we perform backwards 
            iteration from the remaining periods until period t = 0 using rational expectations. """
        par = self.par
        sol = self.sol

        #Unpack solution objects
        w = sol.w #offered wages.
        M = par.M #utility switching cost
        c = sol.c #conditional choice probabilities
        v = sol.v #expected value functions

        #PART I (age iteration within the terminal period, static expectation)
        t = par.T - 1
        a = par.N_ages - 1

        #Start at the retirement age, where there are no continuation values
        #Loop over different lagged sectors (slag)

        for slag in np.arange(0, par.S):
            v[slag, a, t] = self.EV_closedform(par.sigma, w[a, :, t] - M[slag, :])
            c[slag, a, :, t] = self.CCP_closedform(par.sigma, w[a, :, t] - M[slag, :])

        #Then iterate from age 64 to 30.
        # a = 34
        for a in reversed(np.arange(par.N_ages - 1)):
            for slag in np.arange(0, par.S):
                V_alternatives = w[a, :, t] - M[slag, :] + par.rho * v[:, a + 1, t]
                #Nemt her stadig pga. simpel transformation function.
                #når valget i dag påvirker min state i morgen kommer det ind her.
                #Version 3: da mit valg i dag bare er det samme som min state i morgen kan jeg bare tilføje
                #v[:, ...] og stadig få den korrekte continuation value.
                v[slag, a, t] = self.EV_closedform(par.sigma, V_alternatives)
                c[slag, a, :, t] = self.CCP_closedform(par.sigma, V_alternatives)

        #PART II (time iteration from T - 2 to 0.)
        #Iterate from period T-2 (second to last period). Within each period. I don't have to iterate on age yet
        t = par.T - 2
        for t in reversed(np.arange(par.T - 1)):
            #Calculate the value function when age = 65, which has no continuation value
            a = par.N_ages - 1
            for slag in np.arange(0, par.S):
                v[slag, a, t] = self.EV_closedform(par.sigma, w[a, :, t] - M[slag, :])
                c[slag, a, :, t] = self.CCP_closedform(par.sigma, w[a, :, t] - M[slag, :])

            #Calculate the value functions for ages 64 - 30
            #H[:, 0:par.N_ages - 1] * r[:, t][:, np.newaxis] bruger mindre memory men kan ikke precomputes. Hvad er bedst?
            #Løn +M på tværs af valget. Samme dimension er i continuation value. løn i sek 0 skal plusses med cont når jeg kommer fra sek 0.
                V_alternatives =  w[0:par.N_ages - 1, :, t] - M[slag, :] + par.rho*v[:, 1:, t + 1].transpose()
                v[slag, 0:-1, t] = self.EV_closedform(par.sigma, V_alternatives, axis=1)
                c[slag, 0:-1, :, t] = self.CCP_closedform(par.sigma, V_alternatives, axis=1)



    def simulate(self):
        """ Simulate the model forward for a given set of skill prices. 
        We initialize the model in a point where all points of the state space are equally 
        populated. Then, we iterate forward in time and calculate the share of total employment accounted for by each point 
        in the state space. We do not have to draw gumbel shocks since we implicitly invoke a 'law of large numbers' and simply
        use the conditional choice probabilities as "transition" probabilities. """

        c = self.sol.c
        par = self.par

        #Measure how large a fraction of the people are at each point in the {state space x sector}.
        self.sim.density = np.zeros((par.S, par.N_ages, par.S, par.T))
        density = self.sim.density

        #This should be based on data later. Here we assume that all points in the state space are equally populated.
        #
        init_share = 1 / (par.S * par.N_ages * par.S)
        density[:, :, :, 0] = init_share

        #Do something different than equally populated for the sake of illustration and debugging
        density[:, :, :, 0] += np.array([- 1/2 * init_share, 0, 1/2 * init_share])[np.newaxis, np.newaxis, :]
        assert np.around(np.sum(density[:, :, :, 0]), 7) == 1

        #Specify the distribution across the state space for entering cohorts.
        #Eventually this could be based on the distribution of the cohort that entered the last year we have data (e.g. 2016)
        #Right now the only part of the state space that is relevant to specify is the lagged sector choice, since a = 30 by construction.
        #Just insert some non-uniform values
        enter_share = 1 / par.S 
        EnteringStateSpaceDist = np.zeros(shape=(par.S)) + enter_share + np.array([- 1/2 * enter_share, 0, 1/2 * enter_share])

        #loop over t starting from t = 0 going until the second to last period (which inserts values into the last). We only replace values in density in years 1, 2, ... 
        #Since the starting year is data.
        t = 1
        for t in np.arange(1, par.T):
            ## How many retired at the end the previous year?
            # The non-standard part is making sure people come into the model with age 30. 
            # np.sum(density[:, -1, :, t - 1]) sums over all parts of the state except age and time. Because we simply need ALL retiring people.
            retiring = np.sum(density[:, -1, :, t - 1])
            ## Entering cohort:
            density[:, 0, :, t] = retiring * EnteringStateSpaceDist[:, np.newaxis] * c[:, 0, :, t]

            ## What will people of ages 31-65 do in this period.
            # The transpose moves the sector dimension in front of the age dimension. This is intuitive, since the current sector
            # becomes next period's lagged sector. And it needs to line up with the policy function which has dimensions
            # (s_lag, age, s_curr, time) 
            # We sum over period t-2 sectors (previous period's lagged sector) since they do not matter
            # for current choices.
            density[:, 1:, :, t] = np.sum(density[:, 0:-1, :, t-1], axis=0).transpose()[:, :, np.newaxis] * c[:, 1:, :, t]

        assert all(np.around(np.sum(density[:, :, :, :], axis=(0, 1, 2)), 7) == 1)

    def solve_humancap_equilibrium(self, print_out=True):
        """ Calculate the skill price consistent with human capital demand equalizing human capital supply.
            This presumes that the worker's problem has already been solved once and simulated for some value of wages/skill prices.
            Currently, the equilibrium prices are found by successive approximations. """
        idx = pd.IndexSlice
        df = self.diag.tables.skillprice_converge

        #Proposed wages from the first iteration
        r1 = self.par.alpha1 * self.par.pY / (np.sum(self.sim.density, axis=(0,1)) * self.par.MASS)

        #Insert r0 and r1 for the current iteration
        iteration = 0
        df.loc[idx[:, :], idx[:, iteration]] = np.array([self.par.r.reshape(self.par.S * self.par.T), r1.reshape(self.par.S * self.par.T)]).T

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
                r1 = self.par.alpha1 * self.par.pY / (np.sum(self.sim.density, axis=(0,1)) * self.par.MASS)
                #Save the initial and proposed skill prices (for diagnostics plotting)
                df.loc[idx[:, :], idx[:, iteration]] = np.array([self.par.r.reshape(self.par.S * self.par.T), r1.reshape(self.par.S * self.par.T)]).T
                #Calculate deviation
                err_r = np.sum(np.abs(r1 - self.par.r))

    # def post_estimation(self):
    #     """This method calculates post-estimation quantities such as TFP (A), 
    #         physical capital (K), clean energy (E) and dirty energy (O)"""
    #     #This could change if we switch to the abatement form (JBE).
        


    def fig_skillprice_converge(self, save=False):
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
        d = np.sum(self.sim.density, axis=(0, 1))
        sector_names = self.par.sector_names

        fig = plt.figure(figsize=(5, 3.5), dpi=100)
        ax = fig.add_subplot(1, 1, 1)

        for i, name in enumerate(sector_names):
            ax.plot(d[i, :], label=name)

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

        for i, name in enumerate(sector_names):
            ax.plot(w[i, :], label=name)

        ax.legend(loc='upper center', frameon=True, bbox_to_anchor=(0.5, -0.15), ncol=self.par.S)
        ax.set_xlabel("Time")
        ax.set_title("Average wages")
        ax.set_ylabel("Average wage across ages")

        if save:
            self.savefig(fig, "average_wages")

    def uncond_switching_probs(self):
        """ Construct table of unconditional switching probabilities. """
        d = self.sim.density
        #Do not include period 0, since the transition calculated from that means period -1 to period 0, which lies outside the model
        probs = np.mean(np.sum(d[:, :, :, 1:], axis=1) / np.sum(d[:, :, :, 1:], axis=(1, 2))[:, np.newaxis, :], axis=2)
        self.results.tables.uncond_switching_probs = pd.DataFrame(probs, index=self.par.sector_names, columns=self.par.sector_names)

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

        n_params = 0
        for _, index in pars_to_estimate.items():
            if index is None:
                n_params += 1
            else:
                n_params += len(index)
        self.est.n_params = n_params

    def c_theta0(self):
        """ Constructs a single, 1-dimensional array containing all the parameters to be estimated.
            Requires self.pars_to_estimate to be defined (aliased pte here). I use the name theta0 as 
            the starting values of the parameters to be estimated. """
        pte = self.est.pars_to_estimate
        #Initialize theta0 with the correct length.
        theta0 = np.zeros(np.sum([len(pte[k]) if pte[k] is not None else 1 for k in pte]))
        n = 0 #counter for index in theta0 itself

        for par, index in pte.items():
            if index is None:
                theta0[n] = getattr(self.par, par)
                n += 1
            else:
                theta0[n : n + len(index)] = getattr(self.par, par)[index]
                n += len(index)
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

    def estimate(self, method="Nelder-Mead"):
        """ Estimate the parameters stored in self.est.pars_to_estimate. 
            Uses the starting values stored in self.est.theta0. 
        """
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
            self.solve_humancap_equilibrium()
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
        c = self.sol.c
        ll = self.loglik(d, c)
        if self.est.analytic_gradients:
            return (ll, self.score())
        else:
            return ll

    def loglik(self, d, c):
        """Calculate the (negative of the) log likelihood for a sample d with choice probabilities c"""
        return - np.sum(np.log(c[d["slag"], d["a"] - self.par.a_min, d["s"], d["t"]])) / len(d)

    def score(self):
        """Calculate the score of the (negative of the) log likelihood. """
        #unpack
        d = self.est.simulated_data
        dlnP = self.gradients(parameter="beta0")
        return - np.sum(dlnP[d["slag"], d["a"] - self.par.a_min, d["s"], d["t"], :], axis=0) / len(d)

    def score_from_theta(self, theta):
        """ Helper function for self.score. This makes score a function of theta so that it can be used 
            when running functions such as scipy.optimize.check_grad(f, g, theta). 
        """
        #unpack
        self.update_par(self.par, theta, self.est.pars_to_estimate)
        self.solve_and_simulate()
        return self.score()

    def c_simulated_data(self, n_individuals=100_000):
        """ Simulate data from the model currently stored. We simulate 'n_individuals' individuals. 
            The total number of observations (NxT) is random because we simulate individuals of all ages.
            This means that some people only enter the sample 1 year because that is the year they turn 30.
            Other individuals only enter the sample in 3 years because they start at the age of 63. """
        self.rng = default_rng(123456)
        data = pd.DataFrame(-1, 
                            index=pd.MultiIndex.from_product([range(0, n_individuals), range(self.par.T)], names=["i", "t"]), 
                            columns=["s", "a", "slag"]).reset_index()[["i", "slag", "a", "s", "t"]] #i for individual, t for year, s for sector (chosen), a for age, slag for lagged sector

        #Randomly generate starting ages that mean some are younger than 30 in the first period, so they enter model later
        data.loc[data.t == 0, "a"] = self.rng.integers(self.par.a_min - self.par.T + 1, self.par.a_max + 1, n_individuals)
        #Update the ages through time deterministically
        for t in range(1, self.par.T):
            data.loc[(data.t == t), "a"] = (data.loc[data.t == t - 1, "a"] + 1).to_numpy()
        #Delete rows where workers are not active in the labor market anymore
        #For example, this deletes period 0 observations for those that are 29 or younger in that year
        data = data.loc[((data.a <= self.par.a_max) & (data.a >= self.par.a_min))]

        #Randomly generate period -1 sector choices (period 0 sector lagged)
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
                            self.rng.choice(self.par.S, size=size, p=self.sol.c[slag, a - self.par.a_min, :, t])
        self.est.simulated_data = data.reset_index(drop=True)

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

    def gradients(self, parameter):
        """ Calculates analytic gradients of the log likelihood function. """
        
        #Unpack
        w = self.sol.w
        c = self.sol.c
        pte = self.est.pars_to_estimate

        #Containers for results
        dEV = np.zeros(self.sol.v.shape + (self.est.n_params,))
        dw = np.zeros(w.shape + (self.est.n_params,))
        dv = np.zeros(c.shape + (self.est.n_params,))

        if parameter == "beta0":

            #Index in theta (vector of parameters to be estimated)
            idx = np.array([0, 1]) # <--- brug denne korrekt når jeg får flere parametre! Lige nu er den lort

            #Derivative of wage wrt. beta0. Calculated for all state space points once-and-for-all
            dw[:, pte["beta0"], :, idx] = w[:, pte["beta0"], :].swapaxes(0, 1) * self.par.ages[np.newaxis, :, np.newaxis] / getattr(self.par.scale, parameter)

            a = self.par.a_max - self.par.a_min

            #Now we can fill in dv for the last age. There are no continuation values, so only the marginal wage effect matters
            dv[:, a, :, :, :] = dw[np.newaxis, a, :, :, :]

            #Next, the last age, dEV
            #Equation 41: Multiply P and dv for all s, and sum over s afterwards. We get positive values for all combinations of (slag x t),
            #because the expected value of the terminal period goes up no matter which of the wages we change through beta0. 
            dEV[:, a, :, :] = np.sum(c[:, a, :, :, np.newaxis] * dv[:, a, ...], axis=1)

            #ER SIGMA MED KORREKT RUNDT OMKRING? / SKRIV DET I OVERLEAF!

            # Derivative for younger people in the last period: use continuation values from age + 1 in the same period.
            # We have to iterate on age, because the continuation value of someone with age a uses the continuation value of someone
            # with age a + 1. 
            while a > 0:
                a -= 1
                dv[:, a, :, self.par.T - 1, :] = (dw[a, :, self.par.T - 1, :] + self.par.rho * dEV[:, a + 1, self.par.T - 1, :])[np.newaxis, :, :]
                #Here we match s in P with slag in dEV because the choice today becomes the lagged choice tomorrow
                dEV[:, a, self.par.T - 1, :] = np.sum(c[:, a, :, self.par.T - 1, np.newaxis] * dv[:, a, :, self.par.T - 1, :], axis=1) 

            #Now we can perform proper time iteration for the remaining ages.
            t = self.par.T - 1
            while t > 0:
                t -= 1
                #Choice specific value function derivatives
                dv[:, :-1, :, t, :] = (dw[:-1, :, t, :] + self.par.rho * dEV[:, 1:, t + 1, :].swapaxes(0, 1))[np.newaxis, ...]
                dEV[:, :-1, t, :] = np.sum(c[:, :-1, :, t, np.newaxis] * dv[:, :-1, :, t, :], axis=2)

            # This concludes the calculation of the expected value function derivatives (wrt. beta0)
            # Next, calculate the derivatives of the choice probabilities. 
            # This is made easier by the fact that we have already calculated dv.

            dlnP = np.zeros((c.shape) + (self.est.n_params,))
            #this is the dv of the choice minus a log sum, where after summing over k, we create the s dimension again
            dlnP = 1/self.par.sigma * (dv[:, :, :, :, :] - np.sum(c[:, :, :, :, np.newaxis] * dv[:, :, :, :, :], axis=2)[:, :, np.newaxis, :, :])

        return dlnP



# %% Run

gm = GM()
# gm.setup(beta0=np.array([-0.1, -0.4, 0.8]))
gm.setup(simple_statespace=False)
# gm.par.beta0 = gm.par.beta0 - 0.1

gm.c_human_capital_unit_prices()
gm.allocate()
gm.precompute()
gm.solve_worker()
gm.c_simulated_data()
gm.setup_estimation(method="BFGS", partial=True, analytic_gradients=True)
gm.c_pars_to_estimate(parameters=["beta0"], indexes=np.array([1, 2]))
gm.c_theta0()

#%%

method = "BFGS"
gm.est.options["gtol"] = 1e-8
#Make the initial values different from the true ones before testing estimation
# gm.est.theta0 = gm.est.theta0 * np.array([0.5, 1.5])

res = gm.estimate(method=method)
res

#%% Compare analytic and numeric gradients

gm.est.analytic_gradients = False
assert optimize.check_grad(gm.obj_func, gm.score_from_theta, np.array([0.003, 0.001])) < 1e-7, \
       "Differences between analytic and numeric gradients are too large"

#%% 



#%% Figur

#%%
gm.simulate()
gm.fig_avg_wages() #wages are unchanged from this shock
gm.fig_employment_shares(save=True)
# gm.uncond_switching_probs()
# print(gm.results.tables.uncond_switching_probs.to_latex())
#%%

gm = GM()
gm.setup()
gm.c_human_capital_unit_prices()
gm.allocate()
gm.precompute()
gm.solve()
gm.simulate()
# gm.fig_avg_wages()
gm.fig_employment_shares()


#%%


#%% Shocked system

gm = GM(name="Higher dirty wage")
gm.setup()
gm.c_human_capital_unit_prices()

#Add some more wage to the dirty sector
gm.par.r[-1, :] += 1

gm.allocate()
gm.precompute()
gm.solve()
gm.simulate()
gm.fig_avg_wages()
gm.fig_employment_shares()
