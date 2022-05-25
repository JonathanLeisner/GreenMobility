
#%%

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from types import SimpleNamespace
np.set_printoptions(precision=3)

import seaborn as sns
sns.set_theme()
# import warnings
# warnings.filterwarnings("error")


#%%

class GM:
    """ Class to solve the model presented in Leisner (2023)"""
    def __init__(self, name = None):
        self.sol = SimpleNamespace()
        self.par = SimpleNamespace()
        self.sim = SimpleNamespace()
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

    def setup(self):
        """ Initiates parameter values at their default values"""
        
        par = self.par

        #Fundamentals    
        par.T = 8 #periods
        par.sector_names = ["Unemployment", "Clean", "Dirty"]
        par.S = len(par.sector_names)
        par.a_max = 65 #eldest cohort
        par.a_min = 30 #youngest cohort
        par.ages = np.arange(par.a_min, par.a_max + 1)
        par.N_ages = par.a_max + 1 - par.a_min
        
        #Meta-settings
        sol = self.sol
        sol.step_fraction = 0.05
        sol.maxiter = 500
        sol.tolerance_r = 0.0001

        #Not used
        # par.exper_min = 1e-6 # minimum point in grid for experience
        # par.exper_max = 20.0 # maximum point in grid for experience
        # par.gridpoints_exper = 200 #number of gridpoints for experience

        #Fixed quantities from data
        par.MASS = self.load_MASS()
        par.pY = self.load_nominal_output()

        #Parameters that are calibrated
        par.rho = 0.96
        par.sigma = 1

        #Should be loaded later
        par.alpha1 = np.repeat(np.array([0.2, 0.3, 0.4])[:, np.newaxis], par.T, axis=1)
        par.alpha2 = np.repeat(np.array([0.2, 0.2, 0.2])[:, np.newaxis], par.T, axis=1)
        par.alpha3 = 1 - par.alpha1 - par.alpha2

        #Parameters to estimate later
        #Human capital function parameters
        par.beta0 = np.array([0.002, 0.006, 0.012]) #1-dimensional over sectors. 

        #Switching costs (direct utility costs)
        if self.name == "high_switchingcost":
            par.xi_in = 20 * np.array([0.04, 0.05, 0.06])
            par.xi_out = 20 * np.array([0.03, 0.06, 0.09]).transpose()    
        else:
            par.xi_in = np.array([0.04, 0.05, 0.06])
            par.xi_out = np.array([0.03, 0.06, 0.09]).transpose()
        self.compute_switching_cost_matrix()

        self.diag.tables.skillprice_converge = pd.DataFrame(np.nan, 
                                                             index=pd.MultiIndex.from_product([gm.par.sector_names, range(gm.par.T)], 
                                                                                              names=["Sector", "Year"]), 
                                                             columns=pd.MultiIndex.from_product([["r0", "r1"], range(0, gm.sol.maxiter)], 
                                                             names=["Pre/post", "Iteration"]))

        # #Replace default values by those given explicitly to the method
        # for k, v in kwargs.items():
        #     setattr(self, k, v)

    def compute_switching_cost_matrix(self):
        """ Takes the xi-parameters and constructs the matrix of switching costs."""
        par = self.par

        #Switching cost (SC) has dimensions (s_{t-1} x s_{t}) i.e. (s_out x s_in)
        par.SC = sum(np.meshgrid(par.xi_in, par.xi_out)) #add cost of going OUT of a sector with the cost of going IN to a sector.
        np.fill_diagonal(par.SC, 0) #Not switching sector is costless
        par.SC = np.exp(par.SC)

    # def create_grids(self):
    #     """ grids for experience """
    #     par = self.par
    #     par.grid_exper = np.linspace(par.exper_min, par.exper_max, par.gridpoints_exper)

    def create_human_capital_unit_prices(self):
        """ Human capital unit prices, the symbol r in the paper """
        par = self.par
        # SxT
        par.r = np.linspace(np.arange(1, par.S+1), np.arange(1, par.S+1)[::-1], par.T, axis=1)
        #first dimension: sector. Second dimension: time.
        # par.r = np.random.uniform(1, 10, size=(par.T, par.sectors))

    def load_MASS(self):
        #for testing, this makes the population double halfway through the sample
        return np.concatenate([np.ones(int(round(self.par.T/2, 0))), np.ones(self.par.T - int(round(self.par.T/2, 0))) + 1])

    def load_nominal_output(self):
        return np.repeat((np.ones(self.par.S))[:, np.newaxis], self.par.T, axis=1)

    def allocate(self):
        """ allocate empty containers for solution. 
        c contains the CCPs and v contains expected value functions"""
        sol = self.sol
        par = self.par

        sol.c = np.zeros((par.S, par.N_ages, par.S, par.T)) - 1 #4 dimensions: previous sector, age, chosen sector, time
        sol.v = np.zeros((par.S, par.N_ages, par.T)) - 99 #3 dimensions: previous sector, current age, time period

    def precompute(self):
        """ Calculates the wages offered in each sector at each point in the state space. 
        The result is sol.w which we can then use directly, rather than using r and H """
        #Precompute (for all ages) human capital #H has dimensions: s x a
        self.precompute_H()
        #Precompute wages (a x s x t) = (36, 3, 4)
        self.precompute_w()

    def precompute_H(self):
        self.sol.H = np.exp(self.par.beta0[:, np.newaxis] * self.par.ages)

    def precompute_w(self):
        self.sol.w = np.transpose(self.sol.H)[:, :, np.newaxis] * self.par.r[np.newaxis, :, :]

    @staticmethod
    def EV_closedform(sigma, arr, axis=0):
        """ Calculate the expected value using the closed form logsum formula from the logit structure """
        return sigma * np.log(np.sum(np.exp(arr / sigma), axis))

    #måske slå de her to sammen? så CCP og EV regnes samtidig, det er måske hurtigere? Nogle af de samme elementer indgår nemlig. (sum exp())
    @staticmethod
    def CCP_closedform(sigma, arr, axis=0):
        """ Calculate the conditional choice probabilities (CCPs) using the closed form formula from the logit structure"""
        if arr.ndim <= 1:
            return np.exp(arr / sigma) / np.sum(np.exp(arr / sigma))
        else:
            return np.divide(np.exp(arr / sigma), np.sum(np.exp(arr / sigma), axis)[:, np.newaxis])

    def solve_worker(self):
        """ Solve model by backwards induction for a given set of skill prices. First we solve the last period using static
        expectations. Then we perform backwards iteration from the remaining periods until period t = 0
        using rational expectations. """
        par = self.par
        sol = self.sol

        #Unpack solution objects
        w = sol.w #offered wages.
        SC = par.SC #utility switching cost
        c = sol.c #conditional choice probabilities
        v = sol.v #expected value functions

        #PART I (age iteration within the terminal period, static expectation)
        t = par.T - 1
        a = par.N_ages - 1

        #Start at the retirement age, where there are no continuation values
        #Loop over different lagged sectors (slag)

        for slag in np.arange(0, par.S):
            v[slag, a, t] = self.EV_closedform(par.sigma, w[a, :, t] - SC[slag, :])
            c[slag, a, :, t] = self.CCP_closedform(par.sigma, w[a, :, t] - SC[slag, :])

        #Then iterate from age 64 to 30.
        # a = 34
        for a in reversed(np.arange(par.N_ages - 1)):
            for slag in np.arange(0, par.S):
                V_alternatives = w[a, :, t] - SC[slag, :] + par.rho * v[:, a + 1, t]
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
                v[slag, a, t] = self.EV_closedform(par.sigma, w[a, :, t] - SC[slag, :])
                c[slag, a, :, t] = self.CCP_closedform(par.sigma, w[a, :, t] - SC[slag, :])

            #Calculate the value functions for ages 64 - 30
            #H[:, 0:par.N_ages - 1] * r[:, t][:, np.newaxis] bruger mindre memory men kan ikke precomputes. Hvad er bedst?
            #Løn +SC på tværs af valget. Samme dimension er i continuation value. løn i sek 0 skal plusses med cont når jeg kommer fra sek 0.
                V_alternatives =  w[0:par.N_ages - 1, :, t] - SC[slag, :] + par.rho*v[:, 1:, t + 1].transpose()
                v[slag, 0:-1, t] = self.EV_closedform(par.sigma, V_alternatives, axis=1)
                c[slag, 0:-1, :, t] = self.CCP_closedform(par.sigma, V_alternatives, axis=1)



    def simulate(self):
        """ Simulate the model forward for a given set of skill prices. 
        We initialize the model in a point where all points of the state space are equally 
        populated. Then, we iterate forward in time and calculate the share of total employment accounted for by each point 
        in the state space. We do not have to draw gumbel shocks since we implicitly invoke a 'law of large numbers' and simply
        use the conditional choice probabilities as transition probabilities. """

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

    def solve_humancap_equilibrium(self):
        #calculate the skill price consistent with human capital demand equalizing human capital supply.
        #This presumes that the worker's problem has already been solved once for some value of wages/skill prices.
        idx = pd.IndexSlice
        df = self.diag.tables.skillprice_converge

        #Proposed wages from the first iteration
        r1 = self.par.alpha1 * self.par.pY / (np.sum(self.sim.density, axis=(0,1)) * self.par.MASS)

        #Insert r0 and r1 for the current iteration
        iteration = 0
        df.loc[idx[:, :], idx[:, iteration]] = np.array([self.par.r.reshape(self.par.S * self.par.T), r1.reshape(self.par.S * self.par.T)]).T

        err_r = np.sum(np.abs(r1 - self.par.r))
        for iteration in range(1, self.sol.maxiter):
            if err_r < 0.00001:
                print(f"Number of iterations when converged was: {iteration}")
                break
            else:
                if iteration == self.sol.maxiter - 1:
                    raise Exception(f"Human capital equilibrium could not be found after {self.sol.maxiter} iterations.")
                #print(f"Current error of skill prices is: {err_r:.7f}")
                #Make another iteration
                #Update the skill prices (and wages) then resolve and simulate 
                self.par.r = r1 * self.sol.step_fraction + self.par.r * (1 - self.sol.step_fraction)
                self.precompute_w()
                self.solve_worker()
                self.simulate()
                #Proposed skill prices (S x T)
                r1 = self.par.alpha1 * self.par.pY / (np.sum(self.sim.density, axis=(0,1)) * self.par.MASS)
                #Save the initial and proposed skill prices (for diagnostics plotting)
                df.loc[idx[:, :], idx[:, iteration]] = np.array([self.par.r.reshape(self.par.S * self.par.T), r1.reshape(self.par.S * self.par.T)]).T
                #Calculate deviation
                err_r = np.sum(np.abs(r1 - self.par.r))

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
        self.results.tables.uncond_switching_probs = pd.DataFrame(probs, index=gm.par.sector_names, columns=gm.par.sector_names)

#%% Testing



# %% Run

gm = GM()
gm.setup()
gm.create_human_capital_unit_prices()
gm.allocate()
gm.precompute()
gm.solve_worker()
gm.simulate()
gm.solve_humancap_equilibrium()
gm.fig_skillprice_converge(save=True)
#%%




#%% Iterate to find equilibrium skill prices



#%% Figur

#Remove the nans (iterations that were not reached)



#%%
gm.simulate()
gm.fig_avg_wages() #wages are unchanged from this shock
gm.fig_employment_shares(save=True)
# gm.uncond_switching_probs()
# print(gm.results.tables.uncond_switching_probs.to_latex())
#%%

gm = GM()
gm.setup()
gm.create_human_capital_unit_prices()
gm.allocate()
gm.precompute()
gm.solve()
gm.simulate()
# gm.fig_avg_wages()
gm.fig_employment_shares()



#%%


#%%


#%% Shocked system

gm = GM(name="Higher dirty wage")
gm.setup()
gm.create_human_capital_unit_prices()

#Add some more wage to the dirty sector
gm.par.r[-1, :] += 1

gm.allocate()
gm.precompute()
gm.solve()
gm.simulate()
gm.fig_avg_wages()
gm.fig_employment_shares()

#%% Testing Gumbel stuff.

kappa = 2
loc = -0.57721*kappa
vals = np.random.gumbel(loc=loc, scale=kappa, size=100000)

print(f'Mean {vals.mean()}.\nVariance: {vals.var()}')

#%% Old simulation files where we simulate individuals


#initialize 2 persons at each age and simulate 35 periods forward. Do not generate new cohorts.
#This means we create 36 cohorts and let the number of active cohorts decrease by 1 each year.
# N_sim = self.par.N_ages * 2

# ss = np.zeros((N_sim, self.par.T), dtype=int) - 1 #only 1-dimensional apart from N x T since the i-specific part of state space only consists of age
# sol = np.zeros((N_sim, self.par.T), dtype=int) - 1
# #Initialize the state space
# ss[:, 0] = self.par.ages.repeat(2)

# #Loop forward through years
# for t in np.arange(self.par.T):
#     active_i = ss[:, t] <= 65 #only do calculations on people who are active in the labor market
#     i = ss[active_i, t] - self.par.a_min
#     sol[active_i, t] = c[i, t]
#     if t < self.par.T - 1:
#         #Transition function (add 1 to age)
#         ss[:, t + 1] = ss[:, t] + 1

#     self.sim.sol = sol
#     self.sim.ss = ss
