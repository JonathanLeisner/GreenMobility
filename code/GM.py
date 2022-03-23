
#%%

import numpy as np
import matplotlib.pyplot as plt
from types import SimpleNamespace
np.set_printoptions(precision=3)

# import warnings
# warnings.filterwarnings("error")


#%%

class GM:
    """ Class to solve the model presented in Leisner (2023)"""
    def __init__(self, name = None):
        self.sol = SimpleNamespace()
        self.par = SimpleNamespace()
        self.sim = SimpleNamespace()

        if name is None:
            self.name = "Standard"
        else:
            self.name = name

    def setup(self, **kwargs):
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
        
        #Not used
        # par.exper_min = 1e-6 # minimum point in grid for experience
        # par.exper_max = 20.0 # maximum point in grid for experience
        # par.gridpoints_exper = 200 #number of gridpoints for experience

        #Parameters that are calibrated
        par.rho = 0.96
        par.sigma = 1

        #Parameters to estimate later
        par.beta0 = np.array([0.002, 0.006, 0.012]) #1-dimensional over sectors

        #Replace default values by those given explicitly to the method
        for k, v in kwargs.items():
            setattr(self, k, v)

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

    def allocate(self):
        """ allocate empty containers for solution. 
        c contains the CCPs and v contains expected value functions"""
        sol = self.sol
        par = self.par

        sol.c = np.zeros((par.N_ages, par.S, par.T)) - 1
        sol.v = np.zeros((par.N_ages, par.T)) - 99

    def precompute(self):
        """ Calculates the wages offered in each sector at each point in the state space. 
        The result is sol.w which we can then use directly, rather than using r and H """
        #Precompute (for all ages) human capital #H has dimensions: s x a
        self.sol.H = np.exp(self.par.beta0[:, np.newaxis] * self.par.ages)

        #Precompute wages (a x s x t) = (36, 3, 4)
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

    def solve(self):
        """ Solve model by backwards induction. First we solve the last period using static
        expectations. Then we perform backwards iteration from the remaining periods until period t = 0
        using rational expectations. """
        par = self.par
        sol = self.sol

        #Unpack solution objects
        w = sol.w #offered wages.
        c = sol.c #conditional choice probabilities
        v = sol.v #expected value functions

        #PART I (age iteration within the terminal period, static expectation)
        t = par.T - 1
        a = par.N_ages - 1

        #Start at the retirement age, where there are no continuation values
        v[a, t] = self.EV_closedform(par.sigma, w[a, :, t])
        c[a, :, t] = self.CCP_closedform(par.sigma, w[a, :, t])

        #Then iterate from age 64 to 30.
        # a = 34
        for a in reversed(np.arange(par.N_ages - 1)):
            V_alternatives = w[a, :, t] + par.rho * v[a + 1, t] #Nemt her stadig pga. simpel transformation function.
                                                                #når valget i dag påvirker min state i morgen kommer det ind her.
            v[a, t] = self.EV_closedform(par.sigma, V_alternatives)
            c[a, :, t] = self.CCP_closedform(par.sigma, V_alternatives)

        #PART II (time iteration from T - 2 to 0.)
        #Iterate from period T-2 (second to last period). Within each period. I don't have to iterate on age yet
        t = par.T - 2
        for t in reversed(np.arange(par.T - 1)):
            #Calculate the value function when age = 65, which has no continuation value
            a = par.N_ages - 1
            v[a, t] = self.EV_closedform(par.sigma, w[a, :, t])
            c[a, :, t] = self.CCP_closedform(par.sigma, w[a, :, t])

            #Calculate the value functions for ages 64 - 30
            #H[:, 0:par.N_ages - 1] * r[:, t][:, np.newaxis] bruger mindre memory men kan ikke precomputes. Hvad er bedst?
            V_alternatives =  w[0:par.N_ages - 1, :, t] + par.rho*v[1:, t + 1][:, np.newaxis]    
            v[0:-1, t] = self.EV_closedform(par.sigma, V_alternatives, axis=1)
            c[0:-1, :, t] = self.CCP_closedform(par.sigma, V_alternatives, axis=1)

    def simulate(self):
        """ Simulate the model forward. We initialize the model in a point where all points of the state space are equally 
        populated. Then, we iterate forward in time and calculate the share of total employment accounted for by each point 
        in the state space. We do not have to draw gumbel shocks since we implicitly invoke a 'law of large numbers' and simply
        use the conditional choice probabilities as transition probabilities. """
        c = self.sol.c
        par = self.par

        #Measure how large a fraction of the people are at each point in the {state space x sector}.
        self.sol.density = np.zeros((par.N_ages, par.S, par.T))
        density = self.sol.density

        #This should be based on data later. Here we assume that all points in the state space are equally populated
        density[:, :, 0] = 1 / (par.N_ages * par.S)
        assert np.around(np.sum(density[:, :, 0], axis=(0, 1)), 7) == 1

        #loop over t starting from t = 0 going until the second to last period (which inserts values into the last). We only replace values in density in years 1, 2, ... 
        #Since the starting year is data.
        for t in np.arange(1, par.T):
            #The non-standard part is making sure people come into the model with age 30. 
            density[0, :, t] = np.sum(density[-1, :, t-1]) * c[0, :, t]
            #What will people of ages 31-65 do in this period.
            density[1:, :, t] = np.sum(density[0:-1, :, t-1], axis=1)[:, np.newaxis] * c[1:, :, t]

    def fig_employment_shares(self):
        """ Creates a figure of employment shares in each sector, over time"""
        #Unpack
        d = np.sum(self.sol.density, axis=0)
        sector_names = self.par.sector_names

        fig = plt.figure(figsize=(6, 4), dpi=100)
        ax = fig.add_subplot(1, 1, 1)

        for i, name in enumerate(sector_names):
            ax.plot(d[i, :], label=name)

        ax.legend(loc='upper left',frameon=True)
        ax.set_xlabel("Time")
        ax.set_ylabel("Employment share")
        ax.set_title("Title")

        fig.savefig("Employment_shares_" + self.name + ".pdf")

    def fig_avg_wages(self):
        """ Creates a figure of average wages across ages, over time"""
        #Unpack
        w = np.mean(self.sol.w, 0)
        sector_names = self.par.sector_names
        
        fig = plt.figure(figsize=(6, 4), dpi=100)
        ax = fig.add_subplot(1, 1, 1)

        for i, name in enumerate(sector_names):
            ax.plot(w[i, :], label=name)

        ax.legend(loc='upper left',frameon=True)
        ax.set_xlabel("Time")
        ax.set_ylabel("Average wage across ages")

        fig.savefig("Average_wages_" + self.name + ".pdf")

# %% Run

gm = GM()
gm.setup()
gm.create_human_capital_unit_prices()
gm.allocate()
gm.precompute()
gm.solve()
gm.simulate()


#%%

gm.fig_avg_wages()
gm.fig_employment_shares()

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
