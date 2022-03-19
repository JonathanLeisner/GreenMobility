
#%%

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from types import SimpleNamespace


#%%

class GM:
    def __init__(self):
        self.sol = SimpleNamespace()
        self.par = SimpleNamespace()
        self.sim = SimpleNamespace()

    def setup(self):
        par = self.par

        #Fundamentals    
        par.T = 4 #periods
        par.S = 3
        par.a_max = 65 #eldest cohort
        par.a_min = 30 #youngest cohort
        par.ages = np.arange(par.a_min, par.a_max + 1)
        par.N_ages = par.a_max + 1 - par.a_min
        par.exper_min = 1e-6 # minimum point in grid for experience
        par.exper_max = 20.0 # maximum point in grid for experience
        par.gridpoints_exper = 200 #number of gridpoints for experience

        #Parameters that are calibrated
        par.rho = 0.96

        #Parameters to estimate later
        par.beta0 = np.array([0.04, 0.05, 0.06]) #1-dimensional over sectors

    def create_grids(self):
        """ grids for experience """
        par = self.par
        par.grid_exper = np.linspace(par.exper_min, par.exper_max, par.gridpoints_exper)

    def create_human_capital_unit_prices(self):
        """ human capital unit prices. """
        par = self.par
        # np.linspace(t, par.S, par.T)
        # np.linspace(1, par.S, par.T)
        # SxT
        par.r = np.linspace(np.arange(1, par.S+1), np.arange(1, par.S+1)[::-1], par.T, axis=1)
        #first dimension: sector. Second dimension: time.
        # par.r = np.random.uniform(1, 10, size=(par.T, par.sectors))

    def allocate(self):
        """ allocate empty containers for solution"""
        sol = self.sol
        par = self.par

        sol.c = np.zeros((par.N_ages, par.T), dtype=int) - 99
        sol.v = np.zeros((par.N_ages, par.T))

    # def value_of_choice(d, t, a):
    #     par = self.par
    #     assert d in range(par.S)

    def precompute(self):
        #Precompute (for all ages) human capital #H has dimensions: s x a
        self.sol.H = np.exp(self.par.beta0[:, np.newaxis] * self.par.ages)

        #Precompute wages (a x s x t) = (36, 3, 4)
        self.sol.w = np.transpose(self.sol.H)[:, :, np.newaxis] * self.par.r[np.newaxis, :, :]

    def solve(self):
        """ Solve model by backwards induction """
        par = self.par
        sol = self.sol

        #Precompute
        self.precompute()

        #Unpack parameters
        r = par.r

        #Unpack solution objects
        H = sol.H
        w = sol.w
        c = sol.c
        v = sol.v

        #PART I (age iteration within the terminal period, static expectation)
        t = par.T - 1

        #Start at the retirement age, where there are no continuation values
        V_alternatives = H[:, -1] * r[:, t]
        c[-1, -1] = np.argmax(V_alternatives)
        v[-1, -1] = V_alternatives[c[-1, -1]] #np.max(V_alternatives). Hvad er hurtigst af indexering og genregning?

        #Then iterate from age 64 to 30.
        for a in reversed(np.arange(par.N_ages - 1)):
            V_alternatives = w[a, :, t] + par.rho*v[a + 1, t]
            c[a, -1] = np.argmax(V_alternatives)
            v[a, -1] = V_alternatives[c[a, -1]] #np.max(V_alternatives). Hvad er hurtigst?

        #PART II (time iteration from T - 2 to 0.)
        #Iterate from period T-2 (second to last period). Within each period. I don't have to iterate on age yet
        for t in reversed(np.arange(par.T - 1)):
            #Calculate the value function when age = 65, which has no continuation value
            a = par.N_ages - 1
            V_alternatives = w[a, :, t]
            c[a, t] = np.argmax(V_alternatives)
            v[a, t] = V_alternatives[c[a, t]]

            #Calculate the value functions for ages 64 - 30
            #H[:, 0:par.N_ages - 1] * r[:, t][:, np.newaxis] bruger mindre memory men kan ikke precomputes. Hvad er bedst?
            V_alternatives =  w[0:par.N_ages - 1, :, t] + par.rho*v[1:, t + 1][:, np.newaxis]
            c[:-1, t] = np.argmax(V_alternatives, axis=1) #argmax along s dimension
            v[0:-1, t] = V_alternatives[np.arange(V_alternatives.shape[0]), c[:-1, t]]

        assert (c == np.argmax(w, axis=1)).all()

        #experience er en skalar, som kan tage kontinuerte værdier. 
        # 
        # Så derfor skal vi bruge noget interpolation v. regression til at finde steder i statespacet
        # som ikke findes endnu. Det skal nok gå. Måske lav 2 types også, så vi kan få nogle cross-terms
        #?

    def simulate(self):
        """Primitive function for simulating individuals forward"""
        c = self.sol.c

        #initialize 2 persons at each age and simulate 35 periods forward. Do not generate new cohorts.
        #This means we create 36 cohorts and let the number of cohorts decrease by 1 each year.
        N_sim = self.par.N_ages * 2
        ss = np.zeros((N_sim, self.par.T), dtype=int) - 1 #only 1-dimensional apart from N x T since the i-specific part of state space only consists of age
        sol = np.zeros((N_sim, self.par.T), dtype=int) - 1

        #Initialize the state space
        ss[:, 0] = self.par.ages.repeat(2)

        #Loop forward through years
        for t in np.arange(self.par.T):
            active_i = ss[:, t] <= 65 #only do calculations on people who are active in the labor market
            i = ss[active_i, t] - self.par.a_min
            sol[active_i, t] = c[i, t]
            if t < self.par.T - 1:
                #Transition function (add 1 to age)
                ss[:, t + 1] = ss[:, t] + 1

        self.sim.sol = sol
        self.sim.ss = ss

    def fig_employment(self):
        fig = plt.figure()
        ax = fig.add_subplot()
        df = pd.DataFrame(self.sim.sol, dtype=float)
        df[df < 0] = np.nan

        translate = {
            0:"Unemployment",
            1:"Clean",
            2:"Dirty"
        }

        df.replace(translate, inplace=True)
        data = df.apply(pd.value_counts).replace(np.nan, 0).transpose()
        for col in data.columns:
            ax.plot(data.index, data[col], label=col)

        ax.legend()
        fig.savefig("Employment_v1.pdf")



# %% Testing

gm = GM()
gm.setup()
gm.create_human_capital_unit_prices()
gm.allocate()
gm.solve()
gm.simulate()
gm.fig_employment()



#%% Testing Gumbel stuff.

kappa = 2
loc = -0.57721*kappa
vals = np.random.gumbel(loc=loc, scale=kappa, size=100000)

print(f'Mean {vals.mean()}.\nVariance: {vals.var()}')

np.pi**2/6 * kappa**2