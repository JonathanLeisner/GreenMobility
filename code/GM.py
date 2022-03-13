
#%%
import pandas as pd
import numpy as np
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

        sol.c = np.zeros((par.N_age, par.T))
        sol.v = np.zeros((par.N_age, par.T))

    def value_of_choice(d, t, a):
        par = self.par
        assert d in range(par.S)


    def solve(self):
        """ Solve model by backwards induction """
        par = self.par
        sol = self.sol

        #Solve the terminal period, where workers have static expectations.
        #Start at the age of retirement


        # sol.c[-1, :] = par.r[-1, :] * 2
        #sol.v[-1, :] = par.grid_m**(1-par.rho)/(1-par.rho) 

        

        #experience er en skalar, som kan tage kontinuerte værdier. 
        # 
        # Så derfor skal vi bruge noget interpolation v. regression til at finde steder i statespacet
        # som ikke findes endnu. Det skal nok gå. Måske lav 2 types også, så vi kan få nogle cross-terms
        #?

        # def value_of_choice(d, t, a):




# %% Testing

gm = GM()
gm.setup()
gm.create_human_capital_unit_prices()

#unpack
par = gm.par
r = par.r
beta0 = par.beta0
ages = par.ages

c = np.zeros((len(ages), par.T), dtype=int) - 99
v = np.zeros((len(ages), par.T))
#PART I
t = par.T - 1

#precompute (for all ages) human capital #H has dimensions: s x a
H = np.exp(beta0[:, np.newaxis] * ages)

#start at last age, 65, where there are no continuation values
V_alternatives = H[:, -1] * r[:, t]
c[-1, -1] = np.argmax(V_alternatives)
v[-1, -1] = V_alternatives[c[-1, -1]] #np.max(V_alternatives). Hvad er hurtigst af indexering og genregning?

#Calculate wages (a x s x t) = (36, 3, 4)
w = np.transpose(H)[:, :, np.newaxis] * r[np.newaxis, :, :]

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
