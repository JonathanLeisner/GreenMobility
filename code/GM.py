
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
        par.a_min = 25 #youngest cohort
        par.exper_min = 1e-6 # minimum point in grid for experience
        par.exper_max = 20.0 # maximum point in grid for experience
        par.gridpoints_exper = 200 #number of gridpoints for experience

        #Parameters to estimate later
        par.beta0 = np.array([0.4, 0.5, 0.6])


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

        sol.c = np.zeros((par.T, par.gridpoints_exper))
        sol.v = np.zeros((par.T, par.gridpoints_exper))



    def value_of_choice(d, t, a):
        par = self.par
        assert d in range(par.S)


    def solve(self):
        """ Solve model by backwards induction """
        par = self.par
        sol = self.sol

        sol.c[-1, :] = par.r[-1, :] * 2
        sol.v[-1, :] = par.grid_m**(1-par.rho)/(1-par.rho) 



        #experience er en skalar, som kan tage kontinuerte værdier. 
        # 
        # Så derfor skal vi bruge noget interpolation v. regression til at finde steder i statespacet
        # som ikke findes endnu. Det skal nok gå. Måske lav 2 types også, så vi kan få nogle cross-terms
        #?

    # def value_of_choice(d, t, a):




# %% Testing

gm = GM()
gm.setup()
# gm.create_grids()

gm.create_human_capital_unit_prices()
par = gm.par
r = par.r
beta0 = par.beta0


# def value_of_choice(d, t, model):
t = par.T - 1
r[:, t]



ages = np.arange(par.a_min, par.a_max + 1)

beta0 = par.beta0


H = np.exp(beta0[:, np.newaxis] * ages)

par.beta0 = np.array([0.02, 0.03, 0.04])

H.shape

r[:, par.T - 1].shape

Us = r[:, par.T - 1, np.newaxis] * H


Us.argmax(0)

Us.shape


gm.allocate()
par = gm.par

#%%

 np.outer(par.r[-1, :], par.grid_exper).max(axis=0)



x = np.random.uniform(1, 10, size=(40, 2))

