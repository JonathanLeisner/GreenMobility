"""  Defines the class that runs the Green Mobility (GM) project code."""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#%%

class output:
    """ Class to construct output in the form of figures and tables for the model presented in Leisner (2023). 
    """
    def __init__(self, par, sol, sim, est, name=None, version=None):
        if name is None:
            name = "std"
        self.name = name
        if version is None:
            version = "v3"
        self.version = version
        self.resultspath = "../results/"

        self.par = par
        self.sol = sol
        self.sim = sim
        self.est = est

    def savefig(self, fig, figname):
        """Save the figure object 'fig' as a .pdf with name 'figname'. """
        fig.savefig(self.resultspath + figname +  "_" + self.name + "_" + self.version + ".pdf", bbox_inches='tight')
        plt.close(fig)

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
        else:
            plt.show()

    def fig_avg_wages(self, save=False):
        """ Creates a figure of average wages across ages, over time"""
        #Unpack
        w = np.mean(self.sol.w, axis=0)
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
        else:
            plt.show()
            
            
    def avg_yearly_transition_rates(self, data=True):
        """ Construct table of unconditional switching probabilities. This calculates only net-flows. Should be recomputed for gross-flows.
            e.g. such that with just two sectors, if workers are 50/50 dsitributed  and they all switch sectors, the probabilities do not become 0. """ 
        if data:
            d = self.est.simulated_data
            df = d.groupby(["slag", "t"])["s"].value_counts(normalize=True).mean(level=(0, 2)).unstack()
        else:
            d = self.sim.density
            df = pd.DataFrame(np.mean(np.sum(d, axis=1) / np.sum(d, axis=(1, 3))[:, :, np.newaxis], axis=1), 
                              index=self.par.sector_names, 
                              columns=self.par.sector_names)
        return df.to_latex()


    def age_profile_switching(self, save):
        P = self.sol.P

        # 1) Model switching
        d = self.sim.density
        #Fraction of each age-year combination that are in different lagged sectors. So if we sum within age and time, over slag, we get 1.
        D = np.sum(d, axis=-1) / np.sum(d, axis=(0, -1))[np.newaxis, ...] #fraction i hvert state space punkt. 

        #Make a dummy that specifies what elements constitute switches:
        a = np.repeat(np.arange(self.par.S), self.par.S)
        b = np.tile(np.arange(self.par.S), self.par.S)
        slag, s, = a[~(a == b)], b[~(a == b)]
        dummy = np.zeros(P.shape)
        dummy[slag, :, :, s] = 1

        #Switching probabilities in the simulation
        switching_prob = np.sum(P * dummy, axis=-1)
        age_profile_model = np.mean(np.sum(D * switching_prob, axis=0), axis=1)

        # 2) Data switching
        d = self.est.simulated_data
        d["switch"] = False
        d.loc[d["slag"] != d["s"], "switch"] = True
        age_profile_data = d.groupby(["a", "t"])["switch"].value_counts(normalize=True).mean(level=(0, 2)).loc[(slice(None), True)]

        #Construct figure
        fig = plt.figure(figsize=(5, 3.5), dpi=100)
        ax = fig.add_subplot(1, 1, 1)

        ax.plot(self.par.ages, age_profile_model, label="Model prediction")
        ax.plot(self.par.ages, age_profile_data, label="Data")

        ax.legend(loc='upper center', frameon=True, bbox_to_anchor=(0.5, -0.15))
        ax.set_xlabel("Age")
        ax.set_ylabel("Switch rate")
        ax.set_title("Age profile of switching probability")
        
        if save:
            self.savefig(fig, "age_profile_switching")
        else:
            plt.show()

    def time_profile_switching(self, save):
        P = self.sol.P

        # 1) Model switching
        d = self.sim.density
        D = np.sum(d, axis=-1)

        #Make a dummy that specifies what elements constitute switches:
        a = np.repeat(np.arange(self.par.S), self.par.S)
        b = np.tile(np.arange(self.par.S), self.par.S)
        slag, s, = a[~(a == b)], b[~(a == b)]
        dummy = np.zeros(P.shape)
        dummy[slag, :, :, s] = 1

        switching_prob = np.sum(P * dummy, axis=-1)
        time_profile_model = np.sum(D * switching_prob, axis=(0, 1))

        # 2) Data switching
        d = self.est.simulated_data
        d["switch"] = False
        d.loc[d["slag"] != d["s"], "switch"] = True
        time_profile_data = d.groupby("t")["switch"].value_counts(normalize=True).loc[(slice(None), True)]

        #Construct figure
        fig = plt.figure(figsize=(5, 3.5), dpi=100)
        ax = fig.add_subplot(1, 1, 1)

        ax.plot(np.arange(0, self.par.T), time_profile_model, label="Model prediction")
        ax.plot(np.arange(0, self.par.T), time_profile_data, label="Data")

        ax.legend(loc='upper center', frameon=True, bbox_to_anchor=(0.5, -0.15))
        ax.set_xlabel("Year")
        ax.set_ylabel("Switch rate")
        ax.set_title("Time profile of switching probabilities")

        if save:
            self.savefig(fig, "time_profile_switching")
        else:
            plt.show()



# %%
