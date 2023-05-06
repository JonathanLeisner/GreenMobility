"""  Defines the class that runs the Green Mobility (GM) project code."""

from types import SimpleNamespace
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import util
import python_to_latex as pytex

#%%

plt.rcParams.update({
    "legend.framealpha":1,
    "text.usetex":True,
    "text.latex.preview":True, #Makes vertical alignment between marker and label in legend better
    "font.family":"serif",
    "font.serif":  ['Computer Modern Roman',
                    'DejaVu Sans', 
                    'Bitstream Vera Sans',
                    'Lucida Grande',
                    'Verdana',
                    'Geneva',
                    'Lucid',
                    'Arial',
                    'Helvetica',
                    'Avant Garde',
                    'sans-serif']
})

class output_basic:
    """ Class to construct output in the form of figures and tables for the model presented in Leisner (2023).
    """
    # pylint: disable=no-member
    # pylint: disable=multiple-statements
    # pylint: disable=attribute-defined-outside-init
    def c_basic_output_settings(self):
        self.settings = SimpleNamespace()
        settings = self.settings
        settings.resultspath = "../output/"
        settings.translate_sectors = {i:self.par.sector_names[i] for i in range(self.par.S)}
        short_names = {"Non-employment":"Non-emp.", "Food":"Food", "Energy":"Energy", "Transport":"Transport", "Manufacturing":"Manuf.", "Services":"Services"}
        settings.translate_sectors_short = {i:short_names[self.par.sector_names[i]] for i in range(self.par.S)}
        
        settings.note_strings = {}
        settings.note_strings["relative"] = "Each value is the relative difference between the policy" +\
            " scenario and the baseline scenario of \\cref{appsec:counterfactual baseline scenario} where all taxes rates are constant."
        settings.note_strings["difference"] = "Each value is the difference between policy" +\
            " scenario and the baseline scenario of \\cref{appsec:counterfactual baseline scenario} where all taxes rates are constant."

        if self.par.T < 9:
            settings.yeargap = 0
            settings.labelrotation = 25
        elif self.par.T < 20:
            settings.yeargap = 1
            settings.labelrotation = 25
        elif self.par.T < 30:
            settings.yeargap = 2
            settings.labelrotation = 0
        elif self.par.T < 60:
            settings.labelrotation = 0
            settings.yeargap = 5
        else:
            settings.labelrotation = 0
            settings.yeargap = 20
        settings.colors = "#1f77b4 #ff770e #2ca02c #d62728 #9467bd #8c564b #e377c2 #7f7f7f #bcbd22 #17becf".split()
        #Define font sizes
        settings.fontsizes = SimpleNamespace()
        self.c_fontsizes()

    def c_fontsizes(self):
        """ 
        Defines the font sizes used in figures.
        """ 
        self.settings.fontsizes.ylabel = 12
        self.settings.fontsizes.xlabel = 12

    def figures_to_pdf(self, filename="TEST"):
        """
        Write all the figures currently stored in self.results.figures to a PDF using latex rendering. 
        """
        pytex.write_figs_to_latex(self.results.figures, filename)
        self.results.figures = []

    def tables_to_pdf(self, filename="TEST"):
        """
        Write all the tables currently stored in self.results.tables to a PDF using latex rendering.
        """
        pytex.write_tbls_to_latex(self.results.tables, filename)
        self.results.tables = []

class output_counterfactual_simulation(output_basic):
    # pylint: disable=no-member
    # pylint: disable=multiple-statements
    def c_output_settings(self):
        self.c_basic_output_settings()
        self.settings.init_year = 2016

    def f_setup_sector_by_time(self, array):
        fig = plt.figure(figsize=(5, 3.5), dpi=100)
        ax = fig.add_subplot(1, 1, 1)

        if array.shape[-1] == 6:
            names = self.par.sector_names
            sector_add = 0
        elif array.shape[-1] == 5:
            names = self.par.sector_names[1:]
            sector_add = 1
        for s, name in enumerate(names):
            ax.plot(np.arange(0, self.par.T) + self.settings.init_year, array[:, s], 
                    label=name, color=self.settings.colors[s + sector_add])
        return (fig, ax)

    def f_settings_sector_by_time(self, fig, ax, ylabel, note, label, caption, compare):
        if not compare:
            ax.set_ylim([0.0, ax.get_ylim()[1]])
        ax.legend(loc='upper center', frameon=True, bbox_to_anchor=(0.5, -0.15), ncol=int(np.ceil(self.par.S / 2)))
        ax.set_xlabel("Year", fontsize=self.settings.fontsizes.ylabel)
        ax.xaxis.set_major_locator(ticker.IndexLocator(base=1+self.settings.yeargap, offset=0))
        ax.set_ylabel(ylabel, fontsize=self.settings.fontsizes.ylabel)
        ax.set_axisbelow(True) #Makes grey linear appear underneath other graphical elements
        ax.yaxis.grid(color="lightgrey")
        plt.show()
        self.results.figures += [pytex.TexFigure(fig, caption=caption, label=label, note=note)]

    def f_employment_shares(self, latex=None, caption="Employment shares", compare=False):
        """ Creates a figure of employment shares in each sector, over time"""
        emp = np.sum(self.sim.density, axis=tuple(i for i in range(self.sim.density.ndim - 2)))
        if compare: 
            emp = emp - np.sum(self.comparison_scenario.sim.density, axis=tuple(i for i in range(self.sim.density.ndim - 2)))
        (fig, ax) = self.f_setup_sector_by_time(emp)
        label = util.suffix_to_str("cf_employment_shares", latex)
        if compare:
            note = "The figure shows employment share changes for each sector over time."
            ylabel = "Change in employment share"
        else:
            note = "The figure shows employment shares for each sector over time."
            ylabel = "Employment share"
        if compare: note += " " + self.settings.note_strings["difference"]
        note += " " + self.tau_scenario_description
        self.f_settings_sector_by_time(fig, ax, ylabel=ylabel, note=note, label=label, caption=caption, compare=compare)

    def f_r(self, latex=None, caption="Human capital prices", compare=False):
        """ Creates a figure of human capital prices shares in each sector, over time"""
        r = self.par.r
        if compare: r = (r - self.comparison_scenario.par.r) / self.comparison_scenario.par.r
        (fig, ax) = self.f_setup_sector_by_time(r)
        label = util.suffix_to_str("cf_r", latex)
        if compare:
            note = "The figure shows relative changes in human capital prices for each sector over time."
            ylabel = "Relative human capital price change"
        else:
            note = "The figure shows human capital prices for each sector over time."
            ylabel = "Human capital price"
        if compare: note += " " + self.settings.note_strings["relative"]
        note += " " + self.tau_scenario_description 
        self.f_settings_sector_by_time(fig, ax, ylabel=ylabel, note=note, label=label, caption=caption, compare=compare)

    def f_Y1z(self, latex=None, caption="Output produced", compare=False):
        """ Creates a figure of production in each sector, over time"""
        Y1z = self.sol.Y1z 
        if compare: 
            Y1z = (Y1z - self.comparison_scenario.sol.Y1z) / self.comparison_scenario.sol.Y1z
        else:
            Y1z = Y1z * self.par.scale.MASS / 1000
        (fig, ax) = self.f_setup_sector_by_time(Y1z)
        label = util.suffix_to_str("cf_Y1z", latex)
        if compare:
            note = "The figure shows relative changes in production for each sector over time."
            ylabel = "Relative change in production"
        else:
            note = "The figure shows production for each sector over time."
            ylabel = "Production (billion DKK)"
        if compare: note += " " + self.settings.note_strings["relative"]
        note += " " + self.tau_scenario_description 
        self.f_settings_sector_by_time(fig, ax, ylabel=ylabel, note=note, label=label, caption=caption, compare=compare)

    def f_p(self, latex=None, caption="Output prices", compare=False):
        """ Creates a figure of output capital prices in each sector, over time"""
        p = self.sol.p 
        if compare: p = (p - self.comparison_scenario.sol.p) / self.comparison_scenario.sol.p
        (fig, ax) = self.f_setup_sector_by_time(p)
        label = util.suffix_to_str("cf_p", latex)
        if compare:
            note = "The figure shows relative changes in output prices for each sector over time."
            ylabel = "Relative change in output price"
        else:
            note = "The figure shows output prices for each sector over time."
            ylabel = "Output price"
        if compare: note += " " + self.settings.note_strings["relative"]
        note += " " + self.tau_scenario_description
        self.f_settings_sector_by_time(fig, ax, ylabel=ylabel, note=note, label=label, caption=caption, compare=compare)

    def f_avg_offered_wages(self, latex=None, caption="Average offered wages over time", compare=False):
        """ Creates a figure of average offered wages across state space, over time"""
        wage = np.mean(self.sol.w, axis=0) 
        if compare: wage = (wage - np.mean(self.comparison_scenario.sol.w, axis=0)) / np.mean(self.comparison_scenario.sol.w, axis=0)
        (fig, ax) = self.f_setup_sector_by_time(wage)
        label = util.suffix_to_str("cf_avg_offered_wages", latex)
        if compare:
            note = "The figure shows the relative change in the average offered wage for each sector over time."
            ylabel = "Relative change in average offered wage"
        else:
            note = "The figure shows the average offered wage for each sector over time."
            ylabel = "Average offered wage (million DKK)"
        if compare: note += " " + self.settings.note_strings["relative"]
        note += " " + self.tau_scenario_description 
        self.f_settings_sector_by_time(fig, ax, ylabel=ylabel, note=note, label=label, caption=caption, compare=compare)

    def f_avg_realized_wages(self, latex=None, caption="Average realized wages over time", compare=False):
        wage = np.average(self.sol.w, axis=tuple(np.arange(self.sol.w.ndim - 2)), weights=self.sim.density)
        if compare:
            wage0 = np.average(self.comparison_scenario.sol.w, axis=tuple(np.arange(self.sol.w.ndim - 2)), weights=self.comparison_scenario.sim.density)
            wage = (wage - wage0) / wage0
        wage = wage[:, 1:] #Remove unemployment "wage"
        (fig, ax) = self.f_setup_sector_by_time(wage)
        label = util.suffix_to_str("cf_avg_realized_wages", latex)
        if compare:
            note = "The figure shows the relative change in the average wage for each sector over time."
            ylabel = "Relative change in average wage"
        else:
            note = "The figure shows the average wage for each sector over time."
            ylabel = "Average wage (million DKK)"
        if compare: note += " " + self.settings.note_strings["relative"]
        note += " Each average wage is weighted by the number of individuals in the state space point. " + self.tau_scenario_description
        self.f_settings_sector_by_time(fig, ax, ylabel=ylabel, note=note, label=label, caption=caption, compare=compare)

    def f_Hsup(self, latex=None, caption="Aggregate human capital supply", compare=False):
        """ Creates a figure of human capital supply in each sector, over time"""
        Hsup = self.sol.Hsup
        if compare: 
            Hsup = (Hsup - self.comparison_scenario.sol.Hsup) / self.comparison_scenario.sol.Hsup
        else:
            Hsup = Hsup * self.par.scale.MASS / 1_000_000
        (fig, ax) = self.f_setup_sector_by_time(Hsup)
        label = util.suffix_to_str("cf_Hsup", latex)
        if compare:
            note = "The figure shows the relative change in human capital supply for each sector over time."
            ylabel = "Relative change in human capital supply"
        else:
            note = "The figure shows human capital supply for each sector over time."
            ylabel = "Human capital supply (millions)"
        if compare: note += " " + self.settings.note_strings["relative"]
        note += " " + self.tau_scenario_description
        self.f_settings_sector_by_time(fig, ax, ylabel=ylabel, note=note, label=label, caption=caption, compare=compare)

    def f_z(self, latex=None, caption="Abatement share", compare=False):
        """ Creates a figure of the abatement share in each sector, over time"""
        z = self.sol.z 
        if compare: z = (z - self.comparison_scenario.sol.z)
        (fig, ax) = self.f_setup_sector_by_time(z)
        label = util.suffix_to_str("cf_z", latex)
        if compare:
            note = "The figure shows the change in the abatement share for each sector over time."
            ylabel = "Change in abatement share"
        else:
            note = "The figure shows the abatement share for each sector over time."
            ylabel = "Abatement share"
        if compare: note += " " + self.settings.note_strings["difference"]
        note += " " + self.tau_scenario_description
        self.f_settings_sector_by_time(fig, ax, ylabel=ylabel, note=note, label=label, caption=caption, compare=compare)

    def f_Z(self, latex=None, caption="Emissions", compare=False):
        Z = self.sol.g * self.sol.Y1z / (1 - self.sol.z) 
        if compare:
            Z0 = self.comparison_scenario.sol.g * self.comparison_scenario.sol.Y1z / (1 - self.comparison_scenario.sol.z)
            Z = (Z - Z0) / Z0
        else:
            Z = Z * self.par.scale.MASS
        (fig, ax) = self.f_setup_sector_by_time(Z)
        label = util.suffix_to_str("cf_emissions", latex)
        if compare:
            note = "The figure shows the relative change in emissions for each sector over time."
            ylabel = "Relative change in emissions"
        else:
            note = "The figure shows emissions for each sector over time."
            ylabel = "Emissions (1000 tons CO2e)" 
        if compare: note += " " + self.settings.note_strings["relative"]
        note += " " + self.tau_scenario_description
        self.f_settings_sector_by_time(fig, ax, ylabel=ylabel, note=note, label=label, caption=caption, compare=compare)

    def f_tau(self, latex=None, caption="Emission tax rates", compare=False):
        tau = self.par.tau
        if compare: tau = (tau - self.comparison_scenario.par.tau) / self.comparison_scenario.par.tau
        (fig, ax) = self.f_setup_sector_by_time(tau)
        label = util.suffix_to_str("cf_tau", latex)
        if compare:
            note = "The figure shows the relative change in the emission tax rate for each sector over time."
            ylabel = "Relative change in emission tax rate"
        else:
            note = "The figure shows emission tax rates for each sector over time."
            ylabel = r"Emission tax rate $\left(\frac{\mathrm{million\ DKK}}{\mathrm{1000\ tons\ CO2e}}\right)$"
        if compare: 
            note += " " + self.settings.note_strings["relative"]
        note += " " + self.tau_scenario_description
        self.f_settings_sector_by_time(fig, ax, ylabel=ylabel, note=note, label=label, caption=caption, compare=compare)

    def f_emsint(self, latex=None, caption="Emission intensities", compare=False):
        emsint = self.sol.g / (1 - self.sol.z)
        if compare:
            emsint0 = self.comparison_scenario.sol.g / (1 - self.comparison_scenario.sol.z)
            emsint = (emsint - emsint0) / emsint0
        (fig, ax) = self.f_setup_sector_by_time(emsint)
        label = util.suffix_to_str("cf_emsint", latex)
        if compare:
            note = "The figure shows the relative change in the emission intensity for each sector over time."
            ylabel = "Relative change in emission intensity"
        else:
            note = "The figure shows the emission intensity for each sector over time."
            ylabel = r"Emission intensity $\left(\frac{\mathrm{1000\ tons\ CO2e}}{\mathrm{million\ DKK}}\right)$"
        if compare: note += " " + self.settings.note_strings["relative"]
        note += " " + self.tau_scenario_description
        self.f_settings_sector_by_time(fig, ax, ylabel=ylabel, note=note, label=label, caption=caption, compare=compare)

    def f_K(self, latex=None, caption="Physical capital", compare=False):
        K = self.sol.K
        if compare: 
            K = (K - self.comparison_scenario.sol.K) / self.comparison_scenario.sol.K
        else:
            K = K * self.par.scale.MASS / 1_000_000
        (fig, ax) = self.f_setup_sector_by_time(K)
        label = util.suffix_to_str("cf_K", latex)
        if compare:
            note = "The figure shows the relative change in physical capital for each sector over time."
            ylabel = "Relative change in physical capital"
        else:
            note = "The figure shows physical capital for each sector over time."
            ylabel = "Physical capital (million DKK)"
        if compare: note += " " + self.settings.note_strings["relative"]
        note += " Physical capital is assumed to adjust such that the return on capital is fixed " +\
            "at their (sector-specific) 2016-values. " + self.tau_scenario_description
        self.f_settings_sector_by_time(fig, ax, ylabel=ylabel, note=note, label=label, caption=caption, compare=compare)

    def t_EV(self, latex=None, caption="Expected value functions"):
        assert hasattr(self, "comparison_scenario")
        a_grp_size = 6
        deviation = (self.sol.EV[..., 0] - self.comparison_scenario.sol.EV[..., 0]) / np.abs(self.comparison_scenario.sol.EV[..., 0])
        deviation = deviation.reshape(self.sol.EV.shape[:-2] + (-1, a_grp_size), order="C")
        density_weights = np.sum(self.sim.density, axis=-1)[..., 0].reshape(self.sol.EV.shape[:-2] + (-1, a_grp_size), order="C")
        out = np.average(deviation, axis=(0, -1), weights=density_weights)
        col_names = [f"{self.par.a_min + i * a_grp_size} - {self.par.a_min + a_grp_size - 1 + i * a_grp_size}" for i in np.arange(a_grp_size)]
        df = (pd.DataFrame(out, index=self.par.sector_names, columns=col_names).reset_index()
              .rename(columns={"index":"Initial sector $\\downarrow\\;$ / $\\;$ Initial age $\\rightarrow$"}))
        note = f"The table shows a weighted average of the relative change in the expected value function," +\
            " i.e. expected lifetime utility from the labor market (see \\cref{eq:expected lifetime utility}), between" +\
            " the policy scenario and the baseline scenario of \\cref{appsec:counterfactual baseline scenario}, measured in year $t=$ 2016." +\
            " The average is taken over all tenures and ages within the particular group and each relative change is weighted by the number" +\
            " of people in each state space point. " + self.tau_scenario_description
        label = util.suffix_to_str("EV_deviations", latex)
        self.results.tables += [pytex.TexTable(df, caption=caption, decimals=4, label=label, note=note)]
        return df

class output_estimate(output_basic):
    # pylint: disable=no-member
    def c_output_settings(self):
        """
        Populates the attribute 'settings' to control various options for plotting tables and figures.
        """
        self.c_basic_output_settings()
        settings = self.settings
        
        #Used for separating estimated parameters
        settings.non_sector_parameters = ["sigma", "kappa0", "kappa1", "phi0", "phi1"]
        settings.sector_parameters = ["beta0", "beta1", "xi_in", "xi_out", "beta2"]
        settings.translate_parameters = {"kappa0":"$\\kappa_0$", "kappa1":"$\\kappa_1$", "phi0":"$\\phi_0$", 
                                         "phi1":"$\\phi_1$", "sigma":"$\\sigma$", "beta0":"$\\beta_s^0$",
                                         "xi_in":"$\\xi_s^{in}$", "xi_out":"$\\xi_s^{out}$",
                                         "beta2":"$\\beta_s^2$", "beta1":"$\\beta_s^1$"}
        
        settings.descriptions = {"sigma":"Scale parameter of preference shock",
                                 "kappa0":"Switching cost shifter, age",
                                 "kappa1":"Switching cost shifter, age squared",
                                 "phi0":"Non-employment value, age",
                                 "phi1":"Non-employment value, age squared",
                                 "beta0":"Human capital, age",
                                 "beta1":"Human capital, age squared",
                                 "xi_in":"Switching cost, into sectors",
                                 "xi_out":"Switching cost, out of sectors",
                                 "beta2":"Human capital, tenure"
                                 }

        settings.index_corrections = {"xi_in":1, "beta0":1, "beta1":1, "beta2":1} #e.g. add 1 to xi_in's index to get sector (since it has no value for sector 0)
        settings.init_year = 2002

    def f_employment_shares(self, latex=None, caption="Employment shares"):
        """ Creates a figure of employment shares in each sector, over time"""
        #Unpack
        emp = np.sum(self.sim.density, axis=tuple(i for i in range(self.sim.density.ndim - 2)))
        sector_names = self.par.sector_names

        emp_data = self.est.data.groupby(["t"])["s"].value_counts(normalize=True)

        fig = plt.figure(figsize=(5, 3.5), dpi=100)
        ax = fig.add_subplot(1, 1, 1)

        for s, name in enumerate(sector_names):
            ax.plot(np.arange(0, self.par.T) + self.settings.init_year, emp[:, s], 
                    label=name, color=self.settings.colors[s])
            ax.plot(np.arange(0, self.par.T) + self.settings.init_year, emp_data.loc[:, s], 
                    label=name + " (data)", color=self.settings.colors[s], linestyle="dashed")

        ax.legend(loc='upper center', frameon=True, bbox_to_anchor=(0.5, -0.15), ncol=int(np.ceil(self.par.S / 2)))
        ax.set_xlabel("Year", fontsize=self.settings.fontsizes.ylabel)
        ax.xaxis.set_major_locator(ticker.IndexLocator(base=1+self.settings.yeargap, offset=0))
        ax.set_ylabel("Employment share", fontsize=self.settings.fontsizes.ylabel)
        ax.set_axisbelow(True) #Makes grey linear appear underneath other graphical elements
        ax.yaxis.grid(color="lightgrey")
        plt.show()
        note = "The figure shows employment shares for each sector over time." + \
               " The employment shares from the data are compared to those simulated from the model during estimation." +\
               " The employment shares from the model simulation are calculated as" +\
               r" $\sum_{\omega \in \widetilde{\Omega}} D_t(\omega)P\left(d_{\omega s t } = 1 | \;\omega, \bm{r}\right)$ for each $s$."

        label = util.suffix_to_str("employment_shares", latex)
        self.results.figures += [pytex.TexFigure(fig, caption=caption, label=label, note=note)]

    def f_age_employment_shares(self, latex=None, caption="Employment across ages and sectors"):
        """ Creates a figure of employment shares in each sector, across ages"""
        sector_names = self.par.sector_names
        emp = self.est.data.groupby(["a", "s", "t"]).size().unstack().mean(axis=1).unstack().to_numpy()

        fig = plt.figure(figsize=(5, 3.5), dpi=100)
        ax = fig.add_subplot(1, 1, 1)

        for s, name in enumerate(sector_names):
            ax.plot(self.par.ages, emp[:, s], 
                    label=name, color=self.settings.colors[s])

        ax.legend(loc='upper center', frameon=True, bbox_to_anchor=(0.5, -0.15), ncol=int(np.ceil(self.par.S / 2)))
        ax.set_xlabel("Age", fontsize=self.settings.fontsizes.ylabel)
        ax.xaxis.set_major_locator(ticker.IndexLocator(base=2, offset=0))
        ax.set_ylabel("Average employment", fontsize=self.settings.fontsizes.ylabel)
        ax.set_axisbelow(True) #Makes grey linear appear underneath other graphical elements
        ax.yaxis.grid(color="lightgrey")
        plt.show()
        note = "The figure shows the number of people choosing each sector at each age in the estimation sample." + \
               " The values are averages across all years in the estimation sample."

        label = util.suffix_to_str("age_employment_shares", latex)
        self.results.figures += [pytex.TexFigure(fig, caption=caption, label=label, note=note)]

    def f_avg_realized_wages(self, latex=None, caption="Average realized wages"):
        """ Creates a figure of average offered wages across state space, over time"""
        #Unpack
        w = np.average(self.sol.w, axis=tuple(np.arange(self.sol.w.ndim - 2)), weights=self.sim.density)
        sector_names = self.par.sector_names

        fig = plt.figure(figsize=(5, 3.5), dpi=100)
        ax = fig.add_subplot(1, 1, 1)

        for s, name in enumerate(sector_names[1:]):
            ax.plot(np.arange(0, self.par.T) + self.settings.init_year, w[:, s + 1], label=name, color=self.settings.colors[s + 1])

        ax.xaxis.set_major_locator(ticker.IndexLocator(base=1+self.settings.yeargap, offset=0))
        ax.set_ylim([0.0, ax.get_ylim()[1]])

        ax.legend(loc='upper center', frameon=True, bbox_to_anchor=(0.5, -0.20), ncol=int(np.ceil(self.par.S / 2)))
        ax.set_xlabel("Year", fontsize=self.settings.fontsizes.xlabel)
        ax.set_ylabel("Average wage", fontsize=self.settings.fontsizes.ylabel)

        ax.tick_params(axis="x", which="major", labelsize=10, labelrotation=self.settings.labelrotation)

        ax.set_axisbelow(True) #Makes grey linear appear underneath other graphical elements
        ax.yaxis.grid(color="lightgrey")

        plt.show()
        note = "The figure shows average realized wages for the different sectors." +\
               " The values shown are averages over the state space weighted by the number of individuals in that point."
        label = util.suffix_to_str("avg_realized_wages", latex)
        self.results.figures += [pytex.TexFigure(fig, caption=caption, label=label, note=note)]


            
    def t_avg_yearly_transition_rates(self, data=True, caption=None, latex=None):
        """ Construct a table of average yearly transition rates. 
            The argument 'data' determines whether it is calculated from the data or from the model simulation.
        """ 
        if data:
            translate = pd.Series(self.par.sector_names, index=np.arange(self.par.S))
            d = self.est.data
            df = d.groupby(["slag", "t"])["s"].value_counts(normalize=True).mean(level=(0, 2)).unstack().reset_index()
            df["slag"].replace(translate, inplace=True)
            df.rename(columns=dict(translate), inplace=True)
            df.rename(columns={"slag":"From $\\downarrow\\;$ / $\\;$To $\\rightarrow$"}, inplace=True)
        else:
            d = self.sim.density
            df = pd.DataFrame(np.mean(np.sum(d, axis=tuple(np.arange(d.ndim - 4)) + (-3,))/np.sum(d, axis=tuple(np.arange(d.ndim - 4)) + (-3, -1))[..., np.newaxis], 
                                      axis=-2), 
                              index=self.par.sector_names, 
                              columns=self.par.sector_names).reset_index()
            df.rename(columns={"index":"From $\\downarrow\\;$ / $\\;$ To $\\rightarrow$"}, inplace=True)
        data_or_model = {True:'data', False:'model'}.get(data)
        note = f"The table shows the average yearly sector switching rates in the {data_or_model}"
        if not data:
            note += " during estimation."
        else:
            note += "."
        if caption is None:
            caption = "Average yearly sector switching rates"
            if data:
                caption += " (data)"
            else:
                caption += " (model)"
        label = util.suffix_to_str("avg_yearly_transition_rates_" + data_or_model, latex)
        self.results.tables += [pytex.TexTable(df, caption=caption, decimals=4, label=label, note=note)]
        return df

    def f_age_profile_switching(self, latex=None, data_only=False, caption="Age profile switching"):
        """ Construct a figure of the age profile of sector switching."""
        P = self.sol.P

        if not data_only:
            # 1) Model switching
            d = self.sim.density
            #Fraction of each age-year combination that are in different lagged sectors. So if we sum within age and time, over slag, we get 1.
            # D = np.sum(d, axis=-1) / np.sum(d, axis=(0, -1))[np.newaxis, ...] <- delete 
            #Measure each point in the state space as a fraction of all those with the same (age, time)
            D = np.sum(d, axis=-1) / np.sum(d, axis=tuple(np.arange(d.ndim - 3)) + (-1,))      

            #Make a dummy that specifies what elements constitute switches:
            a = np.repeat(np.arange(self.par.S), self.par.S)
            b = np.tile(np.arange(self.par.S), self.par.S)
            slag, s, = a[~(a == b)], b[~(a == b)]
            dummy = np.zeros(P.shape)
            dummy[..., slag, :, :, s] = 1

            #Switching probabilities in the simulation
            switching_prob = np.sum(P * dummy, axis=-1)
            #Density-weighted average of switching probabilities
            age_profile_model = np.mean(np.sum(D * switching_prob, axis=tuple(np.arange(D.ndim - 2))), axis=-1)

        # 2) Data switching
        d = self.est.data.copy()
        d["switch"] = False
        d.loc[d["slag"] != d["s"], "switch"] = True
        age_profile_data = d.groupby(["a", "t"])["switch"].value_counts(normalize=True).mean(level=(0, 2)).loc[:, True]

        #Construct figure
        fig = plt.figure(figsize=(5, 3.5), dpi=100)
        ax = fig.add_subplot(1, 1, 1)

        ax.plot(self.par.ages, age_profile_data, label="Data", color=self.settings.colors[0])
        if not data_only:
            ax.plot(self.par.ages, age_profile_model, label="Model", color=self.settings.colors[1])

        ax.legend(loc='upper center', frameon=True, bbox_to_anchor=(0.5, -0.2), ncol=2)
        ax.set_xlabel("Age", fontsize=self.settings.fontsizes.xlabel)
        ax.set_ylabel("Switching rate", fontsize=self.settings.fontsizes.ylabel)
        ax.tick_params(axis="x", which="major", labelsize=10, labelrotation=self.settings.labelrotation)
        ax.set_axisbelow(True) #Makes grey linear appear underneath other graphical elements
        ax.yaxis.grid(color="lightgrey")
        plt.show()
        note = "The figure shows the average switching rate of all individuals of a given age, in the data and in the model during estimation."

        label_prefix = "age_profile_switching" + {False:"", True:"_data_only"}.get(data_only)
        label = util.suffix_to_str(label_prefix, latex)
        self.results.figures += [pytex.TexFigure(fig, caption=caption, label=label, note=note)]
    

    @staticmethod
    def scale_as_string(parameter, par):
        """
        Defines a string equal to the parameter's scale for use in tables/figures.
        """
        if hasattr(par.scale, parameter) and (getattr(par.scale, parameter) != 1):
            return f" $(\\times {getattr(par.scale, parameter)})$"
        else:
            return ""

    def t_compare_init_to_estimated(self, latex=None, caption="Comparison of initial values to estimated ones"):
        df = self.compare_init_to_estimated().reset_index()
        note = "The table shows maximum likelihood estimates of the model parameters compared to the initially used values." + \
               " 'SE' refers to standard errors."
        label = util.suffix_to_str("init_values_to_estimates", latex)
        self.results.tables += [pytex.TexTable(df, caption=caption, note=note, label=label, decimals=5)]

    def t_estimates_nosectordim(self, latex=None, caption="Estimates of parameters without a sector dimension"):
        """
        Tabulate the estimates without a sector dimension.
        """
        df = self.compare_init_to_estimated()
        df = df.reset_index()
        df = df[df["Parameter"].isin(self.settings.non_sector_parameters)].copy()
        df = df[["Parameter", "Estimated", "SE"]].copy()
        
        df["Description"] = df["Parameter"].replace(pd.Series(self.settings.descriptions))
        df = pd.melt(df, id_vars=["Parameter", "Description"], value_vars=["Estimated", "SE"]).sort_values("Parameter")

        df["scale"] = df["Parameter"].apply(lambda x: self.scale_as_string(x, self.par))
        df["Parameter"] = df["Parameter"].replace(pd.Series(self.settings.translate_parameters))
        df["Parameter"] = df["Parameter"] + df["scale"]

        df.loc[df["variable"] == "SE", ["Parameter", "Description"]] = ""
        # df.loc[df["variable"] == "SE", "Description"] = ""
        df.loc[df["variable"] == "SE", "value"] = df.loc[df["variable"] == "SE", "value"].apply(lambda x: f"({round(x, 5)})")
        df.rename(columns={"value":"Estimate"}, inplace=True)
        df.drop(columns=["variable", "scale"], inplace=True)

        note = "The table shows maximum likelihood estimates of the model parameters that have no sector dimension." + \
               " Standard errors are in parenthesis."
        label = util.suffix_to_str("estimates_nosectordim", latex)
        self.results.tables += [pytex.TexTable(df, caption=caption, note=note, label=label, decimals=3)]
        return df

    def t_estimates_sectordim(self, latex=None, caption="Estimates of parameters with a sector dimension"):
        """
        Constructs a table with estimation results for the parameters that have a sector dimension.
        """
        df = self.compare_init_to_estimated().reset_index()
        df = df[df["Parameter"].isin(self.settings.sector_parameters)].copy()
        for p in self.settings.index_corrections:
            df.loc[df["Parameter"] == p, "Index"] += self.settings.index_corrections[p]
        df.rename(columns={"Index":"Sector"}, inplace=True)
        df = df[["Parameter", "Sector", "Estimated", "SE"]].copy()
        df["Description"] = df["Parameter"].replace(pd.Series(self.settings.descriptions))
        df = pd.melt(df, id_vars=["Parameter", "Sector", "Description"], value_vars=["Estimated", "SE"])
        sectors = list(df["Sector"].unique())
        df = df.pivot_table(index=["Parameter", "variable", "Description"], columns="Sector", values="value").reset_index()
        df.columns.name = None
        df["Parameter"] = df["Parameter"].replace(pd.Series(self.settings.translate_parameters)) + df["Parameter"].apply(lambda x: self.scale_as_string(x, self.par))
        df.loc[df["variable"] == "SE", ["Parameter", "Description"]] = ""
        df.loc[df["variable"] == "SE", sectors] = \
            df.loc[df["variable"] == "SE", sectors].applymap(lambda x: f"({round(x, 2)})" if not np.isnan(x) else np.nan)
        df.drop(columns=["variable"], inplace=True)
        df.rename(columns=self.settings.translate_sectors_short, inplace=True)
        note = "The table shows maximum likelihood estimates of the model parameters that have a sector dimension." + \
                " Standard errors are in parenthesis."
        label = util.suffix_to_str("estimates_sectordim", latex)
        self.results.tables += [pytex.TexTable(df, caption=caption, note=note, label=label, decimals=2)]
        return df

    def t_median_switching_costs(self, latex=None, caption="Direct utility costs of switching sectors"):
        data = self.est.data.copy()
        #Condition on working in the previous period
        data = data[data.slag > 0]
        #Wage when staying in slag
        w_staying = self.sol.w[data.ten - 1, data.slag, data.a - self.par.a_min, data.t, data.slag]
        #Switching costs into all other sectors than slag (i.e. other than staying)
        M_switching = self.par.M.copy()
        M_switching[np.arange(self.par.S), :, np.arange(self.par.S)] = np.nan
        M_switching = M_switching[data.slag, data.a - self.par.a_min, 1:]

        df = pd.DataFrame(index=self.par.sector_names[1:], columns=["median"])
        df.loc[:, "median"] = np.nanmedian(M_switching/w_staying[:, np.newaxis], axis=0)
        df = df.reset_index().rename(columns={"median":"Median cost of entry", "index":"Sector of entry"})
        note = "The table reports median switching costs of entering a sector conditional on being employed in the previous year." + \
            " The costs are calculated as $\\frac{M(s, \\Omega_{it})}{w_{is_{it-1}t}}$" + \
            " with combinations where $s=s_{t-1}$ and hence $M=0$ are left out." + \
            " The wage in the denominator is the wage worker $i$ in year $t$ in the data would receive if she stayed in the sector she worked in in year $t-1$." + \
            " The switching cost is calculated once for each observation in the data, and the table reports the median for each potential sector of entry."
        label = util.suffix_to_str("median_switching_costs", latex)
        self.results.tables += [pytex.TexTable(df, caption=caption, note=note, label=label, decimals=3)]
        return df
# %%
