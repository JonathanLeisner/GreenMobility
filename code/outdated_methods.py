 # def solve_humancap_equilibrium(self, print_out=True):
    #     """ Calculate the skill price consistent with human capital demand equalizing human capital supply.
    #         This presumes that the worker's problem has already been solved once and simulated for some value of wages/skill prices.
    #         Currently, the equilibrium prices are found by successive approximations. """
    #     idx = pd.IndexSlice
    #     df = self.diag.tables.skillprice_converge

    #     #Proposed wages from the first iteration
    #     r1 = self.par.alpha1 * self.par.pY / (np.sum(self.sim.density, axis=(0,1)) * self.par.MASS[:, np.newaxis])

    #     #Insert r0 and r1 for the current iteration
    #     iteration = 0
    #     df.loc[idx[:, :], idx[:, iteration]] = np.array([self.par.r.reshape(self.par.T * self.par.S), r1.reshape(self.par.T * self.par.S)]).T

    #     err_r = np.sum(np.abs(r1 - self.par.r))
    #     for iteration in range(1, self.sol.maxiter):
    #         if err_r < self.sol.tolerance_r:
    #             if print_out:
    #                 print(f"Number of iterations when converged was: {iteration}")
    #                 break
    #         else:
    #             if iteration == self.sol.maxiter - 1:
    #                 raise Exception(f"Human capital equilibrium could not be found after {self.sol.maxiter} iterations.")
    #             #print(f"Current error of skill prices is: {err_r:.7f}")
    #             #Make another iteration
    #             #Update the skill prices (and wages) then resolve and simulate 
    #             self.par.r = r1 * self.sol.step_fraction + self.par.r * (1 - self.sol.step_fraction)
    #             self.precompute_w() #H need not be recomputed within this inner loop, since it is the unaltered given some vector of parameters.
    #             self.solve_worker()
    #             self.simulate()
    #             #Proposed skill prices (S x T)
    #             r1 = self.par.alpha1 * self.par.pY / (np.sum(self.sim.density, axis=(0,1)) * self.par.MASS[:, np.newaxis])
    #             #Save the initial and proposed skill prices (for diagnostics plotting)
    #             df.loc[idx[:, :], idx[:, iteration]] = np.array([self.par.r.reshape(self.par.S * self.par.T), r1.reshape(self.par.S * self.par.T)]).T
    #             #Calculate deviation
    #             err_r = np.sum(np.abs(r1 - self.par.r))

    # def fig_skillprice_converge(self, save=False):
    #     """Figure showing the time series of unit skill prices (r) at various iterations from the human capital equilibrium function. """
    #     idx = pd.IndexSlice
    #     df = self.diag.tables.skillprice_converge
    #     #Remove nans (iterations that were not reached before convergence)
    #     df = df.loc[:, ~df.isna().all(axis=0)]

    #     fig = plt.figure(figsize=(5, 3.5), dpi=100)
    #     ax = fig.add_subplot(1, 1, 1)
    #     r = "r0"
    #     ax.plot(df.loc[idx["Unemployment", :], idx[r, 0]].unstack(level="Sector"), marker="o", color="green")
    #     ax.plot(df.loc[idx["Unemployment", :], idx[r, 5::20]].unstack(level="Sector"), color="orange", alpha=0.7)
    #     ax.plot(df.loc[idx["Unemployment", :], idx[r, df.columns.get_level_values(level="Iteration")[-1]]].unstack(level="Sector"), marker="x", color="red")
    #     #todo: skriv grøn: start etc. i en legend
    #     #todo: latex-skrift for r og t
    #     # ax.legend(loc='upper center', frameon=True, bbox_to_anchor=(0.5, -0.15), ncol=self.par.S)
    #     ax.set_xlabel("Time (t)")
    #     ax.set_title("Convergence of skill prices")
    #     ax.set_ylabel("Skill price (r)")

    #     if save:
    #         self.savefig(fig, "skillprice_converge")