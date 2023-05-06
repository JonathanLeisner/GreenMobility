"""  Defines the class that runs the Green Mobility (GM) project code."""
#%%
import os
os.chdir("I:/Workdata/706651/users/Jonathan/Projects/GreenMobility/code")
import numpy as np
np.set_printoptions(precision=3)
import estimate

#%% Estimate the model on actual data
#Setup
data_name = "standard"
gm = estimate.GM_estimate(name="estimated_v4", endo_D=True)
gm.partial_model = False
# gm.setup(from_data=True, sigma=np.array([5.0]), 
#          beta0=np.array([-65.0, -23.0, -280, -28, -24]),
#          xi_in=np.array([26.0, 25.0, 27, 22, 21]), xi_out=np.array([3.4, 1, 0, -0.5, 0.0, 1.5]),
#          kappa0=np.array([3.4]), kappa1=np.array([-2.8]), phi0=np.array([-22.0]), phi1=np.array([6.0]))
gm.setup(from_data=True, sigma=np.array([5.0]))
gm.l_data(data_name)
gm.precompute()
print(f"The number of observations used for the estimation is {len(gm.est.data)}.")
gm.allocate_worker()
gm.solve_worker()
gm.simulate()
gm.successive_approximations(90)
gm.sim.endo_D = False
gm.find_humcap_equilibrium()
# gm.setup_estimation(parameters=["xi_in", "xi_out", "kappa0", "phi0", "kappa1", "phi1", "beta0", "sigma"])
# gm.setup_estimation(parameters=["beta0", "beta2"], indexes=[np.array([4]), np.array([0])])
# gm.setup_estimation()
# gm.estimate()
# gm.to_pickle()

#%% Estimation in steps

gm.setup_estimation(parameters=["xi_in", "xi_out"])
gm.estimate()
gm.t_compare_init_to_estimated()
gm.tables_to_pdf(gm.name + " Step 1")
gm.setup_estimation(parameters=["phi0", "phi1", "xi_in", "xi_out", "beta0", "beta1"])
gm.estimate()
gm.t_compare_init_to_estimated()
gm.tables_to_pdf(gm.name.replace("_", " ") + " Step 2")
gm.setup_estimation(parameters=["xi_in", "xi_out", "kappa0", "phi0", "kappa1", "phi1"])
gm.estimate()
gm.t_compare_init_to_estimated()
gm.tables_to_pdf(gm.name.replace("_", " ") + " Step 3")
gm.setup_estimation(parameters=["xi_in", "xi_out", "kappa0", "phi0", "kappa1", "phi1", "beta0", "beta1"])
gm.estimate()
gm.t_compare_init_to_estimated()
gm.tables_to_pdf(gm.name.replace("_", " ") + " Step 4")
gm.setup_estimation(parameters=["xi_in", "xi_out", "kappa0", "phi0", "kappa1", "phi1", "beta0", "sigma", "beta1"])
gm.estimate()
gm.t_compare_init_to_estimated()
gm.tables_to_pdf(gm.name.replace("_", " ") + " Step 5 - All but no beta2")
gm.setup_estimation()
gm.estimate()
gm.t_compare_init_to_estimated()
gm.tables_to_pdf(gm.name + "Step 6 - All")
gm.to_pickle()