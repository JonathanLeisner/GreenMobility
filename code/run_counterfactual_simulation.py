"""  Defines the class that runs the Green Mobility (GM) project code."""
#%%
import os
os.chdir("I:/Workdata/706651/users/Jonathan/Projects/GreenMobility/code")
import numpy as np
np.set_printoptions(precision=3)

import warnings
warnings.filterwarnings("error", category=RuntimeWarning)
warnings.filterwarnings("ignore", message="numpy.ufunc size has changed")
warnings.filterwarnings("ignore", message="numpy.ufunc size changed, may indicate binary incompatibility. Expected 192 from C header, got 216 from PyObject")

#User-made modules
import counterfactual_simulation
# import importlib
# import _output
# importlib.reload(_output)

#%%

cf = counterfactual_simulation.counterfactual_simulation(from_estimation_pickle="estimated_v4")
cf.setup()
cf.allocate_worker()
cf.allocate_firm()
cf.precompute()
cf.c_future_tau("fixed tax")
cf.update_p()
cf.solve_worker()
cf.setup_simulation()
cf.simulate_endoD()
cf.solve_ED0()

#%% Baseline scenario figures

cf.f_emsint(compare=False, latex="baseline")
cf.f_employment_shares(compare=False, latex="baseline")
cf.f_p(compare=False, latex="baseline")
cf.f_r(compare=False, latex="baseline")
cf.f_Y1z(compare=False, latex="baseline")
cf.f_Hsup(compare=False, latex="baseline")
cf.f_z(compare=False, latex="baseline")
cf.f_Z(compare=False, latex="baseline")
cf.f_tau(compare=False, latex="baseline")
cf.f_avg_realized_wages(compare=False, latex="baseline")
cf.f_K(compare=False, latex="baseline")
cf.figures_to_pdf("Baseline simulation scenario")

#%% Tax scenario. 

cf.update_comparison_scenario()
cf.c_future_tau("rising tax manufacturing", tau_factor_end=2)
cf.solve_ED0()



#%% Test nye figurer
cf.f_K(compare=True, latex=True)
cf.f_emsint(compare=True, latex=True)
cf.f_employment_shares(compare=True, latex=True)
cf.f_p(compare=True, latex=True)
cf.f_r(compare=True, latex=True)
cf.f_Y1z(compare=True, latex=True)
cf.f_Hsup(compare=True, latex=True)
cf.f_z(compare=True, latex=True)
cf.f_Z(compare=True, latex=True)
cf.f_tau(compare=True, latex=True)
cf.f_avg_realized_wages(compare=True, latex=True)
cf.figures_to_pdf("Counterfactual simulations")
cf.t_EV(latex=True)
cf.tables_to_pdf("Counterfactual simulations")

#%%