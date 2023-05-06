"""  Defines the class that runs the Green Mobility (GM) project code."""
#%%
import os
os.chdir("I:/Workdata/706651/users/Jonathan/Projects/GreenMobility/code")

#User-made modules
import estimate

#%% Load and test
gm = estimate.GM_estimate(from_pickle="estimated_v4")
data_name = "standard"
gm.c_output_settings()


#Tables
gm.t_compare_init_to_estimated()
gm.t_estimates_nosectordim(latex=True)
gm.t_estimates_sectordim(latex=True)
gm.t_median_switching_costs(latex=True)
gm.t_avg_yearly_transition_rates(data=True, latex=True)
gm.t_avg_yearly_transition_rates(data=False, latex=True)
gm.tables_to_pdf("estimation_" + data_name)

#Figures
gm.f_age_employment_shares(latex=True)
gm.f_employment_shares(latex=True)
gm.f_age_profile_switching(latex=True)
gm.f_avg_realized_wages(latex=True)
gm.figures_to_pdf("estimation_" + data_name)


#%%