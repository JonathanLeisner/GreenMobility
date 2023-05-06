"""  Calibration of emission parameters."""
#%%
import os
os.chdir("I:/Workdata/706651/users/Jonathan/Projects/GreenMobility/code")
import numpy as np
import pandas as pd
import python_to_latex as pytex
import matplotlib.pyplot as plt
from scipy import optimize


#%%

def g_funcs(gamma, form):
    if form == "quadratic":
        g = lambda z: gamma[0] * (gamma[1] + gamma[2] * (1 - z) + gamma[3] * (1 - z)**2)
        gprime = lambda z: gamma[0] * (- gamma[2] - 2 * gamma[3] * (1 - z))
        zopt = lambda x: 1 - (x / gamma[0] - gamma[2]) * (1 / (2 * gamma[3]))
    elif form == "quadratic_parameterized":
        g = lambda z: gamma[0] * (gamma[1] + gamma[2] * (1 - z) + gamma[3] * (1 - z)**(gamma[4] + 1))
        gprime = lambda z: gamma[0] * (- gamma[2] - (gamma[4] + 1) * gamma[3] * (1 - z)**gamma[4])
        zopt = lambda x: 1 - ((x / gamma[0] - gamma[2]) * (1/(gamma[4] + 1)) * (1/gamma[3]))**(1/gamma[4])
    return (g, gprime, zopt)

def n_params(form):
    if form == "quadratic":
        return 4
    elif form == "quadratic_parameterized":
        return 5

#Container for results (init_z and gamma)
form = "quadratic_parameterized"
init_z = np.zeros(shape=np.load("../data_prep/output/pY1z.npy")[-1, :].shape)
gamma = np.zeros(shape=(n_params(form), len(np.load("../data_prep/output/pY1z.npy")[-1, :])))

#%% Choose sector to calibrate for
for s in range(len(init_z)):
    # s = 1
    tau0 = np.load("../data_prep/output/tau.npy")[-1, s]
    etarget = np.load("../data_prep/output/e.npy")[-1, s]
    tau_levels = np.array([tau0, tau0*2])
    e_targets = np.array([etarget, etarget*(1-0.35)])

    def def_objfunc(form, tau_levels, e_targets):
        def objfunc(gamma):
            # def_values = np.array([1, 1, 1, 1, 1])[:len(gamma)]
            # if gamma[4] < -1:
            #     print(f"gamma_4 is below -1, with value {gamma[4]}")
            (g, _, zopt) = g_funcs(gamma, form)
            targets_sum = 0
            for i, tau in enumerate(tau_levels):
                targets_sum += np.square(g(zopt(1/tau))/(1-zopt(1/tau)) - e_targets[i])
            # gamma_sum = np.sum(np.square(gamma - def_values))
            # return gamma_sum + 10 * targets_sum
            return targets_sum
        return objfunc

    objfunc = def_objfunc(form, tau_levels, e_targets)

    def constraints_funcs(tau_levels, e_targets, form):
        def z_under_limit0(gamma):
            (_, _, zopt) = g_funcs(gamma, form)
            return 1 - zopt(1/tau_levels[0])
        def z_under_limit1(gamma):
            (_, _, zopt) = g_funcs(gamma, form)
            return 1 - zopt(1/tau_levels[1])
        def z_over_limit0(gamma):
            (_, _, zopt) = g_funcs(gamma, form)
            return -0.01 + zopt(1/tau_levels[0])
        def init_emission_intensity(gamma):
            (g, _, zopt) = g_funcs(gamma, form)
            return g(zopt(1/tau_levels[0]))/(1-zopt(1/tau_levels[0])) - e_targets[0]
        def reduc_target1(gamma):
            (g, _, zopt) = g_funcs(gamma, form)
            return g(zopt(1/tau_levels[1]))/(1-zopt(1/tau_levels[1])) - e_targets[1]
        if form == "quadratic_parameterized":
            def gamma4_over_limit(gamma):
                return gamma[4] + 1
        # return init_emission_intensity
        return (init_emission_intensity, reduc_target1, z_under_limit0, z_under_limit1, z_over_limit0, gamma4_over_limit)

    def constraints_gprime(form, tau_levels):
        def gprime_negative_init(gamma):
            (_, gprime, zopt) = g_funcs(gamma, form)
            return - gprime(zopt(tau_levels[0]))
        def gprime_negative_end(gamma):
            (_, gprime, zopt) = g_funcs(gamma, form)
            return - gprime(zopt(tau_levels[1]))
        def gprime_negative_99(gamma):
            (_, gprime, _) = g_funcs(gamma, form)
            return - gprime(0.99)
        return (gprime_negative_init, gprime_negative_end)

    (init_emission_intensity, reduc_target1, z_under_limit0, z_under_limit1, z_over_limit0, gamma4_over_limit) = \
        constraints_funcs(tau_levels, e_targets, form=form)
    (gprime_negative_init, gprime_negative_end) = constraints_gprime(form=form, tau_levels=tau_levels)

    gamma0 = np.array([np.mean(e_targets), 93, -34, 133, 66])
    objfunc(gamma0)
    # (g, gprime, zopt) = g_funcs(gamma0, form)

    cons_list = [{"type":"ineq", "fun":z_under_limit0},
                {"type":"ineq", "fun":z_over_limit0},
                {"type":"ineq", "fun":gprime_negative_init},
                {"type":"ineq", "fun":gprime_negative_end}]

    res = optimize.minimize(objfunc, gamma0, method="SLSQP", constraints=cons_list, options={"maxiter":10000, "disp":True, "ftol":1e-14})
    gamma[:, s] = res.x
    (g, _, zopt) = g_funcs(res.x, form=form)
    init_z[s] = zopt(1/tau0)

#%%Testing
# e_targets

#%% Post calibration save
np.save("../data_prep/output/init_z", init_z)
np.save("../data_prep/output/gamma", gamma)

#%% Produce summary table of calibrated gamma parameters

gamma = np.load("../data_prep/output/gamma.npy")
df = pd.DataFrame(gamma, 
                  index=[f"$\\gamma^s_{i}$" for i in range(gamma.shape[0])], 
                  columns=["Food", "Manufacturing", "Energy", "Transport", "Services"]).reset_index().rename(columns={"index":" "})
note = "The table reports the calibrated values of the sector-specific abatement function parameters." + \
    " All parameters are calibrated such that" + \
    " the squared deviations of the model's emission intensities to the corresponding data" + \
    " counterparts, see \\cref{eq:squared deviations calibration}, are minimized. The hypothetical calibration targets" +\
    " are chosen to be $(\\tau^s_0, EI^s_0)$ and $(3 \\times \\tau^s_0, 0.75 \\times EI^s_0)$ where $\\tau^s_0$ and $EI^s_0$" +\
    " are initial values of the tax and the emission intensity calculated from sector-level data."
label = "calibrated_gamma"
pytex.write_tbls_to_latex([pytex.TexTable(df, caption="Calibrated abatement function parameters", label=label, decimals=4, note=note)], 
                          "gamma parameter values")

#%% Post calibration testing
form = "quadratic_parameterized"
s = 1
tau0 = np.load("../data_prep/output/tau.npy")[-1, s]
etarget = np.load("../data_prep/output/e.npy")[-1, s]
tau_levels = np.array([tau0, tau0*2])
e_targets = np.array([etarget, etarget*(1-0.35)])
(g, gprime, zopt) = g_funcs(gamma[:, s], form=form)

# (g, gprime, zopt) = g_funcs(gamma0, form)
#Emission intensity
f = lambda tau: g(zopt(1/tau))/(1-zopt(1/tau))
h = lambda tau: g(zopt(1/tau))

Delta_tau_vals = np.linspace(-0.2, 1.4)
plt.plot(Delta_tau_vals, f(tau_levels[0] + tau_levels[0] * Delta_tau_vals))
plt.plot(Delta_tau_vals, h(tau_levels[0] + tau_levels[0] * Delta_tau_vals))
# plt.plot(Delta_tau_vals, f(tau_levels[0] + tau_levels[0] * Delta_tau_vals) / f(tau_levels[0]))
#%%
# plt.axhline(e_targets[0])
# plt.axhline(e_targets[1])

#g as a function of tau
f = lambda tau: g(zopt(1/tau))
plt.plot(f(tau_levels[0] + tau_levels[0] * np.linspace(0, 3)))

#Optimal z
plt.plot(zopt(1/(tau_levels[0] + tau_levels[0] * np.linspace(0, 3))))

#gprime
f = lambda tau: gprime(zopt(1/tau))
plt.plot(f(tau_levels[0] + tau_levels[0] * np.linspace(0, 3)))

#gprime (between 0 and 1)
plt.plot(gprime(np.linspace(0, 1)))

#g (z: between 0 and 1)
plt.plot(g(np.linspace(0, 1)))


#
plt.plot(g(zopt(1 / (tau_levels[0] + tau_levels[0] * np.linspace(0, 3)))))
plt.plot(g(zopt(1 / (tau_levels[0] + tau_levels[0] * np.linspace(0, 3)))))
plt.plot(zopt(1 / (tau_levels[0] + tau_levels[0] * np.linspace(0, 10))))


plt.plot(gprime(np.linspace(0, 1)))





