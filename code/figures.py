
#%%

import numpy as np
import matplotlib.pyplot as plt
np.set_printoptions(precision=3)

def empty_fig(model):
    fig = plt.figure()
    plt.plot(model.sol.density[0, 0, :])
    fig.savefig("filename.pdf")

# %%
