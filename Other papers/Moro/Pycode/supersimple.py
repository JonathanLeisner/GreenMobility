# -*- coding: utf-8 -*-
"""
Created on Tue Jul  5 16:57:57 2016

@author: andreamoro
"""

import numpy as np
import pandas as pd
from scipy.stats import norm
from scipy.optimize import minimize
from math import log

#### generated the data
# parameters
fColl = 1.5
fHigh = 1
sigma = .35

# other parameters
samplesize = 10000

# initialize seed
np.random.seed(1234)

# draw from normal distribution the two epsilon draws
draw1 = np.random.normal(0,sigma,size = samplesize)
draw2 = np.random.normal(0,sigma,size = samplesize)

#initialize vectors of choices and wages
choice = []
wage = []

# construct the data with the given parameter
for i in range(samplesize) :
    if fHigh + draw1[i] > fColl:
        choice.extend([0])
        wage.extend([fHigh + draw1[i]])
    else :
        choice.extend([1])
        wage.extend([fColl + draw2[i]])
        
# put everything into a pandas data frame for nice computations
df = pd.DataFrame({'wage':wage,'choice':choice})

# compute some statistics
print('Fraction choosing schooling: ',df.choice.mean())
print('Mean wage of low schooling: ',df[df.choice==0].wage.mean())
print('Mean wage of high schooling: ',df[df.choice==1].wage.mean())


# plot the conditional densities
df.wage.plot(kind='density')

df.wage[df.choice==1].plot(kind='density')
df.wage[df.choice==0].plot(kind='density')
        
def contribLik(choice,wage,pars) :
    ''' computes contribution to the likelihood of one observation, given
        choice = {1,0}
        wage = float
        pars = list of 3 parameter values
    '''
    
    fC = pars[0]
    fH = pars[1]
    sig = pars[2]
    
    if choice == 0 :
        cL =1/sig* norm.pdf((wage-fH)/sig)#*(1/(1-norm.cdf((fC-fH)/sig)))
    else :
        cL =1/sig* norm.pdf((wage-fC)/sig)   * norm.cdf((fC-fH)/sig)
    return cL

def logLik(pars):
    loglikelihood = 0
    for i in range(samplesize) :
        try: 
            loglikelihood = loglikelihood + log(contribLik(choice[i],wage[i],pars))
        except:
            print('exception',choice[i],wage[i],pars)
    myformat = '{:5.3f} {:5.3f} {:5.3f} {:8.5f}'
    print(myformat.format(*pars,loglikelihood))
    
    # take the negative if you want to minimize
    return -loglikelihood

parv = [1.5,1,1]
print(logLik(parv))

#sigmas = [x / 200 + 0.5 - 0.1  for x in list(range(30))]
#
#liks = [logLik([1.5,1,x]) for x in sigmas]
#
#pd.DataFrame(liks,sigmas).plot()

#result = minimize(logLik,[1.3,1.1,0.3],method='nelder-mead',options={'disp':True})
print(result.x)

samplesize = 1000000
ddd = np.random.normal(0,sigma,size = samplesize)
dfd = pd.DataFrame(ddd)
print(dfd[dfd<=0.26].count()/samplesize)
print(norm.cdf(0.26/sigma))

dfd.mean()

samplesize = 100000
aaa = np.random.uniform(0,5,samplesize)
b = aaa*3
type(b)
np.sum(b)/samplesize
