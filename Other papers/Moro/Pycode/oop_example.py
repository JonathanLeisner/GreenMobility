# -*- coding: utf-8 -*-
import numpy as np

theta = 2

data = [1.3,3.2,2.5,1.1,5.0]

def contrlikelihood(x,theta) :
    return 1/(np.sqrt(2*np.pi)) * np.exp(-(x-theta)**2/2)

likelihood = 0
for datum in data:
    likelihood += contrlikelihood(datum,theta)
    
print(likelihood)



class mymodel :
    def __init__(self,th = 2):
        self.th = th
    
    def contrlikelihood(self,x):
        return 1/(np.sqrt(2*np.pi)) * np.exp(-(x-self.th)**2/2)
    
    def likelihood(self,dd) :
        likelihood = 0
        for datum in dd:
            likelihood += self.contrlikelihood(datum)
        return likelihood
    
    
d1 = mymodel()
d1.contrlikelihood(3.2)
    
    
    
        
    
