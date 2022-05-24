# -*- coding: utf-8 -*-
"""
Created on Fri May  6 17:01:23 2016

@author: Andrea Moro, partially ported and adatted from Fortran code by Meghan Skira
"""

import numpy as np
import csv
from math import exp,sqrt
import pandas as pd
import statsmodels.formula.api as sm

class ddcpmodel() :
    """
    This class defines methods and attributes needed to solve a dynamic discrete
    -choice problem
    """
    
    def __init__(self,     
   
            #non-labor income parameters
            nli_intercept =   7.50,
            nli_work =        -0.85,
            nli_educ  =       [0,0.1,0.2],
            nli_married =     0.25,
    
            # earnings parameters
            w_intercept =     10.1,              # Intercept.
            w_xper=           0.05,              # Coefficient on experience.
            w_xpersq =        -0.0006,           # Coefficient on experience squared.
            w_educ =          [0,0.26,0.65],     # Coefficients on 3 levels of education
            w_lastnw =        -0.13,             # Penalty from not working last period.
            
            # (dis)utility from working
            u_work=           -2.65,             
            
            # job offer probability 
            offerprob_intc =  -0.30,
            offerprob_educ =  [0,0.15,0.26],
            offerprob_married = -0.14,
            
            # var-covariance parameters
            alpha_cov_u_1_1 = 0.08,   # Parameter that governs variance of shock to working.
            alpha_var_wsh =  0.40,    # Parameter that governs variance of wage shock.
                                      # Note, this is not the value of the variance itself (see below).                    
    
            # discount factor
            beta = 0.95,
            
            dt = 2500,              # number of monte-carlo integration draws
            dtdraws = 15000 ,       # number of simulation draws
            nt = 4                  # number of time periods -1 (period 0 is first)
        ):

        self.nli_intercept,self.nli_work,self.nli_educ = nli_intercept,nli_work,nli_educ  
        self.nli_married = nli_married
        self.w_intercept,self.w_xper,self.w_xpersq = w_intercept,w_xper,w_xpersq
        self.w_educ,self.w_lastnw = w_educ,w_lastnw
        self.u_work = u_work
        self.offerprob_intc, self.offerprob_educ, self.offerprob_married = offerprob_intc,offerprob_educ,offerprob_married
        self.alpha_cov_u_1_1 = alpha_cov_u_1_1
        self.alpha_var_wsh = alpha_var_wsh
        self.dt,self.dtdraws,self.nt = dt,dtdraws,nt
        self.beta = beta
        
        # additional variables we need later
        self._maxexp = nt+1
        
        # this multidimensional arrays contains values of eqn. 
#??        self._VF = np.zeros((dimensions))     
        
        # start the object by computing the value function
#        self.ValueFunction()

    def nlincome(self,works,edu,mar) :
        """ returns non labor income """
        return exp(self.nli_intercept + works*self.nli_work + self.nli_educ[edu] + (mar*self.nli_married))

    def detwage(self,x,edu,lw) :
        """ returns the deterministic part of the wage"""
        return 0
    def offerProb(self,edu,mar) :
        """ returns the probability of receiving an offer"""
        
        return 0 

        

if __name__ == "__main__":
    # stuff only to run when not called via 'import' here
    # set up the ddcp object  

    np.random.seed(1234)
    d2 = ddcpmodel(dt=2500)
    
#    
#d2.nlincome(1,1,1)
#Out[27]: 1096.6331584284585
#
#d2.nlincome(1,2,1)
#Out[28]: 1211.9670744925775
#
#d2.detwage(1,1,0)
#Out[34]: 29126.39247425565
#
#d2.detwage(1,0,1)
##Out[35]: 25575.752150842934

#d2.offerProb(1,1)
#Out[36]: 0.42800386706848137
#
#d2.offerProb(1,0)
#Out[37]: 0.46257015465625045