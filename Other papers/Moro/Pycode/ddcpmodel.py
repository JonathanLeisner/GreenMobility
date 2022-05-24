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
        self._VF = np.zeros((self.nt+1,self._maxexp,3,2,2))     
        
        # start the object by computing the value function
        self.ValueFunction()

    def nlincome(self,works,edu,mar) :
        """ returns non labor income """
        return exp(self.nli_intercept + works*self.nli_work + self.nli_educ[edu] + (mar*self.nli_married))

    def detwage(self,x,edu,lw) :
        """ returns the deterministic part of the wage"""
        return exp(self.w_intercept + self.w_xper * x + self.w_xpersq*x**2 + self.w_educ[edu] + self.w_lastnw * (1-lw))
    
    def offerProb(self,edu,mar) :
        """ returns the probability of receiving an offer"""
        
        temp = self.offerprob_intc + self.offerprob_educ[edu] + self.offerprob_married*mar
        return exp(temp)/(1+exp(temp))
        
    def currentUtility(self,work,x,edu,lw,mar,wsh,ush) :
        """ returns the utility from the current period"""
        if work==0 :
            return np.log(self.nlincome(0,edu,mar))
        else :
            return np.log(self.nlincome(1,edu,mar)+ self.detwage(x,edu,lw)*wsh) + self.u_work + ush
    
    def _choiceSpecificValue(self,work,time,x,edu,lw,mar,wsh,ush) :
        """ returns the choice specific value function 
            first argument is the choice        
        """
        if (time == self.nt) :
            return self.currentUtility(0,x,edu,lw,mar,0,0)
        else : # if time from 0 to nt-1
            return self.currentUtility(work,x,edu,lw,mar,wsh,ush) + self.beta * self._VF[time+1,x+work,edu,work,mar]
                                             
    def _valueNoOffer(self,time,x,edu,lw,mar,wsh,ush) :
        """ value for those that do not receive an offer """
        return self._choiceSpecificValue(0,time,x,edu,lw,mar,0,0)
            
    def _valueOffer(self,time,x,edu,lw,mar,wsh,ush) :
        """ value for those that receive an offer """
        return max(self._choiceSpecificValue(0,time,x,edu,lw,mar,wsh,ush),
                   self._choiceSpecificValue(1,time,x,edu,lw,mar,wsh,ush))
    
    def _EV_NoOffer(self,time,x,edu,lw,mar,wsh,ush) :
        """ returns expected value if don't receive an offer 
            for a given vector of state variables
        """
        # nothing to integrate here because if no work, shocks don't affect utility
        # keep it anyway if needed in future
        return self._valueNoOffer(time,x,edu,lw,mar,0,0)
                
    def _EV_Offer(self,time,x,edu,lw,mar,ush,wsh) :
        """ expected value if receives an offer
            for a given vector of state variables        
        """
        if time == self.nt :
            return -float('inf')
        else : 
            sum = 0
            for draw in range(self.dt) :
                sum = sum + self._valueOffer(time,x,edu,lw,mar,wsh[draw],ush[draw])
            return sum/self.dt
    
    def ValueFunction(self) :
        """ returns the vector of value function for all states
            and the vector of shocks used to generate them
            list(VFno,VFof,ush,wsh)
        """
        # draw all shocks now
        np.random.seed(1234)
        periods = self.nt + 1
        ush = np.random.normal( 0, exp(self.alpha_cov_u_1_1),size=self.dt*periods).reshape(self.dt,periods)
        wsh = np.random.lognormal(mean=0.0, sigma=sqrt(exp(self.alpha_var_wsh)),size=self.dt*periods).reshape(self.dt,periods)
        
        for t in reversed(range(self.nt+1)) :  # start from last period
            if t == self.nt :
                print('Computing value function for period ',t,end="")
            else :
                print(' ...', t,end="")
            if t == 0 :
                print('')

            for xper in range(t+1) :           # experience can't be more than t
                for ed in range(3) :           # education 0 to 2
                    for lw in range(2) :       # working last period 0,1
                        for mar in range (2) : # married 0,1
                            if xper<=t :  
                                if lw == 1 and t < self.nt : # receives offer for sure
                                    self._VF[t,xper,ed,lw,mar] = self._EV_Offer(t,xper,ed,lw,mar,ush[:,t],wsh[:,t]) 
                                elif t<self.nt : # may receive offer
                                    self._VF[t,xper,ed,lw,mar] = self.offerProb(ed,mar) * self._EV_Offer(t,xper,ed,lw,mar,ush[:,t],wsh[:,t]) + \
                                        (1-self.offerProb(ed,mar)) * self._EV_NoOffer(t,xper,ed,lw,mar,ush[:,t],wsh[:,t])                              
                                else : # last period
                                    self._VF[t,xper,ed,lw,mar] = self._EV_NoOffer(t,xper,ed,lw,mar,ush[:,t],wsh[:,t]) 
        return (self._VF,ush,wsh)
    
    def history(self,edu,mar, drawshock = 1, wshock=float('inf'),ushock=float('inf'),ofshock=float('inf')) :
        """ returns a single history """

        # initialize data
        x      = [0]*(self.nt+1) #np.zeros(self.nt+1,dtype=np.int)
        lw     = [0]*(self.nt+1) #np.zeros(self.nt+1,dtype=np.int)
        work   = [0]*(self.nt+1) #np.zeros(self.nt+1,dtype=np.int)
        owage  = [0]*(self.nt+1) #np.zeros(self.nt+1)
        offer  = [0]*(self.nt+1) #np.zeros(self.nt+1,dtype=np.int)
        time   = [0]*(self.nt+1) #np.zeros(self.nt+1,dtype=np.int)
        
        if drawshock == 1 :
            ushock = np.random.normal( 0, exp(self.alpha_cov_u_1_1),size=self.nt)
            wshock = np.random.lognormal(mean=0.0, sigma=sqrt(exp(self.alpha_var_wsh)),size=self.nt)
            ofshock = np.random.uniform(0,1,size=self.nt)

        for t in range(self.nt+1) :
            time[t] = t         #want to include time in return function
            if t>0 :
                x[t] = x[t-1]+work[t-1]     
                lw[t] = work[t-1]
        
        # now figure out choices
            if t==self.nt : 
                work[t] = 0 #nobody works last period
                owage[t] = -float('inf')
            else :    
            #worked last time or good draw: receives offer and must choose
                if lw[t] == 1 or ofshock[t] <= self.offerProb(edu,mar)  :
                    offer[t] = 1
                    chvalue =   self._choiceSpecificValue(0,t,x[t],edu,lw[t],mar,wshock[t],ushock[t]) \
                               ,self._choiceSpecificValue(1,t,x[t],edu,lw[t],mar,wshock[t],ushock[t])         
                    for choice in range(2) :
                        if chvalue[choice] == max(chvalue) :
                            work[t] = choice
                else : #no offer, no job
                    offer[t] = 0
                    work[t] == 0
                owage[t] =  self.detwage(x[t],edu,lw[t])*wshock[t]
        
        return {'time':time,'work':work,'exp': x,'offer':offer,'owage':owage}
          
    def many_histories(self,rep=1,drawshock=1,ifile='InitData_BasicDPP_1.txt', 
                       wshock=float('inf'),ushock=float('inf'),ofshock=float('inf')):
        """ returns a set of many histories from initial conditions file """
       
        dfhist = []
        observations = []
        
        #load observations from file
        with open(ifile) as initcond :
            obs = csv.reader(initcond)
            rownum = 0
            for row in obs :
                observations.insert(rownum,{'id': row[0],'exper' : row[1]
                                            ,'edu' : int(row[2]), 'lw': row[3],'mar': int(row[4])})
                rownum= rownum + 1

        # create an empty shell for observation data
        allhists = {'id':[],'time':[],'work':[],'edu':[],
                    'exp':[],'mar':[],'offer':[],'owage': []
                    } 

        n = 0
        for ob in observations :
            # replicate observation as requested
            for reps in range(rep) :
                
                # compute a history
                hist = self.history(ob['edu'],ob['mar'],1) #0,ushock[n,:],wshock[n,:],ofshock[n,:]))

                # append this history to the big list
                for elem in hist.keys() :                   
                    allhists[elem].extend(hist[elem])
                    
                # append nt+1 values of the constant variables
                allhists['id'].extend([n]*(self.nt+1))
                allhists['edu'].extend([ob['edu']]*(self.nt+1))
                allhists['mar'].extend([ob['mar']]*(self.nt+1))
                n =n +1
         
        # now create a nice pandas data frame with data 
        reshapeHist =  {
            'id'   : pd.Series(allhists['id']),
            'time' : pd.Series(allhists['time']),
            'exp'  : pd.Series(allhists['exp']),
            'edu'  : pd.Series(allhists['edu']),
            'mar'  : pd.Series(allhists['mar']),
            'work' : pd.Series(allhists['work']),
            'offer':   pd.Series(allhists['offer']),
            'owage':   pd.Series(allhists['owage'])
                }    
        dfhist = pd.DataFrame(reshapeHist)
 
        return dfhist
    
    def saveData(self,df=1,filename='simdata.csv') : 
        '''saves simulated data to csv file , computes histories if not provided
        '''
        
        if type(df) == int :
            df = self.many_histories(rep=10)

        df.to_csv(filename)    
    
    def stats(self,df=1,level=['time']) :
        """ returns some summary statistics"""   
        if type(df) == int :
            df = self.many_histories(rep=10)

        # print some statistics
        if 'time' in level :       
            myformat = '{:3g}{:8.2f}{:10.2f}{:10.2f}'
            print('Time ','%Work ',' Wage','      Offer')
            for t in range(self.nt) :
                print(myformat.format(t
                    ,df["work"][df['time']==t].mean()
                    ,df["owage"][df['work']==1][df['time']==t].mean()
                    ,df["owage"][df['time']==t].mean()        
                    )
                )
            print()
         
        if 'edu' in level :
            myformat = '{:3g}{:4g}{:7.2f}{:10.2f}{:10.2f}'
            print('Time','Edu',' %Work ',' Wage','      Offer')
            for ed in range(3) :
                for t in range(self.nt) :    
                    print(myformat.format(t,ed
                        ,df["work"][df['time']==t][df['edu']==ed].mean()
                        ,df["owage"][df['work']==1][df['time']==t][df['edu']==ed].mean()
                        ,df["owage"][df['time']==t][df['edu']==ed].mean()
                        )        
                    )
            print()
        
        if 'exp' in level :
            myformat = '{:3g}{:4g}{:4g}{:7.2f}{:10.2f}{:10.2f}  {:3g}/{:3g}'
            print('Time','Edu','Exp',' %Work ',' Wage','      Offer','Working/Tot  ')
            for x in range(self._maxexp) :
                for ed in range(3) :
                    for t in range(self.nt) :    
                        if t>=x :
                            print(myformat.format(t,ed,x
                                ,df["work"][df['time']==t][df['edu']==ed][df['exp']==x].mean()                
                                ,df["owage"][df['work']==1][df['time']==t][df['edu']==ed][df['exp']==x].mean()
                                ,df["owage"][df['time']==t][df['edu']==ed][df['exp']==x].mean()
                                ,df["work"][df['time']==t][df['edu']==ed][df['exp']==x][df['work']==1].count()
                                ,df["work"][df['time']==t][df['edu']==ed][df['exp']==x].count()
                            ))
            print()

        if 'reg' in level :
            # A mincerian regression
            dfwork = df[df['work']==1]
            dfwork['x2'] = dfwork['exp']**2
            dfwork['lw'] = np.log(dfwork['owage'])
            dfwork['e1'] = dfwork['edu'] == 1
            dfwork['e2'] = dfwork['edu'] == 2
            result = sm.ols(formula="lw ~ e1 + e2 + exp + x2", data=dfwork).fit()
            print(result.summary())
    


if __name__ == "__main__":
    # stuff only to run when not called via 'import' here
    # set up the ddcp object  

    np.random.seed(1234)
    d2 = ddcpmodel(dt=2500)
#    
#    # compute the value function
#    start_time = time.clock()
#    print('Value function computed in {:5.3f} seconds'.format(time.clock() - start_time))
#    print('')
    
    # generate simulated histories
    hist = d2.many_histories(rep=10)
    
    d2.stats(hist,'time')
    d2._VF[3,3,:,1,0]
    d2.saveData(hist,filename='simdata.csv')
    
    #
    #d2 = ddcp(beta=0)
    #vf2=d2.ValueFunction()
    #hist2 = d2.many_histories()
    #d2.stats(hist2,'reg')




