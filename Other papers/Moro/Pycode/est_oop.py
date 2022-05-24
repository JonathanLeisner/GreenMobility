# -*- coding: utf-8 -*-

from ddcpmodel import ddcpmodel
import pandas as pd
import numpy as np
from scipy.optimize import minimize
from math import exp,log

# these two functions are needed for my constrained optimization technique
def squeeze(x,bottom,top) :
    if (x<=700 and x>=-700) :    # stupid way to avoid overflow errors
        squeeze=(top-bottom)*1.2**x/(1.2**x+1)+bottom
    elif (x>700) :
        squeeze=top
    else :
        squeeze=bottom
    return squeeze
    
def blowup(x,bottom,top) :
    if (x<bottom+1.E-10) : # again to avoid overflow errors
        print('a')
        blowup=-700.
    elif (x>top-1.E-10) :
        print('b')
        blowup=700
    else :
        print('c')
        blowup=log((x-bottom)/(top-x))/log(1.2)
    return blowup

class est_oop(ddcpmodel):
    ''' Computes estimates parameters starting from data file
    '''
    
    def __init__(self, simulate = 0,*args,**kwargs) :
        
        ddcpmodel.__init__(self, *args,**kwargs)
        self.datafile='simdata.csv'
        self.simulate = simulate
        if self.simulate == 0 :
            # get the real moments only the first time the object is defined,
            # not every time simulated moments are needed
            self.RealMoments()

    def RealMoments(self) :
        realData = pd.read_csv(self.datafile)
        if self.print_screen>=1:
            print('Now computing "real" moments from ', self.datafile)
        self.realMoments = self.ComputeMoments(df=realData)
        
    def ComputeMoments(self,df=1):
        
        if type(df) == int :
            if (self.print_screen>=2) :
                print('computing simulated histories ...')
            df = self.many_histories(rep=10)
        
        moments = []
        
        #working and wages by education,time,experience,marriage status
        for x in range(self._maxexp) :
            for ed in range(3) :
                for t in range(self.nt) :
                    for marr in range(2) :
                        if t>x :
                            # participation by time, education, experience, marriage status
                            moments.extend([df["work"][df['time']==t][df['edu']==ed][df['exp']==x][df['mar']==marr].mean() ])    
                            # wage of workers by time, education, experience, marriage status
                            moments.extend([df["owage"][df['work']==1][df['time']==t][df['edu']==ed][df['exp']==x][df['mar']==marr].mean() ])
    
        return moments
        
    def LossFunction(self,pars) :

        dd = est_oop(simulate=1,print_screen=1,**pars) #already have real moments, so just simulate
        theseMoments = dd.ComputeMoments(df=1)
        losses = np.subtract(theseMoments,self.realMoments)
        loss = np.inner(losses,losses) #not GMM here because using identity matrix (can improve)
        
        myformat= '{:7.5f} {:7.5f} loss: {:10.5f}'
        if self.print_screen >= 1 :
            print(myformat.format(pars['w_intercept'],pars['w_xper'],loss))
        return loss
    
    # now define the optimizaition functions required by the minimization
    # routine. Will try different cases/examples
    def UnconstrOpt(self,pars) :
        return self.LossFunction({'w_intercept':pars[0]})

    def OptFunction(self,pars) :       
        return self.LossFunction({'w_intercept':squeeze(pars,bottom,top)})
        
    def OptFunction2(self,pars) :       
        return self.LossFunction({'w_intercept':squeeze(pars[0],bottom,top),'w_xper':squeeze(pars[1],bottom1,top1)})

    def UnconstrOpt(self,pars) :
        return self.LossFunction({'w_intercept':pars[0],'w_xper':pars[1]})


# generate some data and save them
aa = est_oop(dt=1000)
aa.saveData()

#next, try to estimate them

import time

# 1 dimensional optimization          
#bottom = 9.5
#top = 10.5
#start_time = time.clock()
#result = minimize(aa.OptFunction,blowup(10.3,bottom,top),method='nelder-mead',options={'disp':True})
#print('x:  ',squeeze(result.x,bottom,top))
#print('Minimization convergence in {:5.3f} seconds'.format(time.clock() - start_time))

# 2 dimensions
#bottom = 9.5
#top = 10.5
#bottom1 = 0.01
#top1 = 0.09
#x0 = [blowup(10.3,bottom,top),blowup(0.02,bottom1,top1)]
#start_time = time.clock()
#result = minimize(aa.OptFunction2,x0,method='nelder-mead',options={'disp':True})
#print('x:  ',squeeze(result.x[0],bottom,top),squeeze(result.x[1],bottom1,top1))
#print('Minimization convergence in {:5.3f} seconds'.format(time.clock() - start_time))

# 2 dimensions, unconstrained
#result = minimize(aa.UnconstrOpt,[10.3,0.2],method='nelder-mead',options={'disp':True})
#print(result.x)

# try a more thorough minimization routine
#from scipy.optimize import basinhopping
#aa.print_screen = 1 # print everything
#x0 = [10.3,0.2]
#start_time = time.clock()
#ret = basinhopping(aa.UnconstrOpt, x0, minimizer_kwargs={'method':'nelder-mead'})
#print('Minimization convergence in {:5.3f} seconds'.format(time.clock() - start_time))
