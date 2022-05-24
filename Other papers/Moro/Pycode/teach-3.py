
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
        
#        now fill values into self._VF[...]

 