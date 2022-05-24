    def currentUtility(self,work,x,edu,lw,mar,wsh,ush) :
        """ returns the utility from the current period"""
        if work==0 :
            return ??
        else :
            return ??
    
    def _choiceSpecificValue(self,work,time,x,edu,lw,mar,wsh,ush) :
        """ returns the choice specific value function 
            first argument is the choice        
        """
        equation (8)
                                             
    def _valueNoOffer(self,time,x,edu,lw,mar,wsh,ush) :
        """ value for those that do not receive an offer """
        equation (7)
            
    def _valueOffer(self,time,x,edu,lw,mar,wsh,ush) :
        """ value for those that receive an offer """
        equation (7)
    
    def _EV_NoOffer(self,time,x,edu,lw,mar,wsh,ush) :
        """ returns expected value if don't receive an offer 
            for a given vector of state variables
        """
        # nothing to integrate here because if no work, shocks don't affect utility
        # keep it anyway if needed in future
        EV component of equation (8)
                
    def _EV_Offer(self,time,x,edu,lw,mar,ush,wsh) :
        """ expected value if receives an offer
            for a given vector of state variables        
        """
        EV
