MODULE ParallelCohorts_MOD

Contains
    
subroutine ParallelCohorts(RE_loop, older_coh, newer_coh, year, param, cut, avec, EducSim_in, GenderSim_in, typ, AllocEmp, rsk, eps, eta, PI_Myopic, PI_RE, &
                                  ChoiceSim_in, CohortWgt_in, ChoiceSim_out, ExperSim_out, WageSim_out, WelfareSim_out, &
                                  emp_out,emp_by_skill_out,sk_sup_out,Welfare_out,Welfare_Educ_out,Real_Wages_out,Real_Wages_Educ_out)   

USE Global_Data
USE Emax_MOD

implicit none

integer, intent(in)             :: RE_loop, &
                                   older_coh, &
                                   newer_coh, &
                                   year, &
                                   EducSim_in(:,older_coh:), &
                                   GenderSim_in(:,older_coh:), &
                                   typ(:,FIRST_COH:), &
                                   ChoiceSim_in(:,older_coh:,-9:)
real(KIND=DOUBLE), intent(in)   :: param(:), &
                                   cut(:,0:), &
                                   avec(:,0:,FIRST_YR_SIM:), &
                                   AllocEmp(:,:,:,:), &
                                   rsk(:,0:), &
                                   eps(:,older_coh:,0:), &
                                   eta(:,older_coh:,0:), &
                                   CohortWgt_in(older_coh:,0:), &
                                   PI_Myopic(FIRST_AGE+1:,0:,:,:,:,:), &
                                   PI_RE(FIRST_AGE+1:,0:,:,:,:,:)
integer, intent(out)            :: ChoiceSim_out(:,older_coh:), &
                                   ExperSim_out(:,older_coh:)
real(KIND=DOUBLE), intent(out)  :: WageSim_out(:,older_coh:), &
                                   sk_sup_out(:,0:), &
                                   WelfareSim_out(:,older_coh:), &
                                   Welfare_out, & 
                                   Welfare_Educ_out(0:), & 
                                   Real_Wages_out, & 
                                   Real_Wages_Educ_out(0:), &
                                   emp_out(0:), &
                                   emp_by_skill_out(0:,0:)

real(KIND=DOUBLE) theta(7), beta(NSECTORS,12), sigma(0:NSECTORS), kappa(26), &
                  tau(0:NSECTORS), SigmaPref, omega(2:NTYPES,0:NSECTORS), lambda(1:NTYPES)

real(KIND=DOUBLE) cost(0:NSECTORS,0:NSECTORS), &
                  cost25(0:NSECTORS,0:NSECTORS), &
                  eps_aux(0:NSECTORS), &
                  eta_aux(0:NSECTORS), &
                  Emax(0:NSECTORS), &
                  rsk_tomorrow(NSECTORS), &
                  w(0:NSECTORS), &
                  V(0:NSECTORS), &
                  VMAX, &
                  PI_COEF(NREG)

integer coh, n, ed, s, i, s1, s2, &
        age, gender, educ, tp, lag, exper(NSECTORS), ExperTomorrow(NSECTORS), &
        educ2, educ3, educ4, &
        dummy(NSECTORS,0:1)

theta(1:7)      = param(1:7)
beta(1,1:12)    = param(8:19)
beta(2,1:12)    = param(20:31)
beta(3,1:12)    = param(32:43)
beta(4,1:12)    = param(44:55)
beta(5,1:12)    = param(56:67)
beta(6,1:12)    = param(68:79)
beta(7,1:12)    = param(80:91)
sigma(0:7)      = param(92:99)
kappa(1:26)     = param(100:125)
tau             = 0.0
tau(2:7)        = param(126:131)
SigmaPref       = param(132)
omega(2,0:7)    = param(133:140)
omega(3,0:7)    = param(141:148)
lambda          = 0.0
lambda(2:3)     = param(149:150)

dummy = 0
do s = 1, NSECTORS
    do ed = 0, 1
        if (rsk(s,ed) >= cut(s,ed)) then
            dummy(s,ed) = 1
        end if
    end do
end do

emp_out = 0
emp_by_skill_out = 0
sk_sup_out = 0
Welfare_out = 0.0
Welfare_Educ_out = 0.0
Real_Wages_out = 0.0
Real_Wages_Educ_out = 0.0

cohorts_loop: do coh = older_coh, newer_coh
                      
    individuals_loop: do n = 1, NPOP_SIM
        
        age     = year - coh 
        gender  = GenderSim_in(n,coh)
        educ    = EducSim_in(n,coh)
        tp      = typ(n,coh)
        lag     = ChoiceSim_in(n,coh,-1)
        eps_aux = eps(n,coh,:)
        eta_aux = eta(n,coh,:)
                
                
        if (educ <= 2) then
            ed = 0
        else if (educ >= 3) then
            ed = 1
        end if
                            
          
        ! experience accumulated in the last 8 years
        exper = 0
        do i = 1, 8
            do s = 1, NSECTORS                    
                if (ChoiceSim_in(n,coh,-i) == s) then
                    exper(s) = exper(s) + 1
                end if                    
            end do
        end do

        if (age < LAST_AGE) then

            if (year >= YR_ANNMNT .and. year <= (LAST_YR_SIM-1) .and. RE_loop >= 2) then
                rsk_tomorrow = avec(:,ed,year+1)*rsk(:,ed)
            else
                rsk_tomorrow = rsk(:,ed)
            end if

            do s = 0, NSECTORS
                
                if (year >= YR_ANNMNT .and. year <= (LAST_YR_SIM-2) .and. RE_loop >= 2) then
                    PI_COEF = PI_RE(age+1,s,gender,educ,tp,:)
                else
                    PI_COEF = PI_Myopic(age+1,s,gender,educ,tp,:)
                end if
                
                ExperTomorrow    = exper
                if (s >= 1) then
                    ExperTomorrow(s) = exper(s) + 1
                end if
                
                Emax(s)          = &
                Emax_hat(PI_COEF,rsk_tomorrow,ExperTomorrow,cut(:,ed))
                
            end do

        else
            Emax = 0.0
        end if
        
                    
        ! experience accumulated in the last 9 years
        do s = 1, NSECTORS
        
            if (ChoiceSim_in(n,coh,-9) == s) then
                exper(s) = exper(s) + 1
            end if
        
        end do                
                    
        educ2 = 0
        educ3 = 0
        educ4 = 0
        if (educ == 2) then
            educ2 = 1
        elseif (educ == 3) then
            educ3 = 1
        elseif (educ == 4) then
            educ4 = 1
        end if
                 
                 
        if (tp == 1) then
            w(0) = exp(theta(1) + theta(2)*(gender-1) + theta(3)*educ2 + &
                       theta(4)*educ3 + theta(5)*educ4 + theta(6)*(age-25) + & 
                       theta(7)*(age-25)**2) + eps_aux(0)
        else if (tp == 2 .or. tp == 3) then
            w(0) = exp(omega(tp,0) + theta(1) + theta(2)*(gender-1) + theta(3)*educ2 + &
                       theta(4)*educ3 + theta(5)*educ4 + theta(6)*(age-25) + & 
                       theta(7)*(age-25)**2) + eps_aux(0)
        end if
            
        do s = 1, NSECTORS
            if (tp == 1) then
                w(s) = rsk(s,ed)*exp(beta(s,1)*(gender-1) + beta(s,2)*educ2 + & 
                beta(s,3)*educ4 + beta(s,4)*(age-25) + beta(s,5)*(age-25)**2 + &
                beta(s,6)*exper(1) + beta(s,7)*exper(2) + beta(s,8)*exper(3) + & 
                beta(s,9)*exper(4) + beta(s,10)*exper(5) + beta(s,11)*exper(6) + beta(s,12)*exper(7))*exp(eps_aux(s))           
            else if (tp == 2 .or. tp == 3) then
                w(s) = rsk(s,ed)*exp(omega(tp,s) + beta(s,1)*(gender-1) + beta(s,2)*educ2 + & 
                beta(s,3)*educ4 + beta(s,4)*(age-25) + beta(s,5)*(age-25)**2 + &
                beta(s,6)*exper(1) + beta(s,7)*exper(2) + beta(s,8)*exper(3) + & 
                beta(s,9)*exper(4) + beta(s,10)*exper(5) + beta(s,11)*exper(6) + beta(s,12)*exper(7))*exp(eps_aux(s))
            end if
        end do
                       
        cost = 0.0                                                             
           
        do s1 = 0, NSECTORS
            do s2 = 1, NSECTORS
                if (s1 == 0) then
                    cost(s1,s2) = exp(kappa(s2) + kappa(21)*(gender-1) + kappa(22)*educ2 + & 
                    kappa(23)*educ3 + kappa(24)*educ4 + kappa(25)*(age-25) + & 
                    kappa(26)*(age-25)**2 + lambda(tp))
                else if (s1 == 1 .and. (s1 /= s2)) then
                    cost(s1,s2) = exp(kappa(13+s2) + kappa(21)*(gender-1) + kappa(22)*educ2 + & 
                    kappa(23)*educ3 + kappa(24)*educ4 + kappa(25)*(age-25) + & 
                    kappa(26)*(age-25)**2 + lambda(tp))
                else if (s1 /= 0 .and. s1 /= 1 .and. (s1 /= s2)) then
                    cost(s1,s2) = exp(kappa(6+s1) + kappa(13+s2) + kappa(21)*(gender-1) + kappa(22)*educ2 + & 
                    kappa(23)*educ3 + kappa(24)*educ4 + kappa(25)*(age-25) + & 
                    kappa(26)*(age-25)**2 + lambda(tp))
                end if                 
            end do
        end do
        
        cost25 = 0.0
        if (age == 25) then
            do s2 = 1, NSECTORS
                do s1 = 1, NSECTORS
                    if (s2 /= s1) then
                        cost25(:,s2) = cost25(:,s2) + AllocEmp(s1,tp,gender,educ)*cost(s1,s2)
                    end if
                end do
            end do
            cost = cost25
        end if
 
        
        ExperSim_out(n,coh) = exper(3)
        

                
        do s = 0, NSECTORS
            V(s) = tau(s) + w(s) - cost(lag,s) + eta_aux(s) + rho*Emax(s)
        end do            
                
            
        VMAX = maxval(V)
  
                
        do s = 0, NSECTORS   
            if (VMAX == V(s)) then

                ChoiceSim_out(n,coh)   = s
                emp_out(s)             = emp_out(s) + CohortWgt_in(coh,ed)
                emp_by_skill_out(s,ed) = emp_by_skill_out(s,ed) + CohortWgt_in(coh,ed)

                if (s > 0) then

                    sk_sup_out(s,ed)        = sk_sup_out(s,ed) + &
                                              CohortWgt_in(coh,ed)*(w(s)/rsk(s,ed))
                    Real_Wages_out          = Real_Wages_out + w(s)*CohortWgt_in(coh,ed)
                    Real_Wages_Educ_out(ed) = Real_Wages_Educ_out(ed) + w(s)*CohortWgt_in(coh,ed)

                    Welfare_out          = Welfare_out + &
                                          (tau(s) + w(s) - cost(lag,s) + eta_aux(s))*CohortWgt_in(coh,ed)
                    Welfare_Educ_out(ed) = Welfare_Educ_out(ed) + &
                                          (tau(s) + w(s) - cost(lag,s) + eta_aux(s))*CohortWgt_in(coh,ed)
                    WelfareSim_out(n,coh) = (tau(s) + w(s) - cost(lag,s) + eta_aux(s))
                    WageSim_out(n,coh) = w(s)

                else

                    Welfare_out = Welfare_out + &
                                 (tau(0) + w(0) + eta_aux(0))*CohortWgt_in(coh,ed)
                    Welfare_Educ_out(ed) = Welfare_Educ_out(ed) + &
                                         (tau(0) + w(0) + eta_aux(0))*CohortWgt_in(coh,ed)
                    WelfareSim_out(n,coh) = tau(0) + w(0) + eta_aux(0)

                end if
                        
            end if
        end do
               
    end do individuals_loop 
            
end do cohorts_loop

end subroutine ParallelCohorts            
                                                                              
end module ParallelCohorts_MOD                                   
                           
                           
