MODULE ParallelCohorts_MOD

Contains
    
subroutine ParallelCohorts_Myopic(older_coh, newer_coh, year, param, cut, typ, rsk, eps, eta, PI_Myopic, &
                           ChoiceSim_in, ChoiceSim_out, ExperSim_out, WageSim_out, Cost1Sim_out, Cost2Sim_out, VSim_out, CohortWgt_out, emp_out, sk_sup_out)           

USE Global_Data

implicit none

integer, intent(in)             :: older_coh, &
                                   newer_coh, &
                                   year, &
                                   typ(:,FIRST_COH:), &
                                   ChoiceSim_in(:,older_coh:,year-9:)
real(KIND=DOUBLE), intent(in)   :: param(:), &
                                   cut(:,0:), &
                                   rsk(:,0:), &
                                   eps(:,older_coh:,year:,0:), &
                                   eta(:,older_coh:,year:,0:), &
                                   PI_Myopic(FIRST_AGE+1:,0:,:,:,:,:)
integer, intent(out)            :: ChoiceSim_out(:,older_coh:,year:), &
                                   ExperSim_out(:,:,older_coh:,year:)
real(KIND=DOUBLE), intent(out)  :: WageSim_out(:,older_coh:,year:), &
                                   Cost1Sim_out(:,older_coh:,year:), &
                                   Cost2Sim_out(:,older_coh:,year:), &
                                   VSim_out(:,older_coh:,year:), &
                                   CohortWgt_out(older_coh:,0:), &
                                   sk_sup_out(:,0:)
integer, intent(out)            :: emp_out(0:,0:)

real(KIND=DOUBLE) theta(7), beta(NSECTORS,12), sigma(0:NSECTORS), kappa(26), &
                  tau(0:NSECTORS), SigmaPref, omega(2:NTYPES,0:NSECTORS), lambda(1:NTYPES)

real(KIND=DOUBLE) cost(0:NSECTORS,0:NSECTORS), &
                  eps_aux(0:NSECTORS), &
                  eta_aux(0:NSECTORS), &
                  Emax(0:NSECTORS), &
                  w(0:NSECTORS), &
                  V(0:NSECTORS), &
                  VMAX

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
sk_sup_out = 0

cohorts_loop: do coh = older_coh, newer_coh
        
            
    CohortWgt_out(coh,0) =  & 
    ExpansionFactorGlobal*real(CohortSizeData(coh,0)) / real(CohortSizeData(FIRST_COH,0))
    CohortWgt_out(coh,1) =  & 
    ExpansionFactorGlobal*real(CohortSizeData(coh,1)) / real(CohortSizeData(FIRST_COH,0)) 

            
    individuals_loop: do n = 1, NPOP
        
        age     = year - coh 
        gender  = GenderData(n,coh)
        educ    = EducData(n,coh)
        tp      = typ(n,coh)
        lag     = ChoiceSim_in(n,coh,coh+age-1)
        eps_aux = eps(n,coh,year,:)
        eta_aux = eta(n,coh,year,:)

        if (educ <= 2) then
            ed = 0
        else if (educ >= 3) then
            ed = 1
        end if
        
        exper = 0
        do i = 1, 8
            do s = 1, NSECTORS                    
                if (ChoiceSim_in(n,coh,year-i) == s) then
                    exper(s) = exper(s) + 1
                end if                    
            end do
        end do

        if (age < LAST_AGE) then
                           
            do s = 0, NSECTORS    

                ExperTomorrow    = exper
                if (s >= 1) then
                    ExperTomorrow(s) = exper(s) + 1
                end if 
                                
                Emax(s) = &
                PI_Myopic(age+1,s,gender,educ,tp,1)*1.0 + &
                                
                PI_Myopic(age+1,s,gender,educ,tp,2)*rsk(1,ed) + &
                PI_Myopic(age+1,s,gender,educ,tp,3)*rsk(2,ed) + &
                PI_Myopic(age+1,s,gender,educ,tp,4)*rsk(3,ed) + &
                PI_Myopic(age+1,s,gender,educ,tp,5)*rsk(4,ed) + &
                PI_Myopic(age+1,s,gender,educ,tp,6)*rsk(5,ed) + &
                PI_Myopic(age+1,s,gender,educ,tp,7)*rsk(6,ed) + &
                PI_Myopic(age+1,s,gender,educ,tp,8)*rsk(7,ed) + &
                                
                PI_Myopic(age+1,s,gender,educ,tp,9)*dummy(1,ed)*(rsk(1,ed)-cut(1,ed)) + &
                PI_Myopic(age+1,s,gender,educ,tp,10)*dummy(2,ed)*(rsk(2,ed)-cut(2,ed)) + &
                PI_Myopic(age+1,s,gender,educ,tp,11)*dummy(3,ed)*(rsk(3,ed)-cut(3,ed)) + &
                PI_Myopic(age+1,s,gender,educ,tp,12)*dummy(4,ed)*(rsk(4,ed)-cut(4,ed)) + &
                PI_Myopic(age+1,s,gender,educ,tp,13)*dummy(5,ed)*(rsk(5,ed)-cut(5,ed)) + &
                PI_Myopic(age+1,s,gender,educ,tp,14)*dummy(6,ed)*(rsk(6,ed)-cut(6,ed)) + &
                PI_Myopic(age+1,s,gender,educ,tp,15)*dummy(7,ed)*(rsk(7,ed)-cut(7,ed)) + &
                                
                PI_Myopic(age+1,s,gender,educ,tp,16)*ExperTomorrow(1) + &
                PI_Myopic(age+1,s,gender,educ,tp,17)*ExperTomorrow(2) + &
                PI_Myopic(age+1,s,gender,educ,tp,18)*ExperTomorrow(3) + &
                PI_Myopic(age+1,s,gender,educ,tp,19)*ExperTomorrow(4) + &
                PI_Myopic(age+1,s,gender,educ,tp,20)*ExperTomorrow(5) + &
                PI_Myopic(age+1,s,gender,educ,tp,21)*ExperTomorrow(6) + &
                PI_Myopic(age+1,s,gender,educ,tp,22)*ExperTomorrow(7) + &
                               
                PI_Myopic(age+1,s,gender,educ,tp,23)*rsk(1,ed)**2 + &
                PI_Myopic(age+1,s,gender,educ,tp,24)*rsk(2,ed)**2 + &
                PI_Myopic(age+1,s,gender,educ,tp,25)*rsk(3,ed)**2 + &
                PI_Myopic(age+1,s,gender,educ,tp,26)*rsk(4,ed)**2 + &
                PI_Myopic(age+1,s,gender,educ,tp,27)*rsk(5,ed)**2 + &
                PI_Myopic(age+1,s,gender,educ,tp,28)*rsk(6,ed)**2 + &
                PI_Myopic(age+1,s,gender,educ,tp,29)*rsk(7,ed)**2 + &
                               
                PI_Myopic(age+1,s,gender,educ,tp,30)*dummy(1,ed)*(rsk(1,ed)-cut(1,ed))**2 + &
                PI_Myopic(age+1,s,gender,educ,tp,31)*dummy(2,ed)*(rsk(2,ed)-cut(2,ed))**2 + &
                PI_Myopic(age+1,s,gender,educ,tp,32)*dummy(3,ed)*(rsk(3,ed)-cut(3,ed))**2 + &
                PI_Myopic(age+1,s,gender,educ,tp,33)*dummy(4,ed)*(rsk(4,ed)-cut(4,ed))**2 + &
                PI_Myopic(age+1,s,gender,educ,tp,34)*dummy(5,ed)*(rsk(5,ed)-cut(5,ed))**2 + &
                PI_Myopic(age+1,s,gender,educ,tp,35)*dummy(6,ed)*(rsk(6,ed)-cut(6,ed))**2 + &
                PI_Myopic(age+1,s,gender,educ,tp,36)*dummy(7,ed)*(rsk(7,ed)-cut(7,ed))**2 + &
                               
                PI_Myopic(age+1,s,gender,educ,tp,37)*ExperTomorrow(1)**2 + &
                PI_Myopic(age+1,s,gender,educ,tp,38)*ExperTomorrow(2)**2 + &
                PI_Myopic(age+1,s,gender,educ,tp,39)*ExperTomorrow(3)**2 + &
                PI_Myopic(age+1,s,gender,educ,tp,40)*ExperTomorrow(4)**2 + &
                PI_Myopic(age+1,s,gender,educ,tp,41)*ExperTomorrow(5)**2 + &
                PI_Myopic(age+1,s,gender,educ,tp,42)*ExperTomorrow(6)**2 + &
                PI_Myopic(age+1,s,gender,educ,tp,43)*ExperTomorrow(7)**2 + &
                               
                PI_Myopic(age+1,s,gender,educ,tp,44)*rsk(1,ed)*rsk(2,ed) + &
                PI_Myopic(age+1,s,gender,educ,tp,45)*rsk(1,ed)*rsk(3,ed) + &
                PI_Myopic(age+1,s,gender,educ,tp,46)*rsk(1,ed)*rsk(4,ed) + &
                PI_Myopic(age+1,s,gender,educ,tp,47)*rsk(1,ed)*rsk(5,ed) + &
                PI_Myopic(age+1,s,gender,educ,tp,48)*rsk(1,ed)*rsk(6,ed) + &
                PI_Myopic(age+1,s,gender,educ,tp,49)*rsk(1,ed)*rsk(7,ed) + &
                                
                PI_Myopic(age+1,s,gender,educ,tp,50)*rsk(2,ed)*rsk(3,ed) + &
                PI_Myopic(age+1,s,gender,educ,tp,51)*rsk(2,ed)*rsk(4,ed) + &
                PI_Myopic(age+1,s,gender,educ,tp,52)*rsk(2,ed)*rsk(5,ed) + &
                PI_Myopic(age+1,s,gender,educ,tp,53)*rsk(2,ed)*rsk(6,ed) + &
                PI_Myopic(age+1,s,gender,educ,tp,54)*rsk(2,ed)*rsk(7,ed) + &
                               
                PI_Myopic(age+1,s,gender,educ,tp,55)*rsk(3,ed)*rsk(4,ed) + &
                PI_Myopic(age+1,s,gender,educ,tp,56)*rsk(3,ed)*rsk(5,ed) + &
                PI_Myopic(age+1,s,gender,educ,tp,57)*rsk(3,ed)*rsk(6,ed) + &
                PI_Myopic(age+1,s,gender,educ,tp,58)*rsk(3,ed)*rsk(7,ed) + &
                               
                PI_Myopic(age+1,s,gender,educ,tp,59)*rsk(4,ed)*rsk(5,ed) + &
                PI_Myopic(age+1,s,gender,educ,tp,60)*rsk(4,ed)*rsk(6,ed) + &
                PI_Myopic(age+1,s,gender,educ,tp,61)*rsk(4,ed)*rsk(7,ed) + &
                                
                PI_Myopic(age+1,s,gender,educ,tp,62)*rsk(5,ed)*rsk(6,ed) + &
                PI_Myopic(age+1,s,gender,educ,tp,63)*rsk(5,ed)*rsk(7,ed) + &
                                
                PI_Myopic(age+1,s,gender,educ,tp,64)*rsk(6,ed)*rsk(7,ed) + &
                               
                PI_Myopic(age+1,s,gender,educ,tp,65)*rsk(1,ed)*ExperTomorrow(1) + &
                PI_Myopic(age+1,s,gender,educ,tp,66)*rsk(1,ed)*ExperTomorrow(2) + &
                PI_Myopic(age+1,s,gender,educ,tp,67)*rsk(1,ed)*ExperTomorrow(3) + &
                PI_Myopic(age+1,s,gender,educ,tp,68)*rsk(1,ed)*ExperTomorrow(4) + &
                PI_Myopic(age+1,s,gender,educ,tp,69)*rsk(1,ed)*ExperTomorrow(5) + &
                PI_Myopic(age+1,s,gender,educ,tp,70)*rsk(1,ed)*ExperTomorrow(6) + &
                PI_Myopic(age+1,s,gender,educ,tp,71)*rsk(1,ed)*ExperTomorrow(7) + &
                                
                PI_Myopic(age+1,s,gender,educ,tp,72)*rsk(2,ed)*ExperTomorrow(1) + &
                PI_Myopic(age+1,s,gender,educ,tp,73)*rsk(2,ed)*ExperTomorrow(2) + &
                PI_Myopic(age+1,s,gender,educ,tp,74)*rsk(2,ed)*ExperTomorrow(3) + &
                PI_Myopic(age+1,s,gender,educ,tp,75)*rsk(2,ed)*ExperTomorrow(4) + &
                PI_Myopic(age+1,s,gender,educ,tp,76)*rsk(2,ed)*ExperTomorrow(5) + &
                PI_Myopic(age+1,s,gender,educ,tp,77)*rsk(2,ed)*ExperTomorrow(6) + &
                PI_Myopic(age+1,s,gender,educ,tp,78)*rsk(2,ed)*ExperTomorrow(7) + &
                               
                PI_Myopic(age+1,s,gender,educ,tp,79)*rsk(3,ed)*ExperTomorrow(1) + &
                PI_Myopic(age+1,s,gender,educ,tp,80)*rsk(3,ed)*ExperTomorrow(2) + &
                PI_Myopic(age+1,s,gender,educ,tp,81)*rsk(3,ed)*ExperTomorrow(3) + &
                PI_Myopic(age+1,s,gender,educ,tp,82)*rsk(3,ed)*ExperTomorrow(4) + &
                PI_Myopic(age+1,s,gender,educ,tp,83)*rsk(3,ed)*ExperTomorrow(5) + &
                PI_Myopic(age+1,s,gender,educ,tp,84)*rsk(3,ed)*ExperTomorrow(6) + &
                PI_Myopic(age+1,s,gender,educ,tp,85)*rsk(3,ed)*ExperTomorrow(7) + &
                                
                PI_Myopic(age+1,s,gender,educ,tp,86)*rsk(4,ed)*ExperTomorrow(1) + &
                PI_Myopic(age+1,s,gender,educ,tp,87)*rsk(4,ed)*ExperTomorrow(2) + &
                PI_Myopic(age+1,s,gender,educ,tp,88)*rsk(4,ed)*ExperTomorrow(3) + &
                PI_Myopic(age+1,s,gender,educ,tp,89)*rsk(4,ed)*ExperTomorrow(4) + &
                PI_Myopic(age+1,s,gender,educ,tp,90)*rsk(4,ed)*ExperTomorrow(5) + &
                PI_Myopic(age+1,s,gender,educ,tp,91)*rsk(4,ed)*ExperTomorrow(6) + &
                PI_Myopic(age+1,s,gender,educ,tp,92)*rsk(4,ed)*ExperTomorrow(7) + &
                                
                PI_Myopic(age+1,s,gender,educ,tp,93)*rsk(5,ed)*ExperTomorrow(1) + &
                PI_Myopic(age+1,s,gender,educ,tp,94)*rsk(5,ed)*ExperTomorrow(2) + &
                PI_Myopic(age+1,s,gender,educ,tp,95)*rsk(5,ed)*ExperTomorrow(3) + &
                PI_Myopic(age+1,s,gender,educ,tp,96)*rsk(5,ed)*ExperTomorrow(4) + &
                PI_Myopic(age+1,s,gender,educ,tp,97)*rsk(5,ed)*ExperTomorrow(5) + &
                PI_Myopic(age+1,s,gender,educ,tp,98)*rsk(5,ed)*ExperTomorrow(6) + &
                PI_Myopic(age+1,s,gender,educ,tp,99)*rsk(5,ed)*ExperTomorrow(7) + &
                               
                PI_Myopic(age+1,s,gender,educ,tp,100)*rsk(6,ed)*ExperTomorrow(1) + &
                PI_Myopic(age+1,s,gender,educ,tp,101)*rsk(6,ed)*ExperTomorrow(2) + &
                PI_Myopic(age+1,s,gender,educ,tp,102)*rsk(6,ed)*ExperTomorrow(3) + &
                PI_Myopic(age+1,s,gender,educ,tp,103)*rsk(6,ed)*ExperTomorrow(4) + &
                PI_Myopic(age+1,s,gender,educ,tp,104)*rsk(6,ed)*ExperTomorrow(5) + &
                PI_Myopic(age+1,s,gender,educ,tp,105)*rsk(6,ed)*ExperTomorrow(6) + &
                PI_Myopic(age+1,s,gender,educ,tp,106)*rsk(6,ed)*ExperTomorrow(7) + &
                                
                PI_Myopic(age+1,s,gender,educ,tp,107)*rsk(7,ed)*ExperTomorrow(1) + &
                PI_Myopic(age+1,s,gender,educ,tp,108)*rsk(7,ed)*ExperTomorrow(2) + &
                PI_Myopic(age+1,s,gender,educ,tp,109)*rsk(7,ed)*ExperTomorrow(3) + &
                PI_Myopic(age+1,s,gender,educ,tp,110)*rsk(7,ed)*ExperTomorrow(4) + &
                PI_Myopic(age+1,s,gender,educ,tp,111)*rsk(7,ed)*ExperTomorrow(5) + &
                PI_Myopic(age+1,s,gender,educ,tp,112)*rsk(7,ed)*ExperTomorrow(6) + &
                PI_Myopic(age+1,s,gender,educ,tp,113)*rsk(7,ed)*ExperTomorrow(7) + &
                                
                PI_Myopic(age+1,s,gender,educ,tp,114)*ExperTomorrow(1)*ExperTomorrow(2) + &
                PI_Myopic(age+1,s,gender,educ,tp,115)*ExperTomorrow(1)*ExperTomorrow(3) + &
                PI_Myopic(age+1,s,gender,educ,tp,116)*ExperTomorrow(1)*ExperTomorrow(4) + &
                PI_Myopic(age+1,s,gender,educ,tp,117)*ExperTomorrow(1)*ExperTomorrow(5) + &
                PI_Myopic(age+1,s,gender,educ,tp,118)*ExperTomorrow(1)*ExperTomorrow(6) + &
                PI_Myopic(age+1,s,gender,educ,tp,119)*ExperTomorrow(1)*ExperTomorrow(7) + &
                PI_Myopic(age+1,s,gender,educ,tp,120)*ExperTomorrow(2)*ExperTomorrow(3) + &
                PI_Myopic(age+1,s,gender,educ,tp,121)*ExperTomorrow(2)*ExperTomorrow(4) + &
                PI_Myopic(age+1,s,gender,educ,tp,122)*ExperTomorrow(2)*ExperTomorrow(5) + &
                PI_Myopic(age+1,s,gender,educ,tp,123)*ExperTomorrow(2)*ExperTomorrow(6) + &
                PI_Myopic(age+1,s,gender,educ,tp,124)*ExperTomorrow(2)*ExperTomorrow(7) + &
                PI_Myopic(age+1,s,gender,educ,tp,125)*ExperTomorrow(3)*ExperTomorrow(4) + &
                PI_Myopic(age+1,s,gender,educ,tp,126)*ExperTomorrow(3)*ExperTomorrow(5) + &
                PI_Myopic(age+1,s,gender,educ,tp,127)*ExperTomorrow(3)*ExperTomorrow(6) + &
                PI_Myopic(age+1,s,gender,educ,tp,128)*ExperTomorrow(3)*ExperTomorrow(7) + &
                PI_Myopic(age+1,s,gender,educ,tp,129)*ExperTomorrow(4)*ExperTomorrow(5) + &
                PI_Myopic(age+1,s,gender,educ,tp,130)*ExperTomorrow(4)*ExperTomorrow(6) + &
                PI_Myopic(age+1,s,gender,educ,tp,131)*ExperTomorrow(4)*ExperTomorrow(7) + &
                PI_Myopic(age+1,s,gender,educ,tp,132)*ExperTomorrow(5)*ExperTomorrow(6) + &
                PI_Myopic(age+1,s,gender,educ,tp,133)*ExperTomorrow(5)*ExperTomorrow(7) + &
                PI_Myopic(age+1,s,gender,educ,tp,134)*ExperTomorrow(6)*ExperTomorrow(7)

            end do
                       
        else
            Emax = 0.0
        end if


        ! experience accumulated in the last 9 years
        do s = 1, NSECTORS
        
            if (ChoiceSim_in(n,coh,year-9) == s) then
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


        do s = 1, NSECTORS
            ExperSim_out(s,n,coh,year) = exper(s)
        end do

                
        do s = 0, NSECTORS
            V(s) = tau(s) + w(s) - cost(lag,s) + eta_aux(s) + rho*Emax(s)
        end do            
                
        VMAX = maxval(V)
                
        do s = 0, NSECTORS   
            if (VMAX == V(s)) then
                    
                ChoiceSim_out(n,coh,year) = s
                WageSim_out(n,coh,year)   = w(s)
                if (lag == 0) then
                    Cost1Sim_out(n,coh,year) = cost(lag,s) + (tau(lag) - tau(s)) +  (eta_aux(lag) - eta_aux(s)) + eps_aux(lag)
                else
                    Cost1Sim_out(n,coh,year) = cost(lag,s) + (tau(lag) - tau(s)) +  (eta_aux(lag) - eta_aux(s))
                end if
                Cost2Sim_out(n,coh,year) = cost(lag,s)
                VSim_out(n,coh,year) = VMAX
                emp_out(s,ed)        = emp_out(s,ed) + 1
                                                
                if (s > 0) then
                        
                    sk_sup_out(s,ed) = sk_sup_out(s,ed) + & 
                    CohortWgt_out(coh,ed)*(w(s)/rsk(s,ed))                              
                        
                end if
                        
            end if
        end do

    end do individuals_loop 

end do cohorts_loop


end subroutine ParallelCohorts_Myopic            
                           
                           
subroutine ParallelCohorts_RE(older_coh, newer_coh, year, param, cut, typ, rsk, rsk_lag, eps, eta, COEF_EXP, covar_EXP, PI_RE, &
                              ChoiceSim_in, ChoiceSim_out, ExperSim_out, WageSim_out, Cost1Sim_out, Cost2Sim_out, VSim_out, CohortWgt_out, emp_out, sk_sup_out)   

USE Global_Data

implicit none

integer, intent(in)             :: older_coh, &
                                   newer_coh, &
                                   year, &
                                   typ(:,FIRST_COH:), &
                                   ChoiceSim_in(:,older_coh:,year-9:)
real(KIND=DOUBLE), intent(in)   :: param(:), &
                                   cut(:,0:), &
                                   rsk(:,0:), &
                                   rsk_lag(:,0:), &
                                   eps(:,older_coh:,year:,0:), &
                                   eta(:,older_coh:,year:,0:), &
                                   COEF_EXP(:,0:,:,:), &
                                   covar_EXP(:,:,0:), &
                                   PI_RE(FIRST_AGE+1:,0:,:,:,:,:)
integer, intent(out)            :: ChoiceSim_out(:,older_coh:,year:), &
                                   ExperSim_out(:,:,older_coh:,year:)
real(KIND=DOUBLE), intent(out)  :: WageSim_out(:,older_coh:,year:), &
                                   Cost1Sim_out(:,older_coh:,year:), &
                                   Cost2Sim_out(:,older_coh:,year:), &
                                   VSim_out(:,older_coh:,year:), &
                                   CohortWgt_out(older_coh:,0:), &
                                   sk_sup_out(:,0:)
integer, intent(out)            :: emp_out(0:,0:)

real(KIND=DOUBLE) theta(7), beta(NSECTORS,12), sigma(0:NSECTORS), kappa(26), &
                  tau(0:NSECTORS), SigmaPref, omega(2:NTYPES,0:NSECTORS), lambda(1:NTYPES), &
                  rsk_init(NSECTORS,0:1)

real(KIND=DOUBLE) cost(0:NSECTORS,0:NSECTORS), &
                  eps_aux(0:NSECTORS), &
                  eta_aux(0:NSECTORS), &
                  new_rsk_tom(NSECTORS,0:1), &
                  sigma2_EXP(NSECTORS,0:1), &
                  alpha_EXP(NSECTORS,0:1), &
                  beta_EXP(NSECTORS,0:1), &
                  Emax(0:NSECTORS), &
                  w(0:NSECTORS), &
                  V(0:NSECTORS), &
                  VMAX

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
rsk_init(1:7,0) = 1.0
rsk_init(1:7,1) = 1.0
       
do ed = 0, 1
    do s = 1, NSECTORS
        alpha_EXP(s,ed)  = COEF_EXP(s,ed,1,1)
        beta_EXP(s,ed)   = COEF_EXP(s,ed,2,1)
        sigma2_EXP(s,ed) = covar_EXP(s,s,ed)
    end do
end do


if (year > FIRST_YR) then
    do s = 1, NSECTORS
        do ed = 0, 1
            new_rsk_tom(s,ed) = rsk(s,ed) * ((rsk(s,ed)/rsk_lag(s,ed))**beta_EXP(s,ed)) &
                                          * exp(alpha_EXP(s,ed) + sigma2_EXP(s,ed)/2)
        end do
    end do
else
    do s = 1, NSECTORS
        do ed = 0, 1
            new_rsk_tom(s,ed) = rsk(s,ed) * ((1.0/rsk_init(s,ed))**beta_EXP(s,ed)) &
                                                 * exp(alpha_EXP(s,ed) + sigma2_EXP(s,ed)/2)
        end do
    end do
end if

dummy = 0
do s = 1, NSECTORS
    do ed = 0, 1
        if (new_rsk_tom(s,ed) >= cut(s,ed)) then
            dummy(s,ed) = 1
        end if
    end do
end do

emp_out = 0
sk_sup_out = 0

cohorts_loop: do coh = older_coh, newer_coh
        
            
    CohortWgt_out(coh,0) =  & 
    ExpansionFactorGlobal*real(CohortSizeData(coh,0)) / real(CohortSizeData(FIRST_COH,0))
    CohortWgt_out(coh,1) =  & 
    ExpansionFactorGlobal*real(CohortSizeData(coh,1)) / real(CohortSizeData(FIRST_COH,0)) 

            
    individuals_loop: do n = 1, NPOP
        
        age     = year - coh 
        gender  = GenderData(n,coh)
        educ    = EducData(n,coh)
        tp      = typ(n,coh)
        lag     = ChoiceSim_in(n,coh,coh+age-1)
        eps_aux = eps(n,coh,year,:)
        eta_aux = eta(n,coh,year,:)

        if (educ <= 2) then
            ed = 0
        else if (educ >= 3) then
            ed = 1
        end if
        
        exper = 0
        do i = 1, 8
            do s = 1, NSECTORS                    
                if (ChoiceSim_in(n,coh,year-i) == s) then
                    exper(s) = exper(s) + 1
                end if                    
            end do
        end do

        if (age < LAST_AGE) then
                           
            do s = 0, NSECTORS    

                ExperTomorrow    = exper
                if (s >= 1) then
                    ExperTomorrow(s) = exper(s) + 1
                end if 
                                
                Emax(s) = &
                PI_RE(age+1,s,gender,educ,tp,1)*1.0 + &
                                
                PI_RE(age+1,s,gender,educ,tp,2)*new_rsk_tom(1,ed) + &
                PI_RE(age+1,s,gender,educ,tp,3)*new_rsk_tom(2,ed) + &
                PI_RE(age+1,s,gender,educ,tp,4)*new_rsk_tom(3,ed) + &
                PI_RE(age+1,s,gender,educ,tp,5)*new_rsk_tom(4,ed) + &
                PI_RE(age+1,s,gender,educ,tp,6)*new_rsk_tom(5,ed) + &
                PI_RE(age+1,s,gender,educ,tp,7)*new_rsk_tom(6,ed) + &
                PI_RE(age+1,s,gender,educ,tp,8)*new_rsk_tom(7,ed) + &
                                
                PI_RE(age+1,s,gender,educ,tp,9)*dummy(1,ed)*(new_rsk_tom(1,ed)-cut(1,ed)) + &
                PI_RE(age+1,s,gender,educ,tp,10)*dummy(2,ed)*(new_rsk_tom(2,ed)-cut(2,ed)) + &
                PI_RE(age+1,s,gender,educ,tp,11)*dummy(3,ed)*(new_rsk_tom(3,ed)-cut(3,ed)) + &
                PI_RE(age+1,s,gender,educ,tp,12)*dummy(4,ed)*(new_rsk_tom(4,ed)-cut(4,ed)) + &
                PI_RE(age+1,s,gender,educ,tp,13)*dummy(5,ed)*(new_rsk_tom(5,ed)-cut(5,ed)) + &
                PI_RE(age+1,s,gender,educ,tp,14)*dummy(6,ed)*(new_rsk_tom(6,ed)-cut(6,ed)) + &
                PI_RE(age+1,s,gender,educ,tp,15)*dummy(7,ed)*(new_rsk_tom(7,ed)-cut(7,ed)) + &
                                
                PI_RE(age+1,s,gender,educ,tp,16)*ExperTomorrow(1) + &
                PI_RE(age+1,s,gender,educ,tp,17)*ExperTomorrow(2) + &
                PI_RE(age+1,s,gender,educ,tp,18)*ExperTomorrow(3) + &
                PI_RE(age+1,s,gender,educ,tp,19)*ExperTomorrow(4) + &
                PI_RE(age+1,s,gender,educ,tp,20)*ExperTomorrow(5) + &
                PI_RE(age+1,s,gender,educ,tp,21)*ExperTomorrow(6) + &
                PI_RE(age+1,s,gender,educ,tp,22)*ExperTomorrow(7) + &
                               
                PI_RE(age+1,s,gender,educ,tp,23)*new_rsk_tom(1,ed)**2 + &
                PI_RE(age+1,s,gender,educ,tp,24)*new_rsk_tom(2,ed)**2 + &
                PI_RE(age+1,s,gender,educ,tp,25)*new_rsk_tom(3,ed)**2 + &
                PI_RE(age+1,s,gender,educ,tp,26)*new_rsk_tom(4,ed)**2 + &
                PI_RE(age+1,s,gender,educ,tp,27)*new_rsk_tom(5,ed)**2 + &
                PI_RE(age+1,s,gender,educ,tp,28)*new_rsk_tom(6,ed)**2 + &
                PI_RE(age+1,s,gender,educ,tp,29)*new_rsk_tom(7,ed)**2 + &
                               
                PI_RE(age+1,s,gender,educ,tp,30)*dummy(1,ed)*(new_rsk_tom(1,ed)-cut(1,ed))**2 + &
                PI_RE(age+1,s,gender,educ,tp,31)*dummy(2,ed)*(new_rsk_tom(2,ed)-cut(2,ed))**2 + &
                PI_RE(age+1,s,gender,educ,tp,32)*dummy(3,ed)*(new_rsk_tom(3,ed)-cut(3,ed))**2 + &
                PI_RE(age+1,s,gender,educ,tp,33)*dummy(4,ed)*(new_rsk_tom(4,ed)-cut(4,ed))**2 + &
                PI_RE(age+1,s,gender,educ,tp,34)*dummy(5,ed)*(new_rsk_tom(5,ed)-cut(5,ed))**2 + &
                PI_RE(age+1,s,gender,educ,tp,35)*dummy(6,ed)*(new_rsk_tom(6,ed)-cut(6,ed))**2 + &
                PI_RE(age+1,s,gender,educ,tp,36)*dummy(7,ed)*(new_rsk_tom(7,ed)-cut(7,ed))**2 + &
                               
                PI_RE(age+1,s,gender,educ,tp,37)*ExperTomorrow(1)**2 + &
                PI_RE(age+1,s,gender,educ,tp,38)*ExperTomorrow(2)**2 + &
                PI_RE(age+1,s,gender,educ,tp,39)*ExperTomorrow(3)**2 + &
                PI_RE(age+1,s,gender,educ,tp,40)*ExperTomorrow(4)**2 + &
                PI_RE(age+1,s,gender,educ,tp,41)*ExperTomorrow(5)**2 + &
                PI_RE(age+1,s,gender,educ,tp,42)*ExperTomorrow(6)**2 + &
                PI_RE(age+1,s,gender,educ,tp,43)*ExperTomorrow(7)**2 + &
                               
                PI_RE(age+1,s,gender,educ,tp,44)*new_rsk_tom(1,ed)*new_rsk_tom(2,ed) + &
                PI_RE(age+1,s,gender,educ,tp,45)*new_rsk_tom(1,ed)*new_rsk_tom(3,ed) + &
                PI_RE(age+1,s,gender,educ,tp,46)*new_rsk_tom(1,ed)*new_rsk_tom(4,ed) + &
                PI_RE(age+1,s,gender,educ,tp,47)*new_rsk_tom(1,ed)*new_rsk_tom(5,ed) + &
                PI_RE(age+1,s,gender,educ,tp,48)*new_rsk_tom(1,ed)*new_rsk_tom(6,ed) + &
                PI_RE(age+1,s,gender,educ,tp,49)*new_rsk_tom(1,ed)*new_rsk_tom(7,ed) + &
                                
                PI_RE(age+1,s,gender,educ,tp,50)*new_rsk_tom(2,ed)*new_rsk_tom(3,ed) + &
                PI_RE(age+1,s,gender,educ,tp,51)*new_rsk_tom(2,ed)*new_rsk_tom(4,ed) + &
                PI_RE(age+1,s,gender,educ,tp,52)*new_rsk_tom(2,ed)*new_rsk_tom(5,ed) + &
                PI_RE(age+1,s,gender,educ,tp,53)*new_rsk_tom(2,ed)*new_rsk_tom(6,ed) + &
                PI_RE(age+1,s,gender,educ,tp,54)*new_rsk_tom(2,ed)*new_rsk_tom(7,ed) + &
                               
                PI_RE(age+1,s,gender,educ,tp,55)*new_rsk_tom(3,ed)*new_rsk_tom(4,ed) + &
                PI_RE(age+1,s,gender,educ,tp,56)*new_rsk_tom(3,ed)*new_rsk_tom(5,ed) + &
                PI_RE(age+1,s,gender,educ,tp,57)*new_rsk_tom(3,ed)*new_rsk_tom(6,ed) + &
                PI_RE(age+1,s,gender,educ,tp,58)*new_rsk_tom(3,ed)*new_rsk_tom(7,ed) + &
                               
                PI_RE(age+1,s,gender,educ,tp,59)*new_rsk_tom(4,ed)*new_rsk_tom(5,ed) + &
                PI_RE(age+1,s,gender,educ,tp,60)*new_rsk_tom(4,ed)*new_rsk_tom(6,ed) + &
                PI_RE(age+1,s,gender,educ,tp,61)*new_rsk_tom(4,ed)*new_rsk_tom(7,ed) + &
                                
                PI_RE(age+1,s,gender,educ,tp,62)*new_rsk_tom(5,ed)*new_rsk_tom(6,ed) + &
                PI_RE(age+1,s,gender,educ,tp,63)*new_rsk_tom(5,ed)*new_rsk_tom(7,ed) + &
                                
                PI_RE(age+1,s,gender,educ,tp,64)*new_rsk_tom(6,ed)*new_rsk_tom(7,ed) + &
                               
                PI_RE(age+1,s,gender,educ,tp,65)*new_rsk_tom(1,ed)*ExperTomorrow(1) + &
                PI_RE(age+1,s,gender,educ,tp,66)*new_rsk_tom(1,ed)*ExperTomorrow(2) + &
                PI_RE(age+1,s,gender,educ,tp,67)*new_rsk_tom(1,ed)*ExperTomorrow(3) + &
                PI_RE(age+1,s,gender,educ,tp,68)*new_rsk_tom(1,ed)*ExperTomorrow(4) + &
                PI_RE(age+1,s,gender,educ,tp,69)*new_rsk_tom(1,ed)*ExperTomorrow(5) + &
                PI_RE(age+1,s,gender,educ,tp,70)*new_rsk_tom(1,ed)*ExperTomorrow(6) + &
                PI_RE(age+1,s,gender,educ,tp,71)*new_rsk_tom(1,ed)*ExperTomorrow(7) + &
                                
                PI_RE(age+1,s,gender,educ,tp,72)*new_rsk_tom(2,ed)*ExperTomorrow(1) + &
                PI_RE(age+1,s,gender,educ,tp,73)*new_rsk_tom(2,ed)*ExperTomorrow(2) + &
                PI_RE(age+1,s,gender,educ,tp,74)*new_rsk_tom(2,ed)*ExperTomorrow(3) + &
                PI_RE(age+1,s,gender,educ,tp,75)*new_rsk_tom(2,ed)*ExperTomorrow(4) + &
                PI_RE(age+1,s,gender,educ,tp,76)*new_rsk_tom(2,ed)*ExperTomorrow(5) + &
                PI_RE(age+1,s,gender,educ,tp,77)*new_rsk_tom(2,ed)*ExperTomorrow(6) + &
                PI_RE(age+1,s,gender,educ,tp,78)*new_rsk_tom(2,ed)*ExperTomorrow(7) + &
                               
                PI_RE(age+1,s,gender,educ,tp,79)*new_rsk_tom(3,ed)*ExperTomorrow(1) + &
                PI_RE(age+1,s,gender,educ,tp,80)*new_rsk_tom(3,ed)*ExperTomorrow(2) + &
                PI_RE(age+1,s,gender,educ,tp,81)*new_rsk_tom(3,ed)*ExperTomorrow(3) + &
                PI_RE(age+1,s,gender,educ,tp,82)*new_rsk_tom(3,ed)*ExperTomorrow(4) + &
                PI_RE(age+1,s,gender,educ,tp,83)*new_rsk_tom(3,ed)*ExperTomorrow(5) + &
                PI_RE(age+1,s,gender,educ,tp,84)*new_rsk_tom(3,ed)*ExperTomorrow(6) + &
                PI_RE(age+1,s,gender,educ,tp,85)*new_rsk_tom(3,ed)*ExperTomorrow(7) + &
                                
                PI_RE(age+1,s,gender,educ,tp,86)*new_rsk_tom(4,ed)*ExperTomorrow(1) + &
                PI_RE(age+1,s,gender,educ,tp,87)*new_rsk_tom(4,ed)*ExperTomorrow(2) + &
                PI_RE(age+1,s,gender,educ,tp,88)*new_rsk_tom(4,ed)*ExperTomorrow(3) + &
                PI_RE(age+1,s,gender,educ,tp,89)*new_rsk_tom(4,ed)*ExperTomorrow(4) + &
                PI_RE(age+1,s,gender,educ,tp,90)*new_rsk_tom(4,ed)*ExperTomorrow(5) + &
                PI_RE(age+1,s,gender,educ,tp,91)*new_rsk_tom(4,ed)*ExperTomorrow(6) + &
                PI_RE(age+1,s,gender,educ,tp,92)*new_rsk_tom(4,ed)*ExperTomorrow(7) + &
                                
                PI_RE(age+1,s,gender,educ,tp,93)*new_rsk_tom(5,ed)*ExperTomorrow(1) + &
                PI_RE(age+1,s,gender,educ,tp,94)*new_rsk_tom(5,ed)*ExperTomorrow(2) + &
                PI_RE(age+1,s,gender,educ,tp,95)*new_rsk_tom(5,ed)*ExperTomorrow(3) + &
                PI_RE(age+1,s,gender,educ,tp,96)*new_rsk_tom(5,ed)*ExperTomorrow(4) + &
                PI_RE(age+1,s,gender,educ,tp,97)*new_rsk_tom(5,ed)*ExperTomorrow(5) + &
                PI_RE(age+1,s,gender,educ,tp,98)*new_rsk_tom(5,ed)*ExperTomorrow(6) + &
                PI_RE(age+1,s,gender,educ,tp,99)*new_rsk_tom(5,ed)*ExperTomorrow(7) + &
                               
                PI_RE(age+1,s,gender,educ,tp,100)*new_rsk_tom(6,ed)*ExperTomorrow(1) + &
                PI_RE(age+1,s,gender,educ,tp,101)*new_rsk_tom(6,ed)*ExperTomorrow(2) + &
                PI_RE(age+1,s,gender,educ,tp,102)*new_rsk_tom(6,ed)*ExperTomorrow(3) + &
                PI_RE(age+1,s,gender,educ,tp,103)*new_rsk_tom(6,ed)*ExperTomorrow(4) + &
                PI_RE(age+1,s,gender,educ,tp,104)*new_rsk_tom(6,ed)*ExperTomorrow(5) + &
                PI_RE(age+1,s,gender,educ,tp,105)*new_rsk_tom(6,ed)*ExperTomorrow(6) + &
                PI_RE(age+1,s,gender,educ,tp,106)*new_rsk_tom(6,ed)*ExperTomorrow(7) + &
                                
                PI_RE(age+1,s,gender,educ,tp,107)*new_rsk_tom(7,ed)*ExperTomorrow(1) + &
                PI_RE(age+1,s,gender,educ,tp,108)*new_rsk_tom(7,ed)*ExperTomorrow(2) + &
                PI_RE(age+1,s,gender,educ,tp,109)*new_rsk_tom(7,ed)*ExperTomorrow(3) + &
                PI_RE(age+1,s,gender,educ,tp,110)*new_rsk_tom(7,ed)*ExperTomorrow(4) + &
                PI_RE(age+1,s,gender,educ,tp,111)*new_rsk_tom(7,ed)*ExperTomorrow(5) + &
                PI_RE(age+1,s,gender,educ,tp,112)*new_rsk_tom(7,ed)*ExperTomorrow(6) + &
                PI_RE(age+1,s,gender,educ,tp,113)*new_rsk_tom(7,ed)*ExperTomorrow(7) + &
                                
                PI_RE(age+1,s,gender,educ,tp,114)*ExperTomorrow(1)*ExperTomorrow(2) + &
                PI_RE(age+1,s,gender,educ,tp,115)*ExperTomorrow(1)*ExperTomorrow(3) + &
                PI_RE(age+1,s,gender,educ,tp,116)*ExperTomorrow(1)*ExperTomorrow(4) + &
                PI_RE(age+1,s,gender,educ,tp,117)*ExperTomorrow(1)*ExperTomorrow(5) + &
                PI_RE(age+1,s,gender,educ,tp,118)*ExperTomorrow(1)*ExperTomorrow(6) + &
                PI_RE(age+1,s,gender,educ,tp,119)*ExperTomorrow(1)*ExperTomorrow(7) + &
                PI_RE(age+1,s,gender,educ,tp,120)*ExperTomorrow(2)*ExperTomorrow(3) + &
                PI_RE(age+1,s,gender,educ,tp,121)*ExperTomorrow(2)*ExperTomorrow(4) + &
                PI_RE(age+1,s,gender,educ,tp,122)*ExperTomorrow(2)*ExperTomorrow(5) + &
                PI_RE(age+1,s,gender,educ,tp,123)*ExperTomorrow(2)*ExperTomorrow(6) + &
                PI_RE(age+1,s,gender,educ,tp,124)*ExperTomorrow(2)*ExperTomorrow(7) + &
                PI_RE(age+1,s,gender,educ,tp,125)*ExperTomorrow(3)*ExperTomorrow(4) + &
                PI_RE(age+1,s,gender,educ,tp,126)*ExperTomorrow(3)*ExperTomorrow(5) + &
                PI_RE(age+1,s,gender,educ,tp,127)*ExperTomorrow(3)*ExperTomorrow(6) + &
                PI_RE(age+1,s,gender,educ,tp,128)*ExperTomorrow(3)*ExperTomorrow(7) + &
                PI_RE(age+1,s,gender,educ,tp,129)*ExperTomorrow(4)*ExperTomorrow(5) + &
                PI_RE(age+1,s,gender,educ,tp,130)*ExperTomorrow(4)*ExperTomorrow(6) + &
                PI_RE(age+1,s,gender,educ,tp,131)*ExperTomorrow(4)*ExperTomorrow(7) + &
                PI_RE(age+1,s,gender,educ,tp,132)*ExperTomorrow(5)*ExperTomorrow(6) + &
                PI_RE(age+1,s,gender,educ,tp,133)*ExperTomorrow(5)*ExperTomorrow(7) + &
                PI_RE(age+1,s,gender,educ,tp,134)*ExperTomorrow(6)*ExperTomorrow(7)

            end do
                       
        else
            Emax = 0.0
        end if


        ! experience accumulated in the last 9 years
        do s = 1, NSECTORS
        
            if (ChoiceSim_in(n,coh,year-9) == s) then
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


        do s = 1, NSECTORS
            ExperSim_out(s,n,coh,year) = exper(s)
        end do

                
        do s = 0, NSECTORS
            V(s) = tau(s) + w(s) - cost(lag,s) + eta_aux(s) + rho*Emax(s)
        end do            
                
        VMAX = maxval(V)
                
        do s = 0, NSECTORS   
            if (VMAX == V(s)) then
                    
                ChoiceSim_out(n,coh,year) = s
                WageSim_out(n,coh,year)   = w(s)
                if (lag == 0) then
                    Cost1Sim_out(n,coh,year) = cost(lag,s) + (tau(lag) - tau(s)) +  (eta_aux(lag) - eta_aux(s)) + eps_aux(lag)
                else
                    Cost1Sim_out(n,coh,year) = cost(lag,s) + (tau(lag) - tau(s)) +  (eta_aux(lag) - eta_aux(s))
                end if
                Cost2Sim_out(n,coh,year) = cost(lag,s)
                VSim_out(n,coh,year) = VMAX
                emp_out(s,ed)        = emp_out(s,ed) + 1
                                                
                if (s > 0) then
                        
                    sk_sup_out(s,ed) = sk_sup_out(s,ed) + & 
                    CohortWgt_out(coh,ed)*(w(s)/rsk(s,ed))                              
                        
                end if
                        
            end if
        end do

    end do individuals_loop 

end do cohorts_loop
                           
end subroutine ParallelCohorts_RE                           
                           
end module ParallelCohorts_MOD                                   
                           
                           
