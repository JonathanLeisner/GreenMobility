MODULE WelfareComp

USE Global_Data

Contains

subroutine ComputeWelfareChanges(Educ,Choice,Exper,CohortWgt,Welfare) 

implicit none

integer          , intent(in) :: Educ(:,FIRST_COH:), &
                                 Choice(:,FIRST_COH:,(FIRST_AGE-9):), &
                                 Exper(:,(2155-40-LAST_AGE):,FIRST_AGE:)
real(KIND=DOUBLE), intent(in) :: CohortWgt(FIRST_COH:,0:), &
                                 Welfare(:,(2155-40-LAST_AGE):,FIRST_AGE:)  

integer, allocatable, dimension(:,:) :: n_index, coh_index      

integer coh, n, s, age, i, t, ed
integer Sizemax, cont(0:NSECTORS), ind(0:NSECTORS)
real(KIND=DOUBLE) TotalWelfare(0:NSECTORS), TotalWgt(0:NSECTORS), &
                  TotalWelfarePre(0:NSECTORS), TotalWelfarePost(0:NSECTORS), &
                  WelfareChangeOverall(0:NSECTORS), WelfareChangeOldNonEd(0:NSECTORS), &
                  WelfareChangeOldEd(0:NSECTORS), WelfareChangeYngNonEd(0:NSECTORS), &
                  WelfareChangeYngEd(0:NSECTORS)

                
                                 
! =======
! Overall
! =======

! =========                                 
! Pre shock
! =========

cont = 0
do coh = HALF_YR_SIM - 40 - LAST_AGE, HALF_YR_SIM - 40 - FIRST_AGE
    age = HALF_YR_SIM - 40 - coh
    do n = 1, NPOP_SIM
        do s = 0, NSECTORS
            if (Choice(n,coh,age) == s) then
                cont(s) = cont(s) + 1
            end if
        end do
    end do   
end do    



Sizemax = maxval(cont)



allocate(coh_index(0:NSECTORS,Sizemax))
allocate(n_index(0:NSECTORS,Sizemax))

ind = 0
do coh = HALF_YR_SIM - 40 - LAST_AGE, HALF_YR_SIM - 40 - FIRST_AGE
    age = HALF_YR_SIM - 40 - coh
    do n = 1, NPOP_SIM
        do s = 0, NSECTORS
            if (Choice(n,coh,age) == s) then
                ind(s)              = ind(s) + 1
                n_index(s,ind(s))   = n
                coh_index(s,ind(s)) = coh
            end if
        end do
    end do
end do                           
                        
TotalWelfare = 0.0
TotalWgt = 0.0

do s = 0, NSECTORS
    do i = 1, ind(s) 
        coh = coh_index(s,i)
        n   = n_index(s,i)
        if (Educ(n,coh) <= 2) then
            ed = 0
        else if (Educ(n,coh) >= 3) then
            ed = 1
        end if
        TotalWgt(s) = TotalWgt(s) + CohortWgt(coh,ed)
        do t = HALF_YR_SIM - 40 + 1, HALF_YR_SIM - 40 + 35
            age = t - coh
            if (age >= FIRST_AGE .and. age <= LAST_AGE) then
                TotalWelfare(s) = TotalWelfare(s) + &
                (rho**(t - (HALF_YR_SIM - 40 + 1)))*Welfare(n,coh,age)*CohortWgt(coh,ed)
            end if
        end do
    end do  
end do    

TotalWelfarePre = TotalWelfare / TotalWgt
! TotalWelfarePre = TotalWelfare / ind

deallocate(coh_index, n_index)


! ==========
! Post shock
! ==========
cont = 0
do coh = HALF_YR_SIM - LAST_AGE, HALF_YR_SIM - FIRST_AGE
    age = HALF_YR_SIM - coh
    do n = 1, NPOP_SIM
        do s = 0, NSECTORS
            if (Choice(n,coh,age) == s) then
                cont(s) = cont(s) + 1
            end if
        end do
    end do   
end do    



Sizemax = maxval(cont)



allocate(coh_index(0:NSECTORS,Sizemax))
allocate(n_index(0:NSECTORS,Sizemax))

ind = 0
coh_index = 0
n_index = 0
do coh = HALF_YR_SIM - LAST_AGE, HALF_YR_SIM - FIRST_AGE
    age = HALF_YR_SIM - coh
    do n = 1, NPOP_SIM
        do s = 0, NSECTORS
            if (Choice(n,coh,age) == s) then
                ind(s)              = ind(s) + 1
                coh_index(s,ind(s)) = coh
                n_index(s,ind(s))   = n
            end if
        end do
    end do
end do                           
                        
TotalWelfare = 0.0
TotalWgt = 0.0

do s = 0, NSECTORS
    do i = 1, ind(s)
        n   = n_index(s,i)
        coh = coh_index(s,i)
        if (Educ(n,coh) <= 2) then
            ed = 0
        else if (Educ(n,coh) >= 3) then
            ed = 1
        end if
        TotalWgt(s) = TotalWgt(s) + CohortWgt(coh,ed)
        do t = HALF_YR_SIM + 1, HALF_YR_SIM + 35
            age = t - coh
            if (age >= FIRST_AGE .and. age <= LAST_AGE) then
                TotalWelfare(s) = TotalWelfare(s) + &
                (0.95**(t - (HALF_YR_SIM + 1)))*Welfare(n,coh,age)*CohortWgt(coh,ed)
            end if
        end do
    end do
end do

TotalWelfarePost = TotalWelfare / TotalWgt
!TotalWelfarePost = TotalWelfare / ind

deallocate(coh_index, n_index)

WelfareChangeOverall = (TotalWelfarePost - TotalWelfarePre)/TotalWelfarePre


! **************************************************************************************************
! **************************************************************************************************

! ====================
! Old and Non Educated                        
! ====================

! =========                                 
! Pre shock
! =========

cont = 0
do coh = HALF_YR_SIM - 40 - LAST_AGE, HALF_YR_SIM - 40 - FIRST_AGE
    age = HALF_YR_SIM - 40 - coh
    do n = 1, NPOP_SIM
        do s = 0, NSECTORS
            if (Educ(n,coh) <= 2 .and. coh >= HALF_YR_SIM - 40 - 59 .and. &
                coh <= HALF_YR_SIM - 40 - 45) then
                if (Choice(n,coh,age) == s) then
                    cont(s) = cont(s) + 1
                end if
            end if
        end do
    end do   
end do 

Sizemax = maxval(cont)



allocate(coh_index(0:NSECTORS,Sizemax))
allocate(n_index(0:NSECTORS,Sizemax))

ind = 0
do coh = HALF_YR_SIM - 40 - LAST_AGE, HALF_YR_SIM - 40 - FIRST_AGE
    age = HALF_YR_SIM - 40 - coh
    do n = 1, NPOP_SIM
        do s = 0, NSECTORS
            if (Educ(n,coh) <= 2 .and. coh >= HALF_YR_SIM - 40 - 59 .and. &
                coh <= HALF_YR_SIM - 40 - 45) then
                if (Choice(n,coh,age) == s) then
                    ind(s)              = ind(s) + 1
                    n_index(s,ind(s))   = n
                    coh_index(s,ind(s)) = coh
                end if
            end if
        end do
    end do
end do                           
                        
TotalWelfare = 0.0
TotalWgt = 0.0

do s = 0, NSECTORS
    do i = 1, ind(s) 
        coh = coh_index(s,i)
        n   = n_index(s,i)
        if (Educ(n,coh) <= 2) then
            ed = 0
        else if (Educ(n,coh) >= 3) then
            ed = 1
        end if
        TotalWgt(s) = TotalWgt(s) + CohortWgt(coh,ed)
        do t = HALF_YR_SIM - 40 + 1, HALF_YR_SIM - 40 + 35
            age = t - coh
            if (age >= FIRST_AGE .and. age <= LAST_AGE) then
                TotalWelfare(s) = TotalWelfare(s) + &
                (rho**(t - (HALF_YR_SIM - 40 + 1)))*Welfare(n,coh,age)*CohortWgt(coh,ed)
            end if
        end do
    end do  
end do    

TotalWelfarePre = TotalWelfare / TotalWgt
!TotalWelfarePre = TotalWelfare / ind

deallocate(coh_index, n_index)


! ==========
! Post shock
! ==========

cont = 0
do coh = HALF_YR_SIM - LAST_AGE, HALF_YR_SIM - FIRST_AGE
    age = HALF_YR_SIM - coh
    do n = 1, NPOP_SIM
        do s = 0, NSECTORS
            if (Educ(n,coh) <= 2 .and. coh >= HALF_YR_SIM - 59 .and. &
                coh <= HALF_YR_SIM - 45) then
                if (Choice(n,coh,age) == s) then
                    cont(s) = cont(s) + 1
                end if
            end if
        end do
    end do   
end do    



Sizemax = maxval(cont)



allocate(coh_index(0:NSECTORS,Sizemax))
allocate(n_index(0:NSECTORS,Sizemax))

ind = 0
coh_index = 0
n_index = 0
do coh = HALF_YR_SIM - LAST_AGE, HALF_YR_SIM - FIRST_AGE
    age = HALF_YR_SIM - coh
    do n = 1, NPOP_SIM
        do s = 0, NSECTORS
            if (Educ(n,coh) <= 2 .and. coh >= HALF_YR_SIM - 59 .and. &
                coh <= HALF_YR_SIM - 45) then
                if (Choice(n,coh,age) == s) then
                    ind(s)              = ind(s) + 1
                    coh_index(s,ind(s)) = coh
                    n_index(s,ind(s))   = n
                end if
            end if
        end do
    end do
end do                           
                        
TotalWelfare = 0.0
TotalWgt = 0.0

do s = 0, NSECTORS
    do i = 1, ind(s)
        n   = n_index(s,i)
        coh = coh_index(s,i)
        if (Educ(n,coh) <= 2) then
            ed = 0
        else if (Educ(n,coh) >= 3) then
            ed = 1
        end if
        TotalWgt(s) = TotalWgt(s) + CohortWgt(coh,ed)
        do t = HALF_YR_SIM + 1, HALF_YR_SIM + 35
            age = t - coh
            if (age >= FIRST_AGE .and. age <= LAST_AGE) then
                TotalWelfare(s) = TotalWelfare(s) + &
                (0.95**(t - (HALF_YR_SIM + 1)))*Welfare(n,coh,age)*CohortWgt(coh,ed)
            end if
        end do
    end do
end do

TotalWelfarePost = TotalWelfare / TotalWgt
!TotalWelfarePost = TotalWelfare / ind

deallocate(coh_index, n_index)

WelfareChangeOldNonEd = (TotalWelfarePost - TotalWelfarePre)/TotalWelfarePre


! **************************************************************************************************
! **************************************************************************************************

! ====================
! Old and Educated                        
! ====================

! =========                                 
! Pre shock
! =========

cont = 0
do coh = HALF_YR_SIM - 40 - LAST_AGE, HALF_YR_SIM - 40 - FIRST_AGE
    age = HALF_YR_SIM - 40 - coh
    do n = 1, NPOP_SIM
        do s = 0, NSECTORS
            if (Educ(n,coh) >= 3 .and. coh >= HALF_YR_SIM - 40 - 59 .and. &
                coh <= HALF_YR_SIM - 40 - 45) then
                if (Choice(n,coh,age) == s) then
                    cont(s) = cont(s) + 1
                end if
            end if
        end do
    end do   
end do 

Sizemax = maxval(cont)



allocate(coh_index(0:NSECTORS,Sizemax))
allocate(n_index(0:NSECTORS,Sizemax))

ind = 0
do coh = HALF_YR_SIM - 40 - LAST_AGE, HALF_YR_SIM - 40 - FIRST_AGE
    age = HALF_YR_SIM - 40 - coh
    do n = 1, NPOP_SIM
        do s = 0, NSECTORS
            if (Educ(n,coh) >= 3 .and. coh >= HALF_YR_SIM - 40 - 59 .and. &
                coh <= HALF_YR_SIM - 40 - 45) then
                if (Choice(n,coh,age) == s) then
                    ind(s)              = ind(s) + 1
                    n_index(s,ind(s))   = n
                    coh_index(s,ind(s)) = coh
                end if
            end if
        end do
    end do
end do                           
                        
TotalWelfare = 0.0
TotalWgt = 0.0

do s = 0, NSECTORS
    do i = 1, ind(s) 
        coh = coh_index(s,i)
        n   = n_index(s,i)
        if (Educ(n,coh) <= 2) then
            ed = 0
        else if (Educ(n,coh) >= 3) then
            ed = 1
        end if
        TotalWgt(s) = TotalWgt(s) + CohortWgt(coh,ed)
        do t = HALF_YR_SIM - 40 + 1, HALF_YR_SIM - 40 + 35
            age = t - coh
            if (age >= FIRST_AGE .and. age <= LAST_AGE) then
                TotalWelfare(s) = TotalWelfare(s) + &
                (rho**(t - (HALF_YR_SIM - 40 + 1)))*Welfare(n,coh,age)*CohortWgt(coh,ed)
            end if
        end do
    end do  
end do    

TotalWelfarePre = TotalWelfare / TotalWgt
!TotalWelfarePre = TotalWelfare / ind

deallocate(coh_index, n_index)


! ==========
! Post shock
! ==========

cont = 0
do coh = HALF_YR_SIM - LAST_AGE, HALF_YR_SIM - FIRST_AGE
    age = HALF_YR_SIM - coh
    do n = 1, NPOP_SIM
        do s = 0, NSECTORS
            if (Educ(n,coh) >= 3 .and. coh >= HALF_YR_SIM - 59 .and. &
                coh <= HALF_YR_SIM - 45) then
                if (Choice(n,coh,age) == s) then
                    cont(s) = cont(s) + 1
                end if
            end if
        end do
    end do   
end do    



Sizemax = maxval(cont)



allocate(coh_index(0:NSECTORS,Sizemax))
allocate(n_index(0:NSECTORS,Sizemax))

ind = 0
coh_index = 0
n_index = 0
do coh = HALF_YR_SIM - LAST_AGE, HALF_YR_SIM - FIRST_AGE
    age = HALF_YR_SIM - coh
    do n = 1, NPOP_SIM
        do s = 0, NSECTORS
            if (Educ(n,coh) >= 3 .and. coh >= HALF_YR_SIM - 59 .and. &
                coh <= HALF_YR_SIM - 45) then
                if (Choice(n,coh,age) == s) then
                    ind(s)              = ind(s) + 1
                    coh_index(s,ind(s)) = coh
                    n_index(s,ind(s))   = n
                end if
            end if
        end do
    end do
end do                           
                        
TotalWelfare = 0.0
TotalWgt = 0.0

do s = 0, NSECTORS
    do i = 1, ind(s)
        n   = n_index(s,i)
        coh = coh_index(s,i)
        if (Educ(n,coh) <= 2) then
            ed = 0
        else if (Educ(n,coh) >= 3) then
            ed = 1
        end if
        TotalWgt(s) = TotalWgt(s) + CohortWgt(coh,ed)
        do t = HALF_YR_SIM + 1, HALF_YR_SIM + 35
            age = t - coh
            if (age >= FIRST_AGE .and. age <= LAST_AGE) then
                TotalWelfare(s) = TotalWelfare(s) + &
                (0.95**(t - (HALF_YR_SIM + 1)))*Welfare(n,coh,age)*CohortWgt(coh,ed)
            end if
        end do
    end do
end do

TotalWelfarePost = TotalWelfare / TotalWgt
!TotalWelfarePost = TotalWelfare / ind

deallocate(coh_index, n_index)

WelfareChangeOldEd = (TotalWelfarePost - TotalWelfarePre)/TotalWelfarePre


! **************************************************************************************************
! **************************************************************************************************

! ======================
! Young and Non Educated                        
! ======================

! =========                                 
! Pre shock
! =========

cont = 0
do coh = HALF_YR_SIM - 40 - LAST_AGE, HALF_YR_SIM - 40 - FIRST_AGE
    age = HALF_YR_SIM - 40 - coh
    do n = 1, NPOP_SIM
        do s = 0, NSECTORS
            if (Educ(n,coh) <= 2 .and. coh >= HALF_YR_SIM - 40 - 39 .and. &
                coh <= HALF_YR_SIM - 40 - 25) then
                if (Choice(n,coh,age) == s) then
                    cont(s) = cont(s) + 1
                end if
            end if
        end do
    end do   
end do 

Sizemax = maxval(cont)



allocate(coh_index(0:NSECTORS,Sizemax))
allocate(n_index(0:NSECTORS,Sizemax))

ind = 0
do coh = HALF_YR_SIM - 40 - LAST_AGE, HALF_YR_SIM - 40 - FIRST_AGE
    age = HALF_YR_SIM - 40 - coh
    do n = 1, NPOP_SIM
        do s = 0, NSECTORS
            if (Educ(n,coh) <= 2 .and. coh >= HALF_YR_SIM - 40 - 39 .and. &
                coh <= HALF_YR_SIM - 40 - 25) then
                if (Choice(n,coh,age) == s) then
                    ind(s)              = ind(s) + 1
                    n_index(s,ind(s))   = n
                    coh_index(s,ind(s)) = coh
                end if
            end if
        end do
    end do
end do                           
                        
TotalWelfare = 0.0
TotalWgt = 0.0

do s = 0, NSECTORS
    do i = 1, ind(s) 
        coh = coh_index(s,i)
        n   = n_index(s,i)
        if (Educ(n,coh) <= 2) then
            ed = 0
        else if (Educ(n,coh) >= 3) then
            ed = 1
        end if
        TotalWgt(s) = TotalWgt(s) + CohortWgt(coh,ed)
        do t = HALF_YR_SIM - 40 + 1, HALF_YR_SIM - 40 + 35
            age = t - coh
            if (age >= FIRST_AGE .and. age <= LAST_AGE) then
                TotalWelfare(s) = TotalWelfare(s) + &
                (rho**(t - (HALF_YR_SIM - 40 + 1)))*Welfare(n,coh,age)*CohortWgt(coh,ed)
            end if
        end do
    end do  
end do    

TotalWelfarePre = TotalWelfare / TotalWgt
!TotalWelfarePre = TotalWelfare / ind

deallocate(coh_index, n_index)


! ==========
! Post shock
! ==========

cont = 0
do coh = HALF_YR_SIM - LAST_AGE, HALF_YR_SIM - FIRST_AGE
    age = HALF_YR_SIM - coh
    do n = 1, NPOP_SIM
        do s = 0, NSECTORS
            if (Educ(n,coh) <= 2 .and. coh >= HALF_YR_SIM - 39 .and. &
                coh <= HALF_YR_SIM - 25) then
                if (Choice(n,coh,age) == s) then
                    cont(s) = cont(s) + 1
                end if
            end if
        end do
    end do   
end do    



Sizemax = maxval(cont)



allocate(coh_index(0:NSECTORS,Sizemax))
allocate(n_index(0:NSECTORS,Sizemax))

ind = 0
coh_index = 0
n_index = 0
do coh = HALF_YR_SIM - LAST_AGE, HALF_YR_SIM - FIRST_AGE
    age = HALF_YR_SIM - coh
    do n = 1, NPOP_SIM
        do s = 0, NSECTORS
            if (Educ(n,coh) <= 2 .and. coh >= HALF_YR_SIM - 39 .and. &
                coh <= HALF_YR_SIM - 25) then
                if (Choice(n,coh,age) == s) then
                    ind(s)              = ind(s) + 1
                    coh_index(s,ind(s)) = coh
                    n_index(s,ind(s))   = n
                end if
            end if
        end do
    end do
end do                           
                        
TotalWelfare = 0.0
TotalWgt = 0.0

do s = 0, NSECTORS
    do i = 1, ind(s)
        n   = n_index(s,i)
        coh = coh_index(s,i)
        if (Educ(n,coh) <= 2) then
            ed = 0
        else if (Educ(n,coh) >= 3) then
            ed = 1
        end if
        TotalWgt(s) = TotalWgt(s) + CohortWgt(coh,ed)
        do t = HALF_YR_SIM + 1, HALF_YR_SIM + 35
            age = t - coh
            if (age >= FIRST_AGE .and. age <= LAST_AGE) then
                TotalWelfare(s) = TotalWelfare(s) + &
                (0.95**(t - (HALF_YR_SIM + 1)))*Welfare(n,coh,age)*CohortWgt(coh,ed)
            end if
        end do
    end do
end do

TotalWelfarePost = TotalWelfare / TotalWgt
!TotalWelfarePost = TotalWelfare / ind

deallocate(coh_index, n_index)

WelfareChangeYngNonEd = (TotalWelfarePost - TotalWelfarePre)/TotalWelfarePre



! **************************************************************************************************
! **************************************************************************************************

! ======================
! Young and Educated                        
! ======================

! =========                                 
! Pre shock
! =========

cont = 0
do coh = HALF_YR_SIM - 40 - LAST_AGE, HALF_YR_SIM - 40 - FIRST_AGE
    age = HALF_YR_SIM - 40 - coh
    do n = 1, NPOP_SIM
        do s = 0, NSECTORS
            if (Educ(n,coh) >= 3 .and. coh >= HALF_YR_SIM - 40 - 39 .and. &
                coh <= HALF_YR_SIM - 40 - 25) then
                if (Choice(n,coh,age) == s) then
                    cont(s) = cont(s) + 1
                end if
            end if
        end do
    end do   
end do 

Sizemax = maxval(cont)



allocate(coh_index(0:NSECTORS,Sizemax))
allocate(n_index(0:NSECTORS,Sizemax))

ind = 0
do coh = HALF_YR_SIM - 40 - LAST_AGE, HALF_YR_SIM - 40 - FIRST_AGE
    age = HALF_YR_SIM - 40 - coh
    do n = 1, NPOP_SIM
        do s = 0, NSECTORS
            if (Educ(n,coh) >= 3 .and. coh >= HALF_YR_SIM - 40 - 39 .and. &
                coh <= HALF_YR_SIM - 40 - 25) then
                if (Choice(n,coh,age) == s) then
                    ind(s)              = ind(s) + 1
                    n_index(s,ind(s))   = n
                    coh_index(s,ind(s)) = coh
                end if
            end if
        end do
    end do
end do                           
                        
TotalWelfare = 0.0
TotalWgt = 0.0

do s = 0, NSECTORS
    do i = 1, ind(s) 
        coh = coh_index(s,i)
        n   = n_index(s,i)
        if (Educ(n,coh) <= 2) then
            ed = 0
        else if (Educ(n,coh) >= 3) then
            ed = 1
        end if
        TotalWgt(s) = TotalWgt(s) + CohortWgt(coh,ed)
        do t = HALF_YR_SIM - 40 + 1, HALF_YR_SIM - 40 + 35
            age = t - coh
            if (age >= FIRST_AGE .and. age <= LAST_AGE) then
                TotalWelfare(s) = TotalWelfare(s) + &
                (rho**(t - (HALF_YR_SIM - 40 + 1)))*Welfare(n,coh,age)*CohortWgt(coh,ed)
            end if
        end do
    end do  
end do    

TotalWelfarePre = TotalWelfare / TotalWgt
!TotalWelfarePre = TotalWelfare / ind

deallocate(coh_index, n_index)


! ==========
! Post shock
! ==========

cont = 0
do coh = HALF_YR_SIM - LAST_AGE, HALF_YR_SIM - FIRST_AGE
    age = HALF_YR_SIM - coh
    do n = 1, NPOP_SIM
        do s = 0, NSECTORS
            if (Educ(n,coh) >= 3 .and. coh >= HALF_YR_SIM - 39 .and. &
                coh <= HALF_YR_SIM - 25) then
                if (Choice(n,coh,age) == s) then
                    cont(s) = cont(s) + 1
                end if
            end if
        end do
    end do   
end do    



Sizemax = maxval(cont)



allocate(coh_index(0:NSECTORS,Sizemax))
allocate(n_index(0:NSECTORS,Sizemax))

ind = 0
coh_index = 0
n_index = 0
do coh = HALF_YR_SIM - LAST_AGE, HALF_YR_SIM - FIRST_AGE
    age = HALF_YR_SIM - coh
    do n = 1, NPOP_SIM
        do s = 0, NSECTORS
            if (Educ(n,coh) >= 3 .and. coh >= HALF_YR_SIM - 39 .and. &
                coh <= HALF_YR_SIM - 25) then
                if (Choice(n,coh,age) == s) then
                    ind(s)              = ind(s) + 1
                    coh_index(s,ind(s)) = coh
                    n_index(s,ind(s))   = n
                end if
            end if
        end do
    end do
end do                           
                        
TotalWelfare = 0.0
TotalWgt = 0.0

do s = 0, NSECTORS
    do i = 1, ind(s)
        n   = n_index(s,i)
        coh = coh_index(s,i)
        if (Educ(n,coh) <= 2) then
            ed = 0
        else if (Educ(n,coh) >= 3) then
            ed = 1
        end if
        TotalWgt(s) = TotalWgt(s) + CohortWgt(coh,ed)
        do t = HALF_YR_SIM + 1, HALF_YR_SIM + 35
            age = t - coh
            if (age >= FIRST_AGE .and. age <= LAST_AGE) then
                TotalWelfare(s) = TotalWelfare(s) + &
                (0.95**(t - (HALF_YR_SIM + 1)))*Welfare(n,coh,age)*CohortWgt(coh,ed)
            end if
        end do
    end do
end do

TotalWelfarePost = TotalWelfare / TotalWgt
!TotalWelfarePost = TotalWelfare / ind

deallocate(coh_index, n_index)

WelfareChangeYngEd = (TotalWelfarePost - TotalWelfarePre)/TotalWelfarePre


if (CAPMOBILITY == 1) then
    open(unit = 1, file='WelfareChanges1.csv')
else if (CAPMOBILITY == 2) then
    open(unit = 1, file='WelfareChanges2.csv')
else if (CAPMOBILITY == 3) then
    open(unit = 1, file='WelfareChanges3.csv')
end if

    write(1,1234) ' , Residual, Agriculture/Mining, LT Manufacturing, HT Manufacturing, Construction, Trade, Trans/Util, Service'
    write(1,4321) 'Overall'        , (100*WelfareChangeOverall(s) , s = 0, NSECTORS)
    write(1,4321) 'Old/Unskilled'  , (100*WelfareChangeOldNonEd(s), s = 0, NSECTORS)
    write(1,4321) 'Old/Skilled'    , (100*WelfareChangeOldEd(s)   , s = 0, NSECTORS)
    write(1,4321) 'Young/Unskilled', (100*WelfareChangeYngNonEd(s), s = 0, NSECTORS)
    write(1,4321) 'Young/Skilled'  , (100*WelfareChangeYngEd(s)   , s = 0, NSECTORS)

close(1)

1234 format(a120)
4321 format(a20,',',8(f12.6,','))

                                 
end subroutine ComputeWelfareChanges

end MODULE WelfareComp                                     



