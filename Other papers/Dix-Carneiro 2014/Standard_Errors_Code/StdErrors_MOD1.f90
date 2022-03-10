Module StdErrors_MOD

Contains

Subroutine NumDer(param_in)

! **********************************************************
! This subroutine computes the numerical derivatives (G0) at
! point param_in
! **********************************************************

USE Global_Data
USE Loss_Function_MOD
USE LinReg_MOD
USE MPI

implicit none

 
real(KIND=DOUBLE) , intent(in)  :: param_in(:)

integer, parameter :: NGRID = 10

! Model parameters
real(KIND=DOUBLE) theta(7), beta(NSECTORS,12), & 
                  sigma(0:NSECTORS), kappa(26), tau(0:NSECTORS), SigmaPref, &
                  omega(2:NTYPES,0:NSECTORS), lambda(NTYPES), gamma(2:NTYPES,0:NSECTORS+4), aa(NSECTORS,0:1,0:1), &
                  sigma_prod(NSECTORS), param(NPARAM)

real(KIND=DOUBLE) low_theta(7), up_theta(7), &
                  low_beta(NSECTORS,12), up_beta(NSECTORS,12), &
                  low_sigma(0:NSECTORS), up_sigma(0:NSECTORS), &
                  low_kappa(26), up_kappa(26), &
                  low_tau(0:NSECTORS), up_tau(0:NSECTORS), &     
                  low_SigmaPref, up_SigmaPref, &
                  low_omega(2:NTYPES,0:NSECTORS), up_omega(2:NTYPES,0:NSECTORS), &
                  low_lambda(NTYPES), up_lambda(NTYPES), &
                  low_gamma(2:NTYPES,0:NSECTORS+4), up_gamma(2:NTYPES,0:NSECTORS+4), &
                  low_aa(NSECTORS,0:1,0:1), up_aa(NSECTORS,0:1,0:1), &
                  low_sigma_prod(NSECTORS), up_sigma_prod(NSECTORS), &
                  F, X(NPARAM,2*NGRID+1), &
                  delta(NPARAM,2*NGRID+1,NREG_Wage), &
                  XR(2*NGRID+1,3), Y(2*NGRID+1,1), &
                  COEF(3,1), SST, SSE

real(KIND=DOUBLE), allocatable, dimension(:,:) :: G0

integer i, j, n, s, s1, s2, eq, npar, status, info, tp, ed, start, dim_G0

character char, char1, char2

character*25 filename

! MPI Variables
integer ierr, numproc, rank, tag, source, count1, count2, proc, ierror
integer status2(MPI_STATUS_SIZE)

Call MPI_INIT(ierr)
Call MPI_COMM_SIZE (MPI_COMM_WORLD,numproc,ierr)
Call MPI_COMM_RANK (MPI_COMM_WORLD,rank,ierr)


! *****************************
! Use only the Master processor
! *****************************

if (rank == 0) then 

    dim_G0 = NSECTORS*NREG_Wage + & 
            (1+NSECTORS)*NREG_Emp + &
            (1+NSECTORS)*(1+NSECTORS)*NREG_Tr + &
            NSECTORS + &
            NSECTORS + &
            (1+NSECTORS)*NREG_Pers*3 + &
            (1+NSECTORS)*NREG_Freq + &
            (1+NSECTORS)*NREG_Return  

    theta(1:7)      = param_in(1:7)
    beta(1,1:12)    = param_in(8:19)
    beta(2,1:12)    = param_in(20:31)
    beta(3,1:12)    = param_in(32:43)
    beta(4,1:12)    = param_in(44:55)
    beta(5,1:12)    = param_in(56:67)
    beta(6,1:12)    = param_in(68:79)
    beta(7,1:12)    = param_in(80:91)
    sigma(0:7)      = param_in(92:99)
    kappa(1:26)     = param_in(100:125)
    tau             = 0.0
    tau(2:7)        = param_in(126:131)
    SigmaPref       = param_in(132)
    omega(2,0:7)    = param_in(133:140)
    omega(3,0:7)    = param_in(141:148)
    lambda          = 0.0
    lambda(2:3)     = param_in(149:150)
    gamma(2,0:11)   = param_in(151:162)
    gamma(3,0:11)   = param_in(163:174)
    aa(1,0,0:1)     = param_in(175:176)
    aa(2,0,0:1)     = param_in(177:178)
    aa(3,0,0:1)     = param_in(179:180)
    aa(4,0,0:1)     = param_in(181:182)
    aa(5,0,0:1)     = param_in(183:184)
    aa(6,0,0:1)     = param_in(185:186)
    aa(7,0,0:1)     = param_in(187:188)
    aa(1,1,0:1)     = param_in(189:190)
    aa(2,1,0:1)     = param_in(191:192)
    aa(3,1,0:1)     = param_in(193:194)
    aa(4,1,0:1)     = param_in(195:196)
    aa(5,1,0:1)     = param_in(197:198)
    aa(6,1,0:1)     = param_in(199:200)
    aa(7,1,0:1)     = param_in(201:202)
    sigma_prod(1:7) = param_in(203:209)

    ! Read range of parameters along which the interpolation will be performed

    open(unit = 10, file = 'StdErrors_Range.csv')

        do i = 1, 7
            read(10,*) n, low_theta(i), up_theta(i)
        end do

        do s = 1, NSECTORS
            do i = 1, 12
                read(10,*) n, low_beta(s,i), up_beta(s,i)
            end do
        end do    

        do s = 0, NSECTORS
            read(10,*) n, low_sigma(s), up_sigma(s)
        end do

        do i = 1, 26
            read(10,*) n, low_kappa(i), up_kappa(i)
        end do  
    
        do i = 2, 7
            read(10,*) n, low_tau(i), up_tau(i)
        end do
    
        read(10,*) n, low_SigmaPref, up_SigmaPref
        
        do tp = 2, NTYPES
            do i = 0, NSECTORS
                read(10,*) n, low_omega(tp,i), up_omega(tp,i)
            end do
        end do
        
        do i = 2, NTYPES
            read(10,*) n, low_lambda(i), up_lambda(i)
        end do
    
        do i = 2, NTYPES
            do s = 0, NSECTORS+4
                read(10,*) n, low_gamma(i,s), up_gamma(i,s)
            end do
        end do
        
        do ed = 0, 1
            do s = 1, NSECTORS
                do i = 0, 1
                    read(10,*) n, low_aa(s,ed,i), up_aa(s,ed,i)
                end do
            end do
        end do
        
        do s = 1, NSECTORS
            read(10,*), n, low_sigma_prod(s), up_sigma_prod(s)
        end do

    close(1)


    ! Check whether the ranges are well defined

    do i = 1, 7
        if(low_theta(i) > theta(i) .or. up_theta(i) < theta(i)) then
            print*, "Error: theta", i
            pause
            stop
        end if
    end do

    do s = 1, NSECTORS
        do j = 1, 12
            if(low_beta(s,j) > beta(s,j) .or. up_beta(s,j) < beta(s,j)) then
                print*, "Error: beta", s, j
                pause
                stop
            end if
        end do
    end do

    do s = 0, NSECTORS
        if(low_sigma(s) > sigma(s) .or. up_sigma(s) < sigma(s)) then
            print*, "Error: sigma", s
            pause
            stop
        end if
    end do

    do i = 1, 26
        if(low_kappa(i) > kappa(i) .or. up_kappa(i) < kappa(i)) then
            print*, "Error: kappa", i
            pause
            stop
        end if
    end do

    do s = 2, NSECTORS
        if(low_tau(s) > tau(s) .or. up_tau(s) < tau(s)) then
            print*, "Error: alpha", s
            pause
            stop
        end if
    end do

    if(low_SigmaPref > SigmaPref .or. up_SigmaPref < SigmaPref) then
        print*, "Error: SigmaPref"
        pause
        stop
    end if

    do tp = 2, NTYPES
        do s = 0, NSECTORS
            if(low_omega(tp,s) > omega(tp,s) .or. up_omega(tp,s) < omega(tp,s)) then
                print*, "Error: omega", tp, s
                pause
                stop
            end if
        end do        
    end do

    do tp = 2, NTYPES
        do s = 0, NSECTORS+4
            if(low_gamma(tp,s) > gamma(tp,s) .or. up_gamma(tp,s) < gamma(tp,s)) then
                print*, "Error: gamma", tp, s
                pause
                stop
            end if
        end do        
    end do
    
    do tp = 2, NTYPES
        if(low_lambda(tp) > lambda(tp) .or. up_lambda(tp) < lambda(tp)) then
            print*, "Error: lambda", tp
            pause
            stop
        end if
    end do
    
    do s = 1, NSECTORS
        do ed = 0, 1
            do i = 0, 1
                if(low_aa(s,ed,i) > aa(s,ed,i) .or. up_aa(s,ed,i) < aa(s,ed,i)) then
                    print*, "Error: aa", s, ed, i
                    pause
                    stop
                end if
            end do
        end do        
    end do
    
    do s = 1, NSECTORS
        if(low_sigma_prod(s) > sigma_prod(s) .or. up_sigma_prod(s) < sigma_prod(s)) then
            print*, "Error: sigma_prod", s
                pause
                stop
        end if        
    end do

    open(unit = 101, file = 'stderrors1/beta1.csv')
    open(unit = 102, file = 'stderrors1/beta2.csv')
    open(unit = 103, file = 'stderrors1/beta3.csv')
    open(unit = 104, file = 'stderrors1/beta4.csv')
    open(unit = 105, file = 'stderrors1/beta5.csv')
    open(unit = 106, file = 'stderrors1/beta6.csv')
    open(unit = 107, file = 'stderrors1/beta7.csv')
    
    open(unit = 200, file = 'stderrors1/gamma0.csv')
    open(unit = 201, file = 'stderrors1/gamma1.csv')
    open(unit = 202, file = 'stderrors1/gamma2.csv')
    open(unit = 203, file = 'stderrors1/gamma3.csv')
    open(unit = 204, file = 'stderrors1/gamma4.csv')
    open(unit = 205, file = 'stderrors1/gamma5.csv')
    open(unit = 206, file = 'stderrors1/gamma6.csv')
    open(unit = 207, file = 'stderrors1/gamma7.csv')
    
    open(unit = 300, file = 'stderrors1/phi00.csv')
    open(unit = 301, file = 'stderrors1/phi01.csv')
    open(unit = 302, file = 'stderrors1/phi02.csv')
    open(unit = 303, file = 'stderrors1/phi03.csv')
    open(unit = 304, file = 'stderrors1/phi04.csv')
    open(unit = 305, file = 'stderrors1/phi05.csv')
    open(unit = 306, file = 'stderrors1/phi06.csv')
    open(unit = 307, file = 'stderrors1/phi07.csv')
    
    open(unit = 310, file = 'stderrors1/phi10.csv')
    open(unit = 311, file = 'stderrors1/phi11.csv')
    open(unit = 312, file = 'stderrors1/phi12.csv')
    open(unit = 313, file = 'stderrors1/phi13.csv')
    open(unit = 314, file = 'stderrors1/phi14.csv')
    open(unit = 315, file = 'stderrors1/phi15.csv')
    open(unit = 316, file = 'stderrors1/phi16.csv')
    open(unit = 317, file = 'stderrors1/phi17.csv')
    
    open(unit = 320, file = 'stderrors1/phi20.csv')
    open(unit = 321, file = 'stderrors1/phi21.csv')
    open(unit = 322, file = 'stderrors1/phi22.csv')
    open(unit = 323, file = 'stderrors1/phi23.csv')
    open(unit = 324, file = 'stderrors1/phi24.csv')
    open(unit = 325, file = 'stderrors1/phi25.csv')
    open(unit = 326, file = 'stderrors1/phi26.csv')
    open(unit = 327, file = 'stderrors1/phi27.csv')
    
    open(unit = 330, file = 'stderrors1/phi30.csv')
    open(unit = 331, file = 'stderrors1/phi31.csv')
    open(unit = 332, file = 'stderrors1/phi32.csv')
    open(unit = 333, file = 'stderrors1/phi33.csv')
    open(unit = 334, file = 'stderrors1/phi34.csv')
    open(unit = 335, file = 'stderrors1/phi35.csv')
    open(unit = 336, file = 'stderrors1/phi36.csv')
    open(unit = 337, file = 'stderrors1/phi37.csv')
    
    open(unit = 340, file = 'stderrors1/phi40.csv')
    open(unit = 341, file = 'stderrors1/phi41.csv')
    open(unit = 342, file = 'stderrors1/phi42.csv')
    open(unit = 343, file = 'stderrors1/phi43.csv')
    open(unit = 344, file = 'stderrors1/phi44.csv')
    open(unit = 345, file = 'stderrors1/phi45.csv')
    open(unit = 346, file = 'stderrors1/phi46.csv')
    open(unit = 347, file = 'stderrors1/phi47.csv')
    
    open(unit = 350, file = 'stderrors1/phi50.csv')
    open(unit = 351, file = 'stderrors1/phi51.csv')
    open(unit = 352, file = 'stderrors1/phi52.csv')
    open(unit = 353, file = 'stderrors1/phi53.csv')
    open(unit = 354, file = 'stderrors1/phi54.csv')
    open(unit = 355, file = 'stderrors1/phi55.csv')
    open(unit = 356, file = 'stderrors1/phi56.csv')
    open(unit = 357, file = 'stderrors1/phi57.csv')
    
    open(unit = 360, file = 'stderrors1/phi60.csv')
    open(unit = 361, file = 'stderrors1/phi61.csv')
    open(unit = 362, file = 'stderrors1/phi62.csv')
    open(unit = 363, file = 'stderrors1/phi63.csv')
    open(unit = 364, file = 'stderrors1/phi64.csv')
    open(unit = 365, file = 'stderrors1/phi65.csv')
    open(unit = 366, file = 'stderrors1/phi66.csv')
    open(unit = 367, file = 'stderrors1/phi67.csv')
    
    open(unit = 370, file = 'stderrors1/phi70.csv')
    open(unit = 371, file = 'stderrors1/phi71.csv')
    open(unit = 372, file = 'stderrors1/phi72.csv')
    open(unit = 373, file = 'stderrors1/phi73.csv')
    open(unit = 374, file = 'stderrors1/phi74.csv')
    open(unit = 375, file = 'stderrors1/phi75.csv')
    open(unit = 376, file = 'stderrors1/phi76.csv')
    open(unit = 377, file = 'stderrors1/phi77.csv')
    
    open(unit = 401, file = 'stderrors1/sigma.csv')
    
    open(unit = 501, file = 'stderrors1/wagedifsigma.csv')
    
    open(unit = 600, file = 'stderrors1/xsi1998_0.csv')
    open(unit = 601, file = 'stderrors1/xsi1998_1.csv')
    open(unit = 602, file = 'stderrors1/xsi1998_2.csv')
    open(unit = 603, file = 'stderrors1/xsi1998_3.csv')
    open(unit = 604, file = 'stderrors1/xsi1998_4.csv')
    open(unit = 605, file = 'stderrors1/xsi1998_5.csv')
    open(unit = 606, file = 'stderrors1/xsi1998_6.csv')
    open(unit = 607, file = 'stderrors1/xsi1998_7.csv')
    
    open(unit = 700, file = 'stderrors1/xsi2000_0.csv')
    open(unit = 701, file = 'stderrors1/xsi2000_1.csv')
    open(unit = 702, file = 'stderrors1/xsi2000_2.csv')
    open(unit = 703, file = 'stderrors1/xsi2000_3.csv')
    open(unit = 704, file = 'stderrors1/xsi2000_4.csv')
    open(unit = 705, file = 'stderrors1/xsi2000_5.csv')
    open(unit = 706, file = 'stderrors1/xsi2000_6.csv')
    open(unit = 707, file = 'stderrors1/xsi2000_7.csv')
    
    open(unit = 800, file = 'stderrors1/xsi2005_0.csv')
    open(unit = 801, file = 'stderrors1/xsi2005_1.csv')
    open(unit = 802, file = 'stderrors1/xsi2005_2.csv')
    open(unit = 803, file = 'stderrors1/xsi2005_3.csv')
    open(unit = 804, file = 'stderrors1/xsi2005_4.csv')
    open(unit = 805, file = 'stderrors1/xsi2005_5.csv')
    open(unit = 806, file = 'stderrors1/xsi2005_6.csv')
    open(unit = 807, file = 'stderrors1/xsi2005_7.csv')
    
    open(unit = 900, file = 'stderrors1/eta0.csv')
    open(unit = 901, file = 'stderrors1/eta1.csv')
    open(unit = 902, file = 'stderrors1/eta2.csv')
    open(unit = 903, file = 'stderrors1/eta3.csv')
    open(unit = 904, file = 'stderrors1/eta4.csv')
    open(unit = 905, file = 'stderrors1/eta5.csv')
    open(unit = 906, file = 'stderrors1/eta6.csv')
    open(unit = 907, file = 'stderrors1/eta7.csv')
    
    open(unit = 1000, file = 'stderrors1/rho0.csv')
    open(unit = 1001, file = 'stderrors1/rho1.csv')
    open(unit = 1002, file = 'stderrors1/rho2.csv')
    open(unit = 1003, file = 'stderrors1/rho3.csv')
    open(unit = 1004, file = 'stderrors1/rho4.csv')
    open(unit = 1005, file = 'stderrors1/rho5.csv')
    open(unit = 1006, file = 'stderrors1/rho6.csv')
    open(unit = 1007, file = 'stderrors1/rho7.csv')
    
end if

count1 = size(theta)
Call MPI_BCAST(theta,count1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
Call MPI_BCAST(low_theta,count1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
Call MPI_BCAST(up_theta,count1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

count1 = size(beta)
Call MPI_BCAST(beta,count1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
Call MPI_BCAST(low_beta,count1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
Call MPI_BCAST(up_beta,count1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

count1 = size(sigma)
Call MPI_BCAST(sigma,count1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
Call MPI_BCAST(low_sigma,count1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
Call MPI_BCAST(up_sigma,count1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

count1 = size(kappa)
Call MPI_BCAST(kappa,count1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
Call MPI_BCAST(low_kappa,count1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
Call MPI_BCAST(up_kappa,count1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

count1 = size(tau)
Call MPI_BCAST(tau,count1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
Call MPI_BCAST(low_tau,count1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
Call MPI_BCAST(up_tau,count1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

Call MPI_BCAST(SigmaPref,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
Call MPI_BCAST(low_SigmaPref,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
Call MPI_BCAST(up_SigmaPref,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

count1 = size(omega)
Call MPI_BCAST(omega,count1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
Call MPI_BCAST(low_omega,count1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
Call MPI_BCAST(up_omega,count1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

count1 = size(gamma)
Call MPI_BCAST(gamma,count1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
Call MPI_BCAST(low_gamma,count1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
Call MPI_BCAST(up_gamma,count1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

count1 = size(lambda)
Call MPI_BCAST(lambda,count1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
Call MPI_BCAST(low_lambda,count1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
Call MPI_BCAST(up_lambda,count1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

count1 = size(aa)
Call MPI_BCAST(aa,count1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
Call MPI_BCAST(low_aa,count1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
Call MPI_BCAST(up_aa,count1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

count1 = size(sigma_prod)
Call MPI_BCAST(sigma_prod,count1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
Call MPI_BCAST(low_sigma_prod,count1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
Call MPI_BCAST(up_sigma_prod,count1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

! theta Parameters

do i = 1, 7
    param = param_in
    do n = 0, NGRID
        param(i) = low_theta(i) + ((theta(i) - low_theta(i))/real(NGRID,DOUBLE))*n
        PARAM_NUMBER = i
        PARAM_VALUE = param(i)
        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
        Call SMM_OBJ_FCN(param,F)        
    end do
    do n = 1, NGRID
        param(i) = theta(i) + ((up_theta(i) - theta(i))/real(NGRID,DOUBLE))*n
        PARAM_NUMBER = i
        PARAM_VALUE = param(i)
        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
        Call SMM_OBJ_FCN(param,F)
    end do
end do




! beta parameters

do s = 1, NSECTORS
    do i = 1, 12
        param = param_in
        do n = 0, NGRID
            param(7+(s-1)*12+i) = low_beta(s,i) + ((beta(s,i) - low_beta(s,i))/real(NGRID,DOUBLE))*n
            PARAM_NUMBER = 7+(s-1)*12+i
            PARAM_VALUE = param(7+(s-1)*12+i)
            CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
            Call SMM_OBJ_FCN(param,F)
        end do
        do n = 1, NGRID
            param(7+(s-1)*12+i) = beta(s,i) + ((up_beta(s,i) - beta(s,i))/real(NGRID,DOUBLE))*n
            PARAM_NUMBER = 7+(s-1)*12+i
            PARAM_VALUE = param(7+(s-1)*12+i)
            CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
            Call SMM_OBJ_FCN(param,F)
        end do
    end do
end do    



! sigma parameters

do s = 0, NSECTORS
    param = param_in
    do n = 0, NGRID
        param(92+s) = low_sigma(s) + ((sigma(s) - low_sigma(s))/real(NGRID,DOUBLE))*n
        PARAM_NUMBER = 92+s
        PARAM_VALUE = param(92+s)
        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
        Call SMM_OBJ_FCN(param,F)
    end do
   do n = 1, NGRID
        param(92+s) = sigma(s) + ((up_sigma(s) - sigma(s))/real(NGRID,DOUBLE))*n
        PARAM_NUMBER = 92+s
        PARAM_VALUE = param(92+s)
        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
        Call SMM_OBJ_FCN(param,F)
    end do
end do


! kappa parameters

do i = 1, 26
    param = param_in
    do n = 0, NGRID
        param(99+i) = low_kappa(i) + ((kappa(i) - low_kappa(i))/real(NGRID,DOUBLE))*n
        PARAM_NUMBER = 99+i
        PARAM_VALUE = param(99+i)
        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
        Call SMM_OBJ_FCN(param,F)
    end do
    do n = 1, NGRID
        param(99+i) = kappa(i) + ((up_kappa(i) - kappa(i))/real(NGRID,DOUBLE))*n
        PARAM_NUMBER = 99+i
        PARAM_VALUE = param(99+i)
        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
        Call SMM_OBJ_FCN(param,F)
    end do
end do


! tau parameters

do i = 2, NSECTORS
    param = param_in
    do n = 0, NGRID
        param(125+i-1) = low_tau(i) + ((tau(i) - low_tau(i))/real(NGRID,DOUBLE))*n
        PARAM_NUMBER = 125+i-1
        PARAM_VALUE = param(125+i-1)
        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
        Call SMM_OBJ_FCN(param,F)
    end do
    do n = 1, NGRID
        param(125+i-1) = tau(i) + ((up_tau(i) - tau(i))/real(NGRID,DOUBLE))*n
        PARAM_NUMBER = 125+i-1
        PARAM_VALUE = param(125+i-1)
        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
        Call SMM_OBJ_FCN(param,F)
    end do
end do


! SigmaPref

param = param_in
do n = 0, NGRID
    param(132) = low_SigmaPref + (SigmaPref - low_SigmaPref)/real(NGRID,DOUBLE)*n
    PARAM_NUMBER = 132
    PARAM_VALUE = param(132)
    CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
    Call SMM_OBJ_FCN(param,F)
end do
do n = 1, NGRID
    param(132) = SigmaPref + ((up_SigmaPref - SigmaPref)/real(NGRID,DOUBLE))*n
    PARAM_NUMBER = 132
    PARAM_VALUE = param(132)
    CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
    Call SMM_OBJ_FCN(param,F)
end do


! omega parameters

do tp = 2, NTYPEs
    do s = 0, NSECTORS
        param = param_in
        do n = 0, NGRID
            param(132+(tp-2)*(1+NSECTORS)+s+1) = low_omega(tp,s) + ((omega(tp,s) - low_omega(tp,s))/real(NGRID,DOUBLE))*n
            PARAM_NUMBER = 132+(tp-2)*(1+NSECTORS)+s+1
            PARAM_VALUE = param(132+(tp-2)*(1+NSECTORS)+s+1)
            CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
            Call SMM_OBJ_FCN(param,F)
        end do
        do n = 1, NGRID
            param(132+(tp-2)*(1+NSECTORS)+s+1) = omega(tp,s) + ((up_omega(tp,s) - omega(tp,s))/real(NGRID,DOUBLE))*n
            PARAM_NUMBER = 132+(tp-2)*(1+NSECTORS)+s+1
            PARAM_VALUE = param(132+(tp-2)*(1+NSECTORS)+s+1)
            CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
            Call SMM_OBJ_FCN(param,F)
        end do
    end do
end do    


! lambda parameters

do tp = 2, NTYPES
    param = param_in
    do n = 0, NGRID
        param(148+tp-1) = low_lambda(tp) + ((lambda(tp) - low_lambda(tp))/real(NGRID,DOUBLE))*n
        PARAM_NUMBER = 148+tp-1
        PARAM_VALUE = param(148+tp-1)
        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
        Call SMM_OBJ_FCN(param,F)
    end do
    do n = 1, NGRID
        param(148+tp-1) = lambda(tp) + ((up_lambda(tp) - lambda(tp))/real(NGRID,DOUBLE))*n
        PARAM_NUMBER = 148+tp-1
        PARAM_VALUE = param(148+tp-1)
        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
        Call SMM_OBJ_FCN(param,F)
    end do
end do


! gamma parameters

do tp = 2, NTYPES
    do s = 0, NSECTORS+4
        param = param_in
        do n = 0, NGRID
            param(150+(tp-2)*12+s+1) = low_gamma(tp,s) + ((gamma(tp,s) - low_gamma(tp,s))/real(NGRID,DOUBLE))*n
            PARAM_NUMBER = 150+(tp-2)*12+s+1
            PARAM_VALUE = param(150+(tp-2)*12+s+1)
            CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
            Call SMM_OBJ_FCN(param,F)
        end do
        do n = 1, NGRID
            param(150+(tp-2)*12+s+1) = gamma(tp,s) + ((up_gamma(tp,s) - gamma(tp,s))/real(NGRID,DOUBLE))*n
            PARAM_NUMBER = 150+(tp-2)*12+s+1
            PARAM_VALUE = param(150+(tp-2)*12+s+1)
            CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
            Call SMM_OBJ_FCN(param,F)
        end do
    end do
end do   


! aa parameters

do ed = 0, 1
    do s = 1, NSECTORS
        do i = 0, 1
            param = param_in
            do n = 0, NGRID
                param(174+ed*2*NSECTORS+(s-1)*2+i+1) = low_aa(s,ed,i) + ((aa(s,ed,i) - low_aa(s,ed,i))/real(NGRID,DOUBLE))*n
                PARAM_NUMBER = 174+ed*2*NSECTORS+(s-1)*2+i+1
                PARAM_VALUE = param(174+ed*2*NSECTORS+(s-1)*2+i+1)
                CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
                Call SMM_OBJ_FCN(param,F)
            end do
            do n = 1, NGRID
                param(174+ed*2*NSECTORS+(s-1)*2+i+1) = aa(s,ed,i) + ((up_aa(s,ed,i) - aa(s,ed,i))/real(NGRID,DOUBLE))*n
                PARAM_NUMBER = 174+ed*2*NSECTORS+(s-1)*2+i+1
                PARAM_VALUE = param(174+ed*2*NSECTORS+(s-1)*2+i+1)
                CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
                Call SMM_OBJ_FCN(param,F)
            end do
        end do
    end do
end do 

! sigma_prod parameters

do i = 1, NSECTORS
    param = param_in
    do n = 0, NGRID
        param(202+i) = low_sigma_prod(i) + ((sigma_prod(i) - low_sigma_prod(i))/real(NGRID,DOUBLE))*n
        PARAM_NUMBER = 202+i
        PARAM_VALUE = param(202+i)
        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
        Call SMM_OBJ_FCN(param,F)        
    end do
    do n = 1, NGRID
        param(202+i) = sigma_prod(i) + ((up_sigma_prod(i) - sigma_prod(i))/real(NGRID,DOUBLE))*n
        PARAM_NUMBER = 202+i
        PARAM_VALUE = param(202+i)
        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
        Call SMM_OBJ_FCN(param,F)
    end do
end do



if (rank == 0) then 

    close(101)
    close(102)
    close(103)
    close(104)
    close(105)
    close(106)
    close(107)
    
    close(200)
    close(201)
    close(202)
    close(203)
    close(204)
    close(205)
    close(206)
    close(207)
    
    close(300)
    close(301)
    close(302)
    close(303)
    close(304)
    close(305)
    close(305)
    close(307)
    
    close(310)
    close(311)
    close(312)
    close(313)
    close(314)
    close(315)
    close(316)
    close(317)
    
    close(320)
    close(321)
    close(322)
    close(323)
    close(324)
    close(325)
    close(326)
    close(327)
    
    close(330)
    close(331)
    close(332)
    close(333)
    close(334)
    close(335)
    close(336)
    close(337)
    
    close(340)
    close(341)
    close(342)
    close(343)
    close(344)
    close(345)
    close(346)
    close(347)
    
    close(350)
    close(351)
    close(352)
    close(353)
    close(354)
    close(355)
    close(356)
    close(357)
    
    close(360)
    close(361)
    close(362)
    close(363)
    close(364)
    close(365)
    close(366)
    close(367)
    
    close(370)
    close(371)
    close(372)
    close(373)
    close(374)
    close(375)
    close(376)
    close(377)
    
    close(401)
    
    close(501)
    
    close(600)
    close(601)
    close(602)
    close(603)
    close(604)
    close(605)
    close(606)
    close(607)
    
    close(700)
    close(701)
    close(702)
    close(703)
    close(704)
    close(705)
    close(706)
    close(707)
    
    close(800)
    close(801)
    close(802)
    close(803)
    close(804)
    close(805)
    close(806)
    close(807)
    
    close(900)
    close(901)
    close(902)
    close(903)
    close(904)
    close(905)
    close(906)
    close(907)
    
    close(1000)
    close(1001)
    close(1002)
    close(1003)
    close(1004)
    close(1005)
    close(1006)
    close(1007)

end if


! Compute the numerical derivatives

allocate(G0(dim_G0,NPARAM))
G0 = 0.0

if (rank == 0) then

! wage regressions

start = 0

do eq = 1, NSECTORS
    
    write(char,9999) eq
    filename = 'stderrors1/beta'
    filename = filename(1:15) // char // '.csv'
    
    open(unit = 1, file = filename)
        i = 1
        do 
            read(1,3000,IOSTAT=status) npar, X(npar,i), (delta(npar,i,j), j = 1, NREG_Wage)
            i = i + 1
            if (i == 2*NGRID+2) then
                i = 1
            end if
            if (status /= 0 ) exit
        end do
    close(1)
    
    do n = 1, NPARAM
        do j = 1, NREG_Wage
            Y(:,1) = delta(n,:,j)
            XR(:,1) = 1.0
            XR(:,2) = X(n,:)
            XR(:,3) = X(n,:)**2
            Call LinReg(Y,XR,2*NGRID+1,3,COEF,SST,SSE,info)
            if (info /= 0) then
                print*, 'Problem Running Regression'
                print*, 'Wage', eq
                pause
                stop
            end if
            G0(start+(eq-1)*NREG_Wage+j,n) = COEF(2,1) + 2*COEF(3,1)*param_in(n)
        end do
    end do

end do



! sectoral choice regressions

start = start + NSECTORS*NREG_Wage

do eq = 0, NSECTORS
    
    write(char,9999) eq
    filename = 'stderrors1/gamma'
    filename = filename(1:16) // char // '.csv'
    
    open(unit = 1, file = filename)
        i = 1
        do 
            read(1,3000,IOSTAT=status) npar, X(npar,i), (delta(npar,i,j), j = 1, NREG_Emp)
            i = i + 1
            if (i == 2*NGRID+2) then
                i = 1
            end if
            if (status /= 0 ) exit
        end do
    close(1)
    
    do n = 1, NPARAM
        do j = 1, NREG_Emp
            Y(:,1) = delta(n,:,j)
            XR(:,1) = 1.0
            XR(:,2) = X(n,:)
            XR(:,3) = X(n,:)**2
            Call LinReg(Y,XR,2*NGRID+1,3,COEF,SST,SSE,info)
            if (info /= 0) then
                print*, 'Problem Running Regression'
                print*, 'Employment', eq
                pause
                stop
            end if
            G0(start+eq*NREG_Emp+j,n) = COEF(2,1) + 2*COEF(3,1)*param_in(n)
        end do
    end do

end do



! transition rate regressions

start = start + (1+NSECTORS)*NREG_Emp

do s1 = 0, NSECTORS
    do s2 = 0, NSECTORS
    
        write(char1,9999) s1
        write(char2,9999) s2
        filename = 'stderrors1/phi'
        filename = filename(1:14) // char1 // char2 // '.csv'
    
        open(unit = 1, file = filename)
            i = 1
            do 
                read(1,3000,IOSTAT=status) npar, X(npar,i), (delta(npar,i,j), j = 1, NREG_Tr)
                i = i + 1
                if (i == 2*NGRID+2) then
                    i = 1
                end if
                if (status /= 0 ) exit
            end do
        close(1)
    
        do n = 1, NPARAM
            do j = 1, NREG_Tr
                Y(:,1) = delta(n,:,j)
                XR(:,1) = 1.0
                XR(:,2) = X(n,:)
                XR(:,3) = X(n,:)**2
                Call LinReg(Y,XR,2*NGRID+1,3,COEF,SST,SSE,info)
                if (info /= 0) then
                    print*, 'Problem Running Regression'
                    print*, 'Tr', s1, s2
                    pause
                    stop
                end if
                G0(start+s1*(1+NSECTORS)*NREG_Tr+s2*NREG_Tr+j,n) = COEF(2,1) + 2*COEF(3,1)*param_in(n)
            end do
        end do

    end do
end do

! residuals of the wage regressions

start = start + (1+NSECTORS)*(1+NSECTORS)*NREG_Tr

open(unit = 1, file = 'stderrors1/sigma.csv')
    i = 1
    do 
        read(1,4000,IOSTAT=status) npar, X(npar,i), (delta(npar,i,j), j = 1, NSECTORS)
        i = i + 1
        if (i == 2*NGRID+2) then
            i = 1
        end if
        if (status /= 0 ) exit
    end do
close(1)
    
do n = 1, NPARAM
    do j = 1, NSECTORS
        Y(:,1) = delta(n,:,j)
        XR(:,1) = 1.0
        XR(:,2) = X(n,:)
        XR(:,3) = X(n,:)**2
        Call LinReg(Y,XR,2*NGRID+1,3,COEF,SST,SSE,info)
        if (info /= 0) then
            print*, 'Problem Running Regression'
            print*, 'Sigma'
            pause
            stop
        end if
        G0(start+j,n) = COEF(2,1) + 2*COEF(3,1)*param_in(n)
    end do
end do

! residuals of the within individual wage changes

start = start + NSECTORS

open(unit = 1, file = 'stderrors1/WageDifSigma.csv')
    i = 1
    do 
        read(1,4000,IOSTAT=status) npar, X(npar,i), (delta(npar,i,j), j = 1, NSECTORS)
        i = i + 1
        if (i == 2*NGRID+2) then
            i = 1
        end if
        if (status /= 0 ) exit
    end do
close(1)
    
do n = 1, NPARAM
    do j = 1, NSECTORS
        Y(:,1) = delta(n,:,j)
        XR(:,1) = 1.0
        XR(:,2) = X(n,:)
        XR(:,3) = X(n,:)**2
        Call LinReg(Y,XR,2*NGRID+1,3,COEF,SST,SSE,info)
        if (info /= 0) then
            print*, 'Problem Running Regression'
            print*, 'WageDifSigma'
            pause
            stop
        end if
        G0(start+j, n) = COEF(2,1) + 2*COEF(3,1)*param_in(n)
    end do
end do


! Persistence 1998 regressions

start = start + NSECTORS

do eq = 0, NSECTORS

    write(char,9999) eq
    filename = 'stderrors1/xsi1998_'
    filename = filename(1:19) // char // '.csv'
    
    open(unit = 1, file = filename)
        i = 1
        do 
            read(1,3001,IOSTAT=status) npar, X(npar,i), (delta(npar,i,j), j = 1, NREG_Pers)
            i = i + 1
            if (i == 2*NGRID+2) then
                i = 1
            end if
            if (status /= 0 ) exit
        end do
    close(1)
    
    do n = 1, NPARAM
        do j = 1, NREG_Pers
            Y(:,1) = delta(n,:,j)
            XR(:,1) = 1.0
            XR(:,2) = X(n,:)
            XR(:,3) = X(n,:)**2
            Call LinReg(Y,XR,2*NGRID+1,3,COEF,SST,SSE,info)
            if (info /= 0) then
                print*, 'Problem Running Regression'
                print*, 'Pers1998', eq
                pause
                stop
            end if
            G0(start+eq*NREG_Pers+j,n) = COEF(2,1) + 2*COEF(3,1)*param_in(n)
        end do
    end do

end do

! Persistence 2000 regressions

start = start + (1+NSECTORS)*NREG_Pers

do eq = 0, NSECTORS

    write(char,9999) eq
    filename = 'stderrors1/xsi2000_'
    filename = filename(1:19) // char // '.csv'
    
    open(unit = 1, file = filename)
        i = 1
        do 
            read(1,3001,IOSTAT=status) npar, X(npar,i), (delta(npar,i,j), j = 1, NREG_Pers)
            i = i + 1
            if (i == 2*NGRID+2) then
                i = 1
            end if
            if (status /= 0 ) exit
        end do
    close(1)
    
    do n = 1, NPARAM
        do j = 1, NREG_Pers
            Y(:,1) = delta(n,:,j)
            XR(:,1) = 1.0
            XR(:,2) = X(n,:)
            XR(:,3) = X(n,:)**2
            Call LinReg(Y,XR,2*NGRID+1,3,COEF,SST,SSE,info)
            if (info /= 0) then
                print*, 'Problem Running Regression'
                print*, 'Pers2000', eq
                pause
                stop
            end if
            G0(start+eq*NREG_Pers+j,n) = COEF(2,1) + 2*COEF(3,1)*param_in(n)
        end do
    end do

end do


! Persistence 2005 regressions

start = start + (1+NSECTORS)*NREG_Pers

do eq = 0, NSECTORS

    write(char,9999) eq
    filename = 'stderrors1/xsi2005_'
    filename = filename(1:19) // char // '.csv'
    
    open(unit = 1, file = filename)
        i = 1
        do 
            read(1,3001,IOSTAT=status) npar, X(npar,i), (delta(npar,i,j), j = 1, NREG_Pers)
            i = i + 1
            if (i == 2*NGRID+2) then
                i = 1
            end if
            if (status /= 0 ) exit
        end do
    close(1)
    
    do n = 1, NPARAM
        do j = 1, NREG_Pers
            Y(:,1) = delta(n,:,j)
            XR(:,1) = 1.0
            XR(:,2) = X(n,:)
            XR(:,3) = X(n,:)**2
            Call LinReg(Y,XR,2*NGRID+1,3,COEF,SST,SSE,info)
            if (info /= 0) then
                print*, 'Problem Running Regression'
                print*, 'Pers2005', eq
                pause
                stop
            end if
            G0(start+eq*NREG_Pers+j,n) = COEF(2,1) + 2*COEF(3,1)*param_in(n)
        end do
    end do

end do


! Frequency regressions

start = start + (1+NSECTORS)*NREG_Pers

do eq = 0, NSECTORS

    write(char,9999) eq
    filename = 'stderrors1/eta'
    filename = filename(1:14) // char // '.csv'
    
    open(unit = 1, file = filename)
        i = 1
        do 
            read(1,3001,IOSTAT=status) npar, X(npar,i), (delta(npar,i,j), j = 1, NREG_Freq)
            i = i + 1
            if (i == 2*NGRID+2) then
                i = 1
            end if
            if (status /= 0 ) exit
        end do
    close(1)
    
    do n = 1, NPARAM
        do j = 1, NREG_Freq
            Y(:,1) = delta(n,:,j)
            XR(:,1) = 1.0
            XR(:,2) = X(n,:)
            XR(:,3) = X(n,:)**2
            Call LinReg(Y,XR,2*NGRID+1,3,COEF,SST,SSE,info)
            if (info /= 0) then
                print*, 'Problem Running Regression'
                print*, 'Freq', eq
                pause
                stop
            end if
            G0(start+eq*NREG_Freq+j,n) = COEF(2,1) + 2*COEF(3,1)*param_in(n)
        end do
    end do

end do


! Return regressions

start = start + (1+NSECTORS)*NREG_Freq

do eq = 0, NSECTORS

    write(char,9999) eq
    filename = 'stderrors1/rho'
    filename = filename(1:14) // char // '.csv'
    
    open(unit = 1, file = filename)
        i = 1
        do 
            read(1,3000,IOSTAT=status) npar, X(npar,i), (delta(npar,i,j), j = 1, NREG_Return)
            i = i + 1
            if (i == 2*NGRID+2) then
                i = 1
            end if
            if (status /= 0 ) exit
        end do
    close(1)
    
    do n = 1, NPARAM
        do j = 1, NREG_Return
            Y(:,1) = delta(n,:,j)
            XR(:,1) = 1.0
            XR(:,2) = X(n,:)
            XR(:,3) = X(n,:)**2
            Call LinReg(Y,XR,2*NGRID+1,3,COEF,SST,SSE,info)
            if (info /= 0) then
                print*, 'Problem Running Regression'
                print*, 'Return', eq
                pause
                stop
            end if
            G0(start+eq*NREG_Return+j,n) = COEF(2,1) + 2*COEF(3,1)*param_in(n)
        end do
    end do

end do



! File with the numerical derivatives

open(unit = 1, file = 'stderrors1/G0.csv')    

    start = 0
    ! wage regressions
    do eq = 1, NSECTORS
        do j = 1, NREG_Wage
            write(1,5000) start+(eq-1)*NREG_Wage+j, 'Wage', eq, j, (G0(start+(eq-1)*NREG_Wage+j,n), n=1,NPARAM)
        end do
    end do
    
    start = start + NSECTORS*NREG_Wage  
    ! sectoral choice regressions
    do eq = 0, NSECTORS
        do j = 1, NREG_Emp
            write(1,5000) start+eq*NREG_Emp+j, 'Emp', eq, j, (G0(start+eq*NREG_Emp+j,n), n=1,NPARAM)
        end do
    end do
    
    start = start + (1+NSECTORS)*NREG_Emp
    ! transition rate regressions
    do s1 = 0, NSECTORS
        do s2 = 0, NSECTORS
            do j = 1, NREG_Tr
                write(1,5000) start+s1*5*NREG_Tr+s2*NREG_Tr+j, 'Tr', s1*(1+NSECTORS)+s2, j, (G0(start+s1*5*NREG_Tr+s2*NREG_Tr+j,n), n=1,NPARAM)
            end do
        end do
    end do
    
    
    start = start + (1+NSECTORS)*(1+NSECTORS)*NREG_Tr
    ! residuals of wage regressions
    do eq = 1, NSECTORS
        write(1,5000) start+eq, 'Sigma', eq, 1, (G0(start+eq,n), n=1,NPARAM)
    end do
    
    start = start + NSECTORS
    ! residuals of within individual wage changes
    do eq = 1, NSECTORS
        write(1,5000) start+eq, 'WageDifSigma', eq, 1, (G0(start+eq,n), n=1,NPARAM)
    end do
    
    start = start + NSECTORS
    ! persistence regressions 1998
    do eq = 0, NSECTORS
        do j = 1, NREG_Pers
            write(1,5000) start+eq*NREG_Pers+j, 'Pers1998', eq, j, (G0(start+eq*NREG_Pers+j,n), n=1,NPARAM)
        end do
    end do
    
    start = start + (1+NSECTORS)*NREG_Pers
    ! persistence regressions 2000
    do eq = 0, NSECTORS
        do j = 1, NREG_Pers
            write(1,5000) start+eq*NREG_Pers+j, 'Pers2000', eq, j, (G0(start+eq*NREG_Pers+j,n), n=1,NPARAM)
        end do
    end do
    
    start = start + (1+NSECTORS)*NREG_Pers
    ! persistence regressions 2005
    do eq = 0, NSECTORS
        do j = 1, NREG_Pers
            write(1,5000) start+eq*NREG_Pers+j, 'Pers2005', eq, j, (G0(start+eq*NREG_Pers+j,n), n=1,NPARAM)
        end do
    end do
    
    start = start + (1+NSECTORS)*NREG_Pers
    ! Frequency regressions
    do eq = 0, NSECTORS
        do j = 1, NREG_Freq
            write(1,5000) start+eq*NREG_Freq+j, 'Freq', eq, j, (G0(start+eq*NREG_Freq+j,n), n=1,NPARAM)
        end do
    end do
    
    start = start + (1+NSECTORS)*NREG_Pers
    ! Return regressions
    do eq = 0, NSECTORS
        do j = 1, NREG_Return
            write(1,5000) start+eq*NREG_Return+j, 'Return', eq, j, (G0(start+eq*NREG_Return+j,n), n=1,NPARAM)
        end do
    end do
    

close(1)

deallocate(G0)

end if

Call MPI_FINALIZE(ierr)

9999 format(i1)
3000 format(i4,',',25(f20.10,','))
3001 format(i4,',',22(f20.10,','))
4000 format(i4,',',8(f20.10,','))
5000 format(i5,',',a15,',',2(i5,','), 209(f20.10,','))

End Subroutine NumDer



Subroutine ComputeStdErrors

! This subroutine reads all the matrices relevant for the computation
! Of the standard errors

USE Global_Data

implicit none


real(KIND=DOUBLE), allocatable, dimension(:,:) :: G0, G0_T, BOOT, BOOTSIM, Omega, &
                                                  MAT1, MAT1_T, MAT2, invMAT2, MAT3, &
                                                  MAT4, MAT5, MAT6, U, COVAR

real(KIND=DOUBLE) STD_ERRORS(NPARAM)
                  

integer status, blah, i, j, k, s, s1, s2, info, start, dim_G0


dim_G0 = NSECTORS*NREG_Wage + & 
        (1+NSECTORS)*NREG_Emp + &
        (1+NSECTORS)*(1+NSECTORS)*NREG_Tr + &
        NSECTORS + &
        NSECTORS + &
        (1+NSECTORS)*NREG_Pers*3 + &
        (1+NSECTORS)*NREG_Freq + &
        (1+NSECTORS)*NREG_Return  


allocate(G0(dim_G0,NPARAM)) 
allocate(G0_T(NPARAM,dim_G0))
allocate(BOOT(dim_G0,dim_G0)) 
allocate(BOOTSIM(dim_G0,dim_G0)) 
allocate(Omega(dim_G0,dim_G0))
allocate(MAT1(NPARAM,dim_G0)) 
allocate(MAT1_T(dim_G0,NPARAM)) 
allocate(MAT2(NPARAM,NPARAM))
allocate(invMAT2(NPARAM,NPARAM))
allocate(MAT3(dim_G0,dim_G0))
allocate(MAT4(NPARAM,dim_G0))
allocate(MAT5(dim_G0,dim_G0))
allocate(MAT6(NPARAM,NPARAM))
allocate(U(NPARAM,NPARAM))
allocate(COVAR(NPARAM,NPARAM))
    
! *****************************************************************************
! Reading the matrices that will be used in
! The computation of standard errors
! COVAR = inv(G0'*Omega*G0)*G0'*Omega*(BOOT+BOOTSIM)*Omega*G0*inv(G0'*Omega*G0)
! *****************************************************************************

! ***************************
! Omega, the weighting Matrix
! ***************************

Omega = 0.0

start = 0
do s = 1, NSECTORS
    do i = 1, NREG_Wage
        do j = 1, NREG_Wage
            Omega((s-1)*NREG_Wage+i,(s-1)*NREG_Wage+j) = invCOVbetaData(s,i,j)
        end do
    end do
end do  

start = start + NSECTORS*NREG_Wage
do s = 0, NSECTORS
    do i = 1, NREG_Emp
        do j = 1, NREG_Emp
            Omega(start+s*NREG_Emp+i,start+s*NREG_Emp+j) = invCOVgammaData(s,i,j)
        end do
    end do
end do

start = start + (1+NSECTORS)*NREG_Emp
do s1 = 0, NSECTORS
    do s2 = 0, NSECTORS
        do i = 1, 21
            do j = 1, 21
                Omega(start+s1*5*NREG_Tr+s2*NREG_Tr+i,&
                      start+s1*5*NREG_Tr+s2*NREG_Tr+j) = invCOVphiData(s1,s2,i,j)
            end do
        end do
    end do
end do     

start = start + (1+NSECTORS)*(1+NSECTORS)*NREG_Tr
do s = 1, NSECTORS
    Omega(start+s,start+s) = 1/VarSigma2Data(s)
end do

start = start + NSECTORS
do s = 1, NSECTORS
    Omega(start+s,start+s) = 1/VarSigma2WageDifData(s)
end do

start = start + NSECTORS
do s = 0, NSECTORS
    do i = 1, NREG_Pers
        do j = 1, NREG_Pers
            Omega(start+s*NREG_Pers+i,start+s*NREG_Pers+j) = invCOVxsi1998Data(s,i,j)
        end do
    end do
end do

start = start + (1+NSECTORS)*NREG_Pers
do s = 0, NSECTORS
    do i = 1, NREG_Pers
        do j = 1, NREG_Pers
            Omega(start+s*NREG_Pers+i,start+s*NREG_Pers+j) = invCOVxsi2000Data(s,i,j)
        end do
    end do
end do

start = start + (1+NSECTORS)*NREG_Pers
do s = 0, NSECTORS
    do i = 1, NREG_Pers
        do j = 1, NREG_Pers
            Omega(start+s*NREG_Pers+i,start+s*NREG_Pers+j) = invCOVxsi2005Data(s,i,j)
        end do
    end do
end do

start = start + (1+NSECTORS)*NREG_Pers
do s = 0, NSECTORS
    do i = 1, NREG_Freq
        do j = 1, NREG_Freq
            Omega(start+s*NREG_Freq+i,start+s*NREG_Freq+j) = invCOVetaData(s,i,j)
        end do
    end do
end do

start = start + (1+NSECTORS)*NREG_Freq
do s = 0, NSECTORS
    do i = 1, NREG_Return
        do j = 1, NREG_Return
            Omega(start+s*NREG_Return+i,start+s*NREG_Return+j) = invCOVrhoData(s,i,j)
        end do
    end do
end do


! *******************************
! G0, the matrix with derivatives
! *******************************

open (unit = 1, file = 'stderrors1/G0.csv')
    do
        read(1,*,IOSTAT=status) i, blah, blah, (G0(i,j), j = 1, NPARAM)
        if (status /= 0 ) exit
    end do
close(1)



! **************************************************
! BOOT, the variance of delta, computed by bootstrap
! **************************************************

open(unit = 1, file = 'stderrors1/boot_covar.csv')
    read(1,*)
    do i = 1, dim_G0
        read(1,*) (BOOT(i,j), j = 1, dim_G0)
    end do
close(1)


! *******************************************************
! BOOTSIM, the variance of delta_S, computed by bootstrap
! *******************************************************

open(unit = 1, file = 'stderrors1/boot_covar_sim.csv')
    read(1,*)
    do i = 1, dim_G0
        read(1,*) (BOOTSIM(i,j), j = 1, dim_G0)
    end do
close(1)



! ************
! G0 transpose
! ************

do i = 1, NPARAM 
    do j = 1, dim_G0
        G0_T(i,j) = G0(j,i)
    end do
end do    


! ************************
! Compute MAT1 = G0'*Omega
! ************************

MAT1 = 0.0
do i = 1, NPARAM
    do j = 1, dim_G0
        do k = 1, dim_G0
            MAT1(i,j) = MAT1(i,j) + G0_T(i,k)*Omega(k,j)
        end do    
    end do
end do


! *************************************
! Compute MAT2 = MAT1*G0 = G0'*Omega*G0
! *************************************

MAT2 = 0.0
do i = 1, NPARAM
    do j = 1, NPARAM
        do k = 1, dim_G0
            MAT2(i,j) = MAT2(i,j) + MAT1(i,k)*G0(k,j)
        end do
    end do
end do  

! *****************************
! Inverting MAT2 = G0'*Omega*G0
! *****************************

U = MAT2
    
call dpotrf('U',NPARAM,U,NPARAM,info)    
    
invMAT2 = U
    
call dpotri('U',NPARAM,invMAT2,NPARAM,info)
    
do i = 1, NPARAM
    do j = 1, NPARAM
        if (i > j) then
            invMAT2(i,j) = invMAT2(j,i)
        end if
    end do
end do



! ******************************
! MAT3 = Var_delta + Var_delta_S         
! ******************************

MAT3 = BOOT + BOOTSIM


! *************************
! MAT1_T = MAT1' = Omega*G0
! *************************

do i = 1, dim_G0
    do j = 1, NPARAM
        MAT1_T(i,j) = MAT1(j,i)
    end do
end do    

! *******************************************
! MAT4 = MAT1*MAT3 = G0'*Omega*(BOOT+BOOTSIM)
! *******************************************

MAT4 = 0.0
do i = 1, NPARAM
    do j = 1, dim_G0
        do k = 1, dim_G0
            MAT4(i,j) = MAT4(i,j) + MAT1(i,k)*MAT3(k,j)
        end do
    end do
end do    



! ********************************************************
! MAT5 = MAT4*Omega*G0 = G0'*Omega*(BOOT+BOOTSIM)*Omega*G0
! ********************************************************

MAT5 = 0.0
do i = 1, NPARAM
    do j = 1, NPARAM
        do k = 1, dim_G0
            MAT5(i,j) = MAT5(i,j) + MAT4(i,k)*MAT1_T(k,j)
        end do
    end do
end do     



! ***************************************************************************
! MAT6 = inv(MAT2)*MAT5 = inv(G0'*Omega*G0)*G0'*Omega*(BOOT+BOOTSIM)*Omega*G0
! ***************************************************************************

MAT6 = 0.0
do i = 1, NPARAM
    do j = 1, NPARAM
        do k = 1, NPARAM
            MAT6(i,j) = MAT6(i,j) + invMAT2(i,k)*MAT5(k,j)
        end do
    end do
end do


! **********************************************************************************************
! COVAR = MAT6*inv(MAT2) = inv(G0'*Omega*G0)*G0'*Omega*(BOOT+BOOTSIM)*Omega*G0*inv(G0'*Omega*G0)
! **********************************************************************************************

COVAR = 0.0
do i = 1, NPARAM
    do j = 1, NPARAM
        do k = 1, NPARAM
            COVAR(i,j) = COVAR(i,j) + MAT6(i,k)*invMAT2(k,j)
        end do
    end do
end do


open(unit = 1, file = 'stderrors1/COVAR.csv')
    do i = 1, NPARAM
        write(1,1001) (COVAR(i,j), j=1,NPARAM)
    end do
close(1)


! ***************
! Standard Errors
! ***************

do i = 1, NPARAM
    STD_ERRORS(i) = sqrt(COVAR(i,i))
end do

open(unit = 1, file = 'stderrors1/STD_ERRORS.csv')
    do i = 1, NPARAM
        write(1,*) STD_ERRORS(i)
    end do
close(1)

deallocate(G0, G0_T, BOOT, BOOTSIM, Omega, MAT1, MAT1_T, MAT2, invMAT2, MAT3, &
           MAT4, MAT5, MAT6, U, COVAR)

1001 format(209(f36.15,','))

end subroutine ComputeStdErrors



End Module StdErrors_MOD
