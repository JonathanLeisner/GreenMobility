For questions regarding these codes, please contact
Rafael Dix-Carneiro at rafael.dix.carneiro@duke.edu or
dix-carneiro@econ.umd.edu


================================
MPI/Fortran Codes for estimation
================================

Main.f90 (Main program, reads coefficients from auxiliary models and other data)
Global_Data.f90 (Global variables are defined)
Loss_Function_MOD.f90 (Simulate sthe economy and computes the indirect inference loss function)
Parallel_Emax_MOD.f90 (Computes the Emax function)
ParallelCohorts_MOD.f90 (Computes aggregate human capital supply)
NewtonSolver.f90 (Nonlinear solver)
LinReg_MOD.f90 (Linear regression routine)

Nelder-Mead Optimization Routine
minim.f90 

NEWUOA Optimization Routine
newoa.f
newob.f
bigden.f
biglag.f
calfun.f
trsapp.f
update.f

=========================================
MPI and Intel Fortran Compilation command
=========================================
mpif90 NewtonSolver.f90 Global_Data.f90 LinReg_MOD.f90 Parallel_Emax_MOD.f90 ParallelCohorts_MOD.f90 
Loss_Function_MOD.f90 bigden.f biglag.f calfun.f newuoa.f newuob.f trsapp.f update.f minim.f90 
Main.f90 -o estimation -L$LIBRARY_PATH -I$INCLUDE -Wl,--start-group -lmkl_intel_lp64 -lmkl_intel_thread 
-lmkl_core -Wl,--end-group -liomp5 -lpthread -openmp -heap-arrays

========
PBS file
========
#PBS -N Job_Estimation

#PBS -o PBS_out

#PBS -e PBS_error

#PBS -q special

#PBS -l nodes=2:ppn=12,walltime = 72:00:00

#PBS -m abe

#PBS -M dix-carneiro@econ.umd.edu



export OMP_NUM_THREADS=1



cd path



mpiexec -np 24 ./estimation > estimation_out

==============
IMPORTANT NOTE
==============
The parallelization is set in a way that it will only work when the program is executed over 24 cores, 
as is illustrated in the PBS file above.

======
Output
======
The output of the program is illustrated below (one iteration of the optimization algorithm).

 ************************
 Computing Emax Myopic...
 ************************
 
 
 Emax Time     :   12.92188    
 
 ********************
 Emax Myopic Computed
 ********************
 
 
 ***************************************
 Computing Emax Rational Expectations...
 ***************************************
 
 
 SUM_FLAG_EMAX =           0
 
 
 Emax Time     :   100.6875    
 
 ***********************************
 Emax Rational Expectations Computed
 ***********************************
 
 
 Time Computing Equilibrium:      131.0273    
 
 
 Loss Function Wage:             13.4732315057877     
 Loss Function Emp:              6.71126339650336     
 Loss Function Tr:               24.1366654491354     
 Loss Function Return:           1.67202279849637     
 Loss Function Sigma:           0.110467126360453     
 Loss Function Sigma Dif:       6.215321503083153E-002
 Loss Function Pers 1998:        1.87648992524519     
 Loss Function Pers 2000:        1.59774716136451     
 Loss Function Pers 2005:        1.85462799168799     
 Loss Function Freq:             2.87286931531720     
 Loss Function Labor Shares:     6.74373350334372     
 Loss Function K Shares:         6.40024973015558     
 
 BATCH =           2
 
 ALGORITHM: NEWUOA
 
 FUNCTION ITERATION:            1
 
 theta:
 0.8458  -0.2898   0.3434   0.6503   1.6509   0.0472  -0.0014  
 
 wage equation 1
-0.3474   0.1591   0.9600   0.0369  -0.0009  
 0.1332   0.0534   0.0620   0.0406   0.0344   0.0648   0.0520  
 
 wage equation 2
-0.2864   0.3101   0.8023   0.0204  -0.0007  
 0.0226   0.1028   0.0636   0.0250   0.0479   0.0450   0.0583  
 
 wage equation 3
-0.3304   0.3177   0.7978   0.0203  -0.0006  
 0.0534   0.0423   0.0979   0.0398   0.0414   0.0449   0.0404  
 
 wage equation 4
-0.5155   0.2450   0.5529   0.0171  -0.0005  
 0.0114   0.0262   0.0492   0.0991   0.0205   0.0435   0.0453  
 
 wage equation 5
-0.2792   0.3665   0.7210   0.0105  -0.0003  
 0.0094   0.0337   0.0563   0.0101   0.0916   0.0369   0.0551  
 
 wage equation 6
-0.3509   0.3093   0.7202   0.0244  -0.0005  
 0.0678   0.0371   0.0657   0.0430   0.0572   0.1014   0.0533  
 
 wage equation 7
-0.2485   0.3074   0.8443   0.0189  -0.0004  
 0.0056  -0.0083   0.0205  -0.0008   0.0056   0.0145   0.0932  
 
 sigma:
18.48   0.20   0.17   0.15   0.19   0.20   0.17   0.22  
 
 cost:
 3.24   3.28   3.41   3.27   3.24   3.35   3.27  
 0.00  -0.07   0.21  -0.22  -0.05   0.18  -0.18  
 1.85   1.93   2.46   1.81   1.84   2.11   2.09  
 0.1130   0.0104   0.0467   0.2925   0.0262  -0.0006  
 
 lambda:
 0.0000  -0.1888   0.1885  
 
 tau:
 0.0000   0.0000   0.4641   0.4139  -0.0699   0.1174   0.2450   0.7142  
 
 preference shock:
   2.16  
 
 omega 2:
 0.8612   0.4986   0.7983   0.5795   0.8753   0.8356   0.5216   0.7881  
 
 omega 3:
-0.8163  -0.6406  -0.5304  -0.9979  -0.2922  -0.4150  -0.9111  -0.6568  
 
 gamma 2:
 0.2413   0.3588   0.0837   0.3276   0.1701   0.1006   0.0667   0.1945  
 0.0420   0.0224   0.5071   0.3379  
 
 gamma 3:
 0.5011   0.4156   0.1251   0.1571   0.0871   0.2211  -0.0574   0.2478  
 0.2489   0.1035   0.6521   0.1986  
 
 sigma_prod:
 0.12   0.07   0.45   0.08   0.29   0.39   0.12  
 
 alpha_prod 1:
 0.25   0.23  
 0.14   0.12  
 
 alpha_prod 2:
 0.45   0.33  
 0.31   0.29  
 
 alpha_prod 3:
 0.34   0.30  
 0.60   0.57  
 
 alpha_prod 4:
 0.34   0.34  
 0.17   0.14  
 
 alpha_prod 5:
 0.32   0.36  
 0.42   0.42  
 
 alpha_prod 6:
 0.41   0.27  
 0.52   0.50  
 
 alpha_prod 7:
 0.18   0.21  
 0.58   0.68  
 
 rsk_init Ed = 0:
 1.000   1.000   1.000   1.000   1.000   1.000   1.000  
 
 rsk_init Ed = 1:
 1.000   1.000   1.000   1.000   1.000   1.000   1.000  
 
 Skill Prices Ed = 0:
 1995   1.64   1.84   2.28   1.76   1.62   2.21   1.84  
 1996   1.57   1.82   2.13   1.76   1.55   2.17   1.94  
 1997   1.53   1.81   2.28   1.82   1.60   2.12   1.97  
 1998   1.42   1.71   2.13   1.78   1.55   2.07   1.93  
 1999   1.37   1.65   2.08   1.64   1.48   1.92   1.83  
 2000   1.44   1.69   2.20   1.66   1.53   1.96   1.79  
 2001   1.39   1.63   2.09   1.62   1.49   1.82   1.77  
 2002   1.34   1.53   2.00   1.53   1.41   1.75   1.72  
 2003   1.41   1.55   2.07   1.54   1.45   1.72   1.70  
 2004   1.37   1.65   2.07   1.61   1.50   1.77   1.69  
 2005   1.28   1.58   2.03   1.59   1.49   1.85   1.72  
 
 Skill Prices Ed = 1:
 1995   2.59   3.54   4.54   3.39   3.52   4.06   3.75  
 1996   2.59   3.45   4.27   3.45   3.43   4.06   3.96  
 1997   2.50   3.44   4.54   3.52   3.38   4.11   4.03  
 1998   2.47   3.28   4.19   3.42   3.26   4.01   3.98  
 1999   2.32   3.20   4.14   3.17   3.16   3.72   3.82  
 2000   2.37   3.31   4.34   3.12   3.21   3.86   3.70  
 2001   2.34   3.20   4.09   3.10   3.12   3.55   3.64  
 2002   2.31   2.97   3.90   2.86   2.89   3.45   3.51  
 2003   2.33   3.04   4.11   2.75   2.98   3.35   3.47  
 2004   2.29   3.17   4.13   2.91   2.98   3.46   3.45  
 2005   2.20   3.05   4.03   2.96   2.98   3.66   3.52  
 
 Iterations: 
          40          40          40          31          25          40
          27          40          40          40          17
 
 Check Ed = 0: 
 1995   0.0002   0.0001   0.0002   0.0001   0.0001   0.0001   0.0001  
 1996   0.0003   0.0002   0.0000   0.0000   0.0000   0.0002   0.0001  
 1997   0.0010   0.0000   0.0001   0.0008   0.0000   0.0011   0.0000  
 1998   0.0002   0.0007   0.0006   0.0003   0.0003   0.0000   0.0002  
 1999   0.0009   0.0006   0.0005   0.0002   0.0000   0.0000   0.0001  
 2000   0.0000   0.0007   0.0000   0.0000   0.0004   0.0000   0.0000  
 2001   0.0001   0.0001   0.0005   0.0002   0.0001   0.0002   0.0001  
 2002   0.0008   0.0000   0.0001   0.0001   0.0002   0.0002   0.0000  
 2003   0.0000   0.0004   0.0017   0.0000   0.0007   0.0000   0.0002  
 2004   0.0001   0.0002   0.0002   0.0006   0.0001   0.0002   0.0001  
 2005   0.0010   0.0006   0.0002   0.0002   0.0002   0.0001   0.0001  
 
 Check Ed = 1: 
 1995   0.0001   0.0000   0.0003   0.0020   0.0019   0.0001   0.0001  
 1996   0.0027   0.0004   0.0000   0.0021   0.0000   0.0000   0.0000  
 1997   0.0000   0.0001   0.0007   0.0045   0.0006   0.0002   0.0000  
 1998   0.0004   0.0001   0.0001   0.0001   0.0000   0.0000   0.0000  
 1999   0.0003   0.0010   0.0005   0.0002   0.0002   0.0003   0.0002  
 2000   0.0011   0.0005   0.0000   0.0026   0.0001   0.0000   0.0001  
 2001   0.0010   0.0002   0.0005   0.0004   0.0004   0.0004   0.0001  
 2002   0.0007   0.0009   0.0003   0.0046   0.0009   0.0000   0.0001  
 2003   0.0001   0.0002   0.0001   0.0023   0.0000   0.0001   0.0001  
 2004   0.0024   0.0003   0.0000   0.0039   0.0006   0.0001   0.0000  
 2005   0.0001   0.0000   0.0004   0.0002   0.0004   0.0009   0.0002  
 
 Min and Max R_sq for Emax:
 0.94    33   4   2   4   2
 1.00    35   3   2   1   3
 
 Min and Max rsk_ratio:
 0.93     3     1     2  
 1.09     6     2     6  
 
 Sectoral Choices:
40.45   3.56   8.54   2.60   3.15   8.61   4.23  28.85  
 
 Log Wages:
 1.00   1.39   2.11   1.32   1.25   1.86   1.73  
 
 Transitions:
80.25   1.86   2.69   0.59   1.62   3.54   0.98   8.46  
18.30  75.12   1.27   0.38   1.10   1.19   0.81   1.83  
14.89   0.51  79.87   0.64   0.63   1.27   0.64   1.54  
13.37   0.28   1.66  81.21   0.56   1.04   0.48   1.42  
19.04   1.15   2.60   1.02  68.98   1.96   1.48   3.77  
16.23   0.47   1.62   0.44   0.93  77.04   0.88   2.39  
13.86   0.45   0.49   0.16   0.50   0.79  82.67   1.08  
10.34   0.23   0.56   0.14   0.45   0.76   0.37  87.15  
 
 Labor Shares 0 (2005):
 0.17   0.31   0.20   0.29   0.29   0.14   0.20  
 
 Labor Shares 1 (2005):
 0.07   0.26   0.36   0.10   0.28   0.24   0.67  
 
 Capital Shares (2005):
 0.76   0.44   0.44   0.60   0.43   0.63   0.13  
 
 *****************
 Loss Function
   69.6866363412582     
 *****************
 
 Iteration Time     :   151.7500    
 
 =====================================================================
 =====================================================================
