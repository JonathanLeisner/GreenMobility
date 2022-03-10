clear
clear matrix

set matsize 800



cd "C:\Users\Dix-Carneiro\Documents\Trade and Labor Dynamics\Revision_Econometrica\CostsOfMobility"


******************************************
* Getting the equilibrium returns to skill
******************************************

insheet using "Outcomes.csv"

rename year    ano
rename rsk_eq1_0 rsk1_0
rename rsk_eq2_0 rsk2_0
rename rsk_eq3_0 rsk3_0
rename rsk_eq4_0 rsk4_0
rename rsk_eq5_0 rsk5_0
rename rsk_eq6_0 rsk6_0
rename rsk_eq7_0 rsk7_0
rename rsk_eq1_1 rsk1_1
rename rsk_eq2_1 rsk2_1
rename rsk_eq3_1 rsk3_1
rename rsk_eq4_1 rsk4_1
rename rsk_eq5_1 rsk5_1
rename rsk_eq6_1 rsk6_1
rename rsk_eq7_1 rsk7_1

keep ano rsk1_* rsk2_* rsk3_* rsk4_* rsk5_* rsk6_* rsk7_*

matrix LOW = J(7,2,0)
matrix UP = J(7,2,0)

matrix LOW[1,1] = 0.8
matrix LOW[2,1] = 0.8
matrix LOW[3,1] = 1.0
matrix LOW[4,1] = 1.0
matrix LOW[5,1] = 0.8
matrix LOW[6,1] = 1.0
matrix LOW[7,1] = 0.8

matrix UP[1,1] = 2.9
matrix UP[2,1] = 3.0
matrix UP[3,1] = 3.8
matrix UP[4,1] = 3.8
matrix UP[5,1] = 3.0
matrix UP[6,1] = 3.8
matrix UP[7,1] = 3.0

matrix LOW[1,2] = 1.4
matrix LOW[2,2] = 1.6
matrix LOW[3,2] = 2.0
matrix LOW[4,2] = 2.0
matrix LOW[5,2] = 1.8
matrix LOW[6,2] = 2.0
matrix LOW[7,2] = 1.8

matrix UP[1,2] = 5.1
matrix UP[2,2] = 5.6
matrix UP[3,2] = 6.8
matrix UP[4,2] = 6.8
matrix UP[5,2] = 5.6
matrix UP[6,2] = 6.8
matrix UP[7,2] = 5.6

gen cut1_0 = LOW[1,1] + (UP[1,1]-LOW[1,1])/2.0
gen cut2_0 = LOW[2,1] + (UP[2,1]-LOW[2,1])/2.0
gen cut3_0 = LOW[3,1] + (UP[3,1]-LOW[3,1])/2.0
gen cut4_0 = LOW[4,1] + (UP[4,1]-LOW[4,1])/2.0
gen cut5_0 = LOW[5,1] + (UP[5,1]-LOW[5,1])/2.0
gen cut6_0 = LOW[6,1] + (UP[6,1]-LOW[6,1])/2.0
gen cut7_0 = LOW[7,1] + (UP[7,1]-LOW[7,1])/2.0

gen cut1_1 = LOW[1,2] + (UP[1,2]-LOW[1,2])/2.0
gen cut2_1 = LOW[2,2] + (UP[2,2]-LOW[2,2])/2.0
gen cut3_1 = LOW[3,2] + (UP[3,2]-LOW[3,2])/2.0
gen cut4_1 = LOW[4,2] + (UP[4,2]-LOW[4,2])/2.0
gen cut5_1 = LOW[5,2] + (UP[5,2]-LOW[5,2])/2.0
gen cut6_1 = LOW[6,2] + (UP[6,2]-LOW[6,2])/2.0
gen cut7_1 = LOW[7,2] + (UP[7,2]-LOW[7,2])/2.0

gen dummy1_0 = (rsk1_0 >= cut1_0)
gen dummy2_0 = (rsk2_0 >= cut2_0)
gen dummy3_0 = (rsk3_0 >= cut3_0)
gen dummy4_0 = (rsk4_0 >= cut4_0)
gen dummy5_0 = (rsk5_0 >= cut5_0)
gen dummy6_0 = (rsk6_0 >= cut6_0)
gen dummy7_0 = (rsk7_0 >= cut7_0)

gen dummy1_1 = (rsk1_1 >= cut1_1)
gen dummy2_1 = (rsk2_1 >= cut2_1)
gen dummy3_1 = (rsk3_1 >= cut3_1)
gen dummy4_1 = (rsk4_1 >= cut4_1)
gen dummy5_1 = (rsk5_1 >= cut5_1)
gen dummy6_1 = (rsk6_1 >= cut6_1)
gen dummy7_1 = (rsk7_1 >= cut7_1)

sort ano

save rsk, replace




********************************************************
* Getting the coefficients for the approximation of EMAX
********************************************************

clear

insheet using "Emax_Coef.csv"

rename lag sector
rename reg n
rename pi_myopic PI

forvalues a = 26(1)60{
	forvalues s = 0(1)7{
		forvalues g = 1(1)2{
			forvalues e = 1(1)4{
				forvalues tp = 1(1)3{
					matrix PI_`a'_`s'_`g'_`e'_`tp' = J(134,1,0)
					preserve
					qui keep if age == `a' & sector == `s' & gender == `g' & educ == `e' & type == `tp'
					qui keep PI
					qui mkmat PI
					matrix PI_`a'_`s'_`g'_`e'_`tp' = PI
					restore
				}
			}
		}
	}
}




************************************************************************************
* We read the optimal parameter vector, in order to construct individual entry costs
************************************************************************************
 
cd "C:\Users\Dix-Carneiro\Documents\Trade and Labor Dynamics\Revision_Econometrica\CostsOfMobility"

clear

insheet using "Starting_Point.csv"

keep v2

mkmat v2, matrix(A)

* Costs of mobility function parameters
matrix costs = A[100..125,1]
* Type Probabilities parameters
matrix aux = A[151..174,1]
* Type Probabilities parameters
matrix type_coef = J(12,2,0)
* Value Residual Sector Parameters
matrix vhome = A[1..7,1]
* Permanent Unobserved Heterogeneity in VResidual
matrix omega = J(2,1,0)
matrix omega[1,1] = A[133,1]
matrix omega[2,1] = A[141,1]
* Permanent Unobserved Heterogeneity in Costs
matrix lambda = J(2,1,0)
matrix lambda[1,1] = A[149,1]
matrix lambda[2,1] = A[150,1]


forvalues i = 1(1)2{
	forvalues j = 1(1)12{
		matrix type_coef[`j',`i'] = aux[(`i'-1)*12+`j',1]
	}
}


clear 

use New_Painel

sort ano

merge ano using rsk

gen dummy_educ2 = (educ_new == 2)
gen dummy_educ3 = (educ_new == 3)
gen dummy_educ4 = (educ_new == 4)

gen age = idade_new - 25
gen age_2 = (idade_new - 25)^2

keep if idade_new >= 25 & idade_new <= 60

reg log_w year95-year05 gender dummy_educ2-dummy_educ4 age age_2 exper1-exper7 if sector ~= 0, noc
predict LogWageData
gen WageData  = exp(LogWageData)*exp(0.5*e(rmse)^2)

reg log_w year95-year05 gender dummy_educ2-dummy_educ4 age age_2 exper1-exper7 if sector == 1, noc
predict LogWage1Data
gen Wage1Data = exp(LogWage1Data)*exp(0.5*e(rmse)^2)

reg log_w year95-year05 gender dummy_educ2-dummy_educ4 age age_2 exper1-exper7 if sector == 2, noc
predict LogWage2Data
gen Wage2Data = exp(LogWage2Data)*exp(0.5*e(rmse)^2)

reg log_w year95-year05 gender dummy_educ2-dummy_educ4 age age_2 exper1-exper7 if sector == 3, noc
predict LogWage3Data
gen Wage3Data = exp(LogWage3Data)*exp(0.5*e(rmse)^2)

reg log_w year95-year05 gender dummy_educ2-dummy_educ4 age age_2 exper1-exper7 if sector == 4, noc
predict LogWage4Data
gen Wage4Data = exp(LogWage4Data)*exp(0.5*e(rmse)^2)

reg log_w year95-year05 gender dummy_educ2-dummy_educ4 age age_2 exper1-exper7 if sector == 5, noc
predict LogWage5Data
gen Wage5Data = exp(LogWage5Data)*exp(0.5*e(rmse)^2)

reg log_w year95-year05 gender dummy_educ2-dummy_educ4 age age_2 exper1-exper7 if sector == 6, noc
predict LogWage6Data
gen Wage6Data = exp(LogWage6Data)*exp(0.5*e(rmse)^2)

reg log_w year95-year05 gender dummy_educ2-dummy_educ4 age age_2 exper1-exper7 if sector == 7, noc
predict LogWage7Data
gen Wage7Data = exp(LogWage7Data)*exp(0.5*e(rmse)^2)



gen Wage0Data_tp1 = exp(vhome[1,1] + vhome[2,1]*gender + vhome[3,1]*dummy_educ2 + ///
                        vhome[4,1]*dummy_educ3 + vhome[5,1]*dummy_educ4 + ///
                        vhome[6,1]*age + vhome[7,1]*age^2)					  
gen Wage0Data_tp2 = exp(omega[1,1] + vhome[1,1] + vhome[2,1]*gender + vhome[3,1]*dummy_educ2 + ///
                        vhome[4,1]*dummy_educ3 + vhome[5,1]*dummy_educ4 + ///
                        vhome[6,1]*age + vhome[7,1]*age^2)
gen Wage0Data_tp3 = exp(omega[2,1] + vhome[1,1] + vhome[2,1]*gender + vhome[3,1]*dummy_educ2 + ///
                        vhome[4,1]*dummy_educ3 + vhome[5,1]*dummy_educ4 + ///
                        vhome[6,1]*age + vhome[7,1]*age^2)				
						
gen Wage1Data_tp1 =	Wage1Data				
gen Wage1Data_tp2 = Wage1Data
gen Wage1Data_tp3 = Wage1Data

gen Wage2Data_tp1 = Wage2Data
gen Wage2Data_tp2 = Wage2Data
gen Wage2Data_tp3 = Wage2Data

gen Wage3Data_tp1 = Wage3Data
gen Wage3Data_tp2 = Wage3Data
gen Wage3Data_tp3 = Wage3Data

gen Wage4Data_tp1 = Wage4Data
gen Wage4Data_tp2 = Wage4Data
gen Wage4Data_tp3 = Wage4Data
						
gen Wage5Data_tp1 = Wage5Data
gen Wage5Data_tp2 = Wage5Data
gen Wage5Data_tp3 = Wage5Data

gen Wage6Data_tp1 = Wage6Data
gen Wage6Data_tp2 = Wage6Data
gen Wage6Data_tp3 = Wage6Data

gen Wage7Data_tp1 = Wage7Data
gen Wage7Data_tp2 = Wage7Data
gen Wage7Data_tp3 = Wage7Data
						
gen PV1 = 0
gen PV2 = 0
gen PV3 = 0						  
						  
gen exper1_tom = exper1
gen exper2_tom = exper2
gen exper3_tom = exper3
gen exper4_tom = exper4
gen exper5_tom = exper5
gen exper6_tom = exper6
gen exper7_tom = exper7

replace exper1_tom = max(exper1 + 1,9) if sector == 1
replace exper2_tom = max(exper2 + 1,9) if sector == 2
replace exper3_tom = max(exper3 + 1,9) if sector == 3
replace exper4_tom = max(exper4 + 1,9) if sector == 4
replace exper5_tom = max(exper5 + 1,9) if sector == 5
replace exper6_tom = max(exper6 + 1,9) if sector == 6
replace exper7_tom = max(exper7 + 1,9) if sector == 7

gen rsk1 = 0
gen rsk2 = 0
gen rsk3 = 0
gen rsk4 = 0
gen rsk5 = 0
gen rsk6 = 0
gen rsk7 = 0

gen cut1 = 0
gen cut2 = 0
gen cut3 = 0
gen cut4 = 0
gen cut5 = 0
gen cut6 = 0
gen cut7 = 0

gen dummy1 = 0
gen dummy2 = 0
gen dummy3 = 0
gen dummy4 = 0
gen dummy5 = 0
gen dummy6 = 0
gen dummy7 = 0

replace rsk1 = rsk1_0 if educ_new == 1 | educ_new == 2
replace rsk2 = rsk2_0 if educ_new == 1 | educ_new == 2
replace rsk3 = rsk3_0 if educ_new == 1 | educ_new == 2
replace rsk4 = rsk4_0 if educ_new == 1 | educ_new == 2
replace rsk5 = rsk5_0 if educ_new == 1 | educ_new == 2
replace rsk6 = rsk6_0 if educ_new == 1 | educ_new == 2
replace rsk7 = rsk7_0 if educ_new == 1 | educ_new == 2

replace rsk1 = rsk1_1 if educ_new == 3 | educ_new == 4
replace rsk2 = rsk2_1 if educ_new == 3 | educ_new == 4
replace rsk3 = rsk3_1 if educ_new == 3 | educ_new == 4
replace rsk4 = rsk4_1 if educ_new == 3 | educ_new == 4
replace rsk5 = rsk5_1 if educ_new == 3 | educ_new == 4
replace rsk6 = rsk6_1 if educ_new == 3 | educ_new == 4
replace rsk7 = rsk7_1 if educ_new == 3 | educ_new == 4

replace dummy1 = dummy1_0 if educ_new == 1 | educ_new == 2
replace dummy2 = dummy2_0 if educ_new == 1 | educ_new == 2
replace dummy3 = dummy3_0 if educ_new == 1 | educ_new == 2
replace dummy4 = dummy4_0 if educ_new == 1 | educ_new == 2
replace dummy5 = dummy5_0 if educ_new == 1 | educ_new == 2
replace dummy6 = dummy6_0 if educ_new == 1 | educ_new == 2
replace dummy7 = dummy7_0 if educ_new == 1 | educ_new == 2

replace dummy1 = dummy1_1 if educ_new == 3 | educ_new == 4
replace dummy2 = dummy2_1 if educ_new == 3 | educ_new == 4
replace dummy3 = dummy3_1 if educ_new == 3 | educ_new == 4
replace dummy4 = dummy4_1 if educ_new == 3 | educ_new == 4
replace dummy5 = dummy5_1 if educ_new == 3 | educ_new == 4
replace dummy6 = dummy6_1 if educ_new == 3 | educ_new == 4
replace dummy7 = dummy7_1 if educ_new == 3 | educ_new == 4

replace cut1 = cut1_0 if educ_new == 1 | educ_new == 2
replace cut2 = cut2_0 if educ_new == 1 | educ_new == 2
replace cut3 = cut3_0 if educ_new == 1 | educ_new == 2
replace cut4 = cut4_0 if educ_new == 1 | educ_new == 2
replace cut5 = cut5_0 if educ_new == 1 | educ_new == 2
replace cut6 = cut6_0 if educ_new == 1 | educ_new == 2
replace cut7 = cut7_0 if educ_new == 1 | educ_new == 2

replace cut1 = cut1_1 if educ_new == 3 | educ_new == 4
replace cut2 = cut2_1 if educ_new == 3 | educ_new == 4
replace cut3 = cut3_1 if educ_new == 3 | educ_new == 4
replace cut4 = cut4_1 if educ_new == 3 | educ_new == 4
replace cut5 = cut5_1 if educ_new == 3 | educ_new == 4
replace cut6 = cut6_1 if educ_new == 3 | educ_new == 4
replace cut7 = cut7_1 if educ_new == 3 | educ_new == 4

drop rsk1_* rsk2_* rsk3_* rsk4_* rsk5_* rsk6_* rsk7_* cut1_* cut2_* cut3_* cut4_* cut5_* cut6_* cut7_* ///
     dummy1_* dummy2_* dummy3_* dummy4_* dummy5_* dummy6_* dummy7_*
	 

forvalues a = 26(1)60{
	forvalues s = 0(1)7{
		forvalues g = 1(1)2{
			forvalues e = 1(1)4{
				forvalues tp = 1(1)3{
					qui replace PV`tp' = Wage`s'Data_tp`tp' + 0.95*(PI_`a'_`s'_`g'_`e'_`tp'[1,1]*1.0 + ///
						                                    PI_`a'_`s'_`g'_`e'_`tp'[2,1]*rsk1 + ///
															PI_`a'_`s'_`g'_`e'_`tp'[3,1]*rsk2 + ///
															PI_`a'_`s'_`g'_`e'_`tp'[4,1]*rsk3 + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[5,1]*rsk4 + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[6,1]*rsk5 + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[7,1]*rsk6 + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[8,1]*rsk7 + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[9,1]*dummy1*(rsk1-cut1) + ///
															PI_`a'_`s'_`g'_`e'_`tp'[10,1]*dummy2*(rsk2-cut2) + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[11,1]*dummy3*(rsk3-cut3) + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[12,1]*dummy4*(rsk4-cut4) + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[13,1]*dummy5*(rsk5-cut5) + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[14,1]*dummy6*(rsk6-cut6) + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[15,1]*dummy7*(rsk7-cut7) + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[16,1]*exper1_tom + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[17,1]*exper2_tom + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[18,1]*exper3_tom + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[19,1]*exper4_tom + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[20,1]*exper5_tom + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[21,1]*exper6_tom + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[22,1]*exper7_tom + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[23,1]*rsk1^2 + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[24,1]*rsk2^2 + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[25,1]*rsk3^2 + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[26,1]*rsk4^2 + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[27,1]*rsk5^2 + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[28,1]*rsk6^2 + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[29,1]*rsk7^2 + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[30,1]*dummy1*(rsk1-cut1)^2 + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[31,1]*dummy2*(rsk2-cut2)^2 + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[32,1]*dummy3*(rsk3-cut3)^2 + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[33,1]*dummy4*(rsk4-cut4)^2 + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[34,1]*dummy5*(rsk5-cut5)^2 + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[35,1]*dummy6*(rsk6-cut6)^2 + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[36,1]*dummy7*(rsk7-cut7)^2 + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[37,1]*exper1_tom^2 + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[38,1]*exper2_tom^2 + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[39,1]*exper3_tom^2 + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[40,1]*exper4_tom^2 + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[41,1]*exper5_tom^2 + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[42,1]*exper6_tom^2 + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[43,1]*exper7_tom^2 + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[44,1]*rsk1*rsk2 + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[45,1]*rsk1*rsk3 + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[46,1]*rsk1*rsk4 + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[47,1]*rsk1*rsk5 + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[48,1]*rsk1*rsk6 + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[49,1]*rsk1*rsk7 + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[50,1]*rsk2*rsk3 + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[51,1]*rsk2*rsk4 + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[52,1]*rsk2*rsk5 + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[53,1]*rsk2*rsk6 + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[54,1]*rsk2*rsk7 + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[55,1]*rsk3*rsk4 + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[56,1]*rsk3*rsk5 + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[57,1]*rsk3*rsk6 + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[58,1]*rsk3*rsk7 + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[59,1]*rsk4*rsk5 + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[60,1]*rsk4*rsk6 + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[61,1]*rsk4*rsk7 + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[62,1]*rsk5*rsk6 + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[63,1]*rsk5*rsk7 + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[64,1]*rsk6*rsk7 + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[65,1]*rsk1*exper1_tom + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[66,1]*rsk1*exper2_tom + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[67,1]*rsk1*exper3_tom + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[68,1]*rsk1*exper4_tom + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[69,1]*rsk1*exper5_tom + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[70,1]*rsk1*exper6_tom + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[71,1]*rsk1*exper7_tom + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[72,1]*rsk2*exper1_tom + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[73,1]*rsk2*exper2_tom + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[74,1]*rsk2*exper3_tom + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[75,1]*rsk2*exper4_tom + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[76,1]*rsk2*exper5_tom + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[77,1]*rsk2*exper6_tom + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[78,1]*rsk2*exper7_tom + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[79,1]*rsk3*exper1_tom + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[80,1]*rsk3*exper2_tom + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[81,1]*rsk3*exper3_tom + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[82,1]*rsk3*exper4_tom + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[83,1]*rsk3*exper5_tom + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[84,1]*rsk3*exper6_tom + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[85,1]*rsk3*exper7_tom + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[86,1]*rsk4*exper1_tom + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[87,1]*rsk4*exper2_tom + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[88,1]*rsk4*exper3_tom + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[89,1]*rsk4*exper4_tom + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[90,1]*rsk4*exper5_tom + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[91,1]*rsk4*exper6_tom + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[92,1]*rsk4*exper7_tom + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[93,1]*rsk5*exper1_tom + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[94,1]*rsk5*exper2_tom + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[95,1]*rsk5*exper3_tom + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[96,1]*rsk5*exper4_tom + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[97,1]*rsk5*exper5_tom + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[98,1]*rsk5*exper6_tom + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[99,1]*rsk5*exper7_tom + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[100,1]*rsk6*exper1_tom + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[101,1]*rsk6*exper2_tom + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[102,1]*rsk6*exper3_tom + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[103,1]*rsk6*exper4_tom + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[104,1]*rsk6*exper5_tom + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[105,1]*rsk6*exper6_tom + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[106,1]*rsk6*exper7_tom + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[107,1]*rsk7*exper1_tom + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[108,1]*rsk7*exper2_tom + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[109,1]*rsk7*exper3_tom + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[110,1]*rsk7*exper4_tom + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[111,1]*rsk7*exper5_tom + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[112,1]*rsk7*exper6_tom + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[113,1]*rsk7*exper7_tom + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[114,1]*exper1_tom*exper2_tom + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[115,1]*exper1_tom*exper3_tom + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[116,1]*exper1_tom*exper4_tom + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[117,1]*exper1_tom*exper5_tom + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[118,1]*exper1_tom*exper6_tom + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[119,1]*exper1_tom*exper7_tom + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[120,1]*exper2_tom*exper3_tom + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[121,1]*exper2_tom*exper4_tom + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[122,1]*exper2_tom*exper5_tom + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[123,1]*exper2_tom*exper6_tom + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[124,1]*exper2_tom*exper7_tom + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[125,1]*exper3_tom*exper4_tom + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[126,1]*exper3_tom*exper5_tom + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[127,1]*exper3_tom*exper6_tom + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[128,1]*exper3_tom*exper7_tom + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[129,1]*exper4_tom*exper5_tom + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[130,1]*exper4_tom*exper6_tom + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[131,1]*exper4_tom*exper7_tom + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[132,1]*exper5_tom*exper6_tom + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[133,1]*exper5_tom*exper7_tom + /// 
															PI_`a'_`s'_`g'_`e'_`tp'[134,1]*exper6_tom*exper7_tom) ///
															if idade_new == `a'-1 & sector == `s' & gender == `g'-1 & educ == `e'
					qui replace PV`tp' = Wage`s'Data_tp`tp' if idade_new == 60			
				}										 
			}
		}
	}										
}
					

					
drop Wage0Data* Wage1Data* Wage2Data* Wage3Data* Wage4Data* Wage5Data* Wage6Data* Wage7Data*


gen term2 = exp(type_coef[1,1] + type_coef[2,1]*exper1 + type_coef[3,1]*exper2 + ///
                type_coef[4,1]*exper3 + type_coef[5,1]*exper4 + ///
				type_coef[6,1]*exper5 + type_coef[7,1]*exper6 + type_coef[8,1]*exper7 + ///
				type_coef[9,1]*gender + type_coef[10,1]*dummy_educ2 + type_coef[11,1]*dummy_educ3 + ///
				type_coef[12,1]*dummy_educ4)
gen term3 = exp(type_coef[1,2] + type_coef[2,2]*exper1 + type_coef[3,2]*exper2 + ///
                type_coef[4,2]*exper3 + type_coef[5,2]*exper4 + ///
				type_coef[6,2]*exper5 + type_coef[7,2]*exper6 + type_coef[8,2]*exper7 + ///
				type_coef[9,2]*gender + type_coef[10,2]*dummy_educ2 + type_coef[11,2]*dummy_educ3 + ///
				type_coef[12,2]*dummy_educ4)		
				
gen prob2 = term2  / (1.0 + term2 + term3)
gen prob3 = term3  / (1.0 + term2 + term3)
gen prob1 = 1 - prob2 - prob3

gen unif = runiform()

gen tp = .
replace tp = 1 if unif < prob1
replace tp = 2 if unif >= prob1 & unif < prob1 + prob2
replace tp = 3 if unif >= prob1 + prob2

gen PV = prob1*PV1 + prob2*PV2 + prob3*PV3

save DataSetWithPV, replace

**************************************************************************************************
**************************************************************************************************

clear

cd "C:\Users\Dix-Carneiro\Documents\Trade and Labor Dynamics\Revision_Econometrica\CostsOfMobility"

use DataSetWithPV

reg log_w gender dummy_educ2-dummy_educ4 age age_2 exper1-exper7 if sector ~= 0
predict LogWageHat
gen WageHat = exp(LogWageHat)*exp(0.5*e(rmse)^2)
matrix wage_coef = e(b)
matrix wage_rmse = e(rmse)

gen lambda = .
replace lambda = 0 if tp == 1
replace lambda = lambda[1,1] if tp == 2
replace lambda = lambda[2,1] if tp == 3


gen cost_entry1 = exp(costs[1,1]+costs[21,1]*gender+costs[22,1]*dummy_educ2+costs[23,1]*dummy_educ3+costs[24,1]*dummy_educ4+costs[25,1]*age+costs[26,1]*age_2 + lambda)  ///
if sector == 0
gen cost_entry2 = exp(costs[2,1]+costs[21,1]*gender+costs[22,1]*dummy_educ2+costs[23,1]*dummy_educ3+costs[24,1]*dummy_educ4+costs[25,1]*age+costs[26,1]*age_2 + lambda)  ///
if sector == 0
gen cost_entry3 = exp(costs[3,1]+costs[21,1]*gender+costs[22,1]*dummy_educ2+costs[23,1]*dummy_educ3+costs[24,1]*dummy_educ4+costs[25,1]*age+costs[26,1]*age_2 + lambda)  ///
if sector == 0
gen cost_entry4 = exp(costs[4,1]+costs[21,1]*gender+costs[22,1]*dummy_educ2+costs[23,1]*dummy_educ3+costs[24,1]*dummy_educ4+costs[25,1]*age+costs[26,1]*age_2 + lambda)  ///
if sector == 0
gen cost_entry5 = exp(costs[5,1]+costs[21,1]*gender+costs[22,1]*dummy_educ2+costs[23,1]*dummy_educ3+costs[24,1]*dummy_educ4+costs[25,1]*age+costs[26,1]*age_2 + lambda)  ///
if sector == 0
gen cost_entry6 = exp(costs[6,1]+costs[21,1]*gender+costs[22,1]*dummy_educ2+costs[23,1]*dummy_educ3+costs[24,1]*dummy_educ4+costs[25,1]*age+costs[26,1]*age_2 + lambda)  ///
if sector == 0
gen cost_entry7 = exp(costs[7,1]+costs[21,1]*gender+costs[22,1]*dummy_educ2+costs[23,1]*dummy_educ3+costs[24,1]*dummy_educ4+costs[25,1]*age+costs[26,1]*age_2 + lambda)  ///
if sector == 0


replace cost_entry1 = . ///
if sector == 1
replace cost_entry2 = exp(costs[13+2,1]+costs[21,1]*gender+costs[22,1]*dummy_educ2+costs[23,1]*dummy_educ3+costs[24,1]*dummy_educ4+costs[25,1]*age+costs[26,1]*age_2 + lambda)  ///
if sector == 1
replace cost_entry3 = exp(costs[13+3,1]+costs[21,1]*gender+costs[22,1]*dummy_educ2+costs[23,1]*dummy_educ3+costs[24,1]*dummy_educ4+costs[25,1]*age+costs[26,1]*age_2 + lambda)  ///
if sector == 1
replace cost_entry4 = exp(costs[13+4,1]+costs[21,1]*gender+costs[22,1]*dummy_educ2+costs[23,1]*dummy_educ3+costs[24,1]*dummy_educ4+costs[25,1]*age+costs[26,1]*age_2 + lambda)  ///
if sector == 1
replace cost_entry5 = exp(costs[13+5,1]+costs[21,1]*gender+costs[22,1]*dummy_educ2+costs[23,1]*dummy_educ3+costs[24,1]*dummy_educ4+costs[25,1]*age+costs[26,1]*age_2 + lambda)  ///
if sector == 1
replace cost_entry6 = exp(costs[13+6,1]+costs[21,1]*gender+costs[22,1]*dummy_educ2+costs[23,1]*dummy_educ3+costs[24,1]*dummy_educ4+costs[25,1]*age+costs[26,1]*age_2 + lambda)  ///
if sector == 1
replace cost_entry7 = exp(costs[13+7,1]+costs[21,1]*gender+costs[22,1]*dummy_educ2+costs[23,1]*dummy_educ3+costs[24,1]*dummy_educ4+costs[25,1]*age+costs[26,1]*age_2 + lambda)  ///
if sector == 1


replace cost_entry1 = exp(costs[6+2,1]+costs[13+1,1]+costs[21,1]*gender+costs[22,1]*dummy_educ2+costs[23,1]*dummy_educ3+costs[24,1]*dummy_educ4+costs[25,1]*age+costs[26,1]*age_2 + lambda)  ///
if sector == 2
replace cost_entry2 = . ///
if sector == 2
replace cost_entry3 = exp(costs[6+2,1]+costs[13+3,1]+costs[21,1]*gender+costs[22,1]*dummy_educ2+costs[23,1]*dummy_educ3+costs[24,1]*dummy_educ4+costs[25,1]*age+costs[26,1]*age_2 + lambda)  ///
if sector == 2
replace cost_entry4 = exp(costs[6+2,1]+costs[13+4,1]+costs[21,1]*gender+costs[22,1]*dummy_educ2+costs[23,1]*dummy_educ3+costs[24,1]*dummy_educ4+costs[25,1]*age+costs[26,1]*age_2 + lambda)  ///
if sector == 2
replace cost_entry5 = exp(costs[6+2,1]+costs[13+5,1]+costs[21,1]*gender+costs[22,1]*dummy_educ2+costs[23,1]*dummy_educ3+costs[24,1]*dummy_educ4+costs[25,1]*age+costs[26,1]*age_2 + lambda)  ///
if sector == 2
replace cost_entry6 = exp(costs[6+2,1]+costs[13+6,1]+costs[21,1]*gender+costs[22,1]*dummy_educ2+costs[23,1]*dummy_educ3+costs[24,1]*dummy_educ4+costs[25,1]*age+costs[26,1]*age_2 + lambda)  ///
if sector == 2
replace cost_entry7 = exp(costs[6+2,1]+costs[13+7,1]+costs[21,1]*gender+costs[22,1]*dummy_educ2+costs[23,1]*dummy_educ3+costs[24,1]*dummy_educ4+costs[25,1]*age+costs[26,1]*age_2 + lambda)  ///
if sector == 2


replace cost_entry1 = exp(costs[6+3,1]+costs[13+1,1]+costs[21,1]*gender+costs[22,1]*dummy_educ2+costs[23,1]*dummy_educ3+costs[24,1]*dummy_educ4+costs[25,1]*age+costs[26,1]*age_2 + lambda)  ///
if sector == 3
replace cost_entry2 = exp(costs[6+3,1]+costs[13+2,1]+costs[21,1]*gender+costs[22,1]*dummy_educ2+costs[23,1]*dummy_educ3+costs[24,1]*dummy_educ4+costs[25,1]*age+costs[26,1]*age_2 + lambda)  ///
if sector == 3
replace cost_entry3 = . ///
if sector == 3
replace cost_entry4 = exp(costs[6+3,1]+costs[13+4,1]+costs[21,1]*gender+costs[22,1]*dummy_educ2+costs[23,1]*dummy_educ3+costs[24,1]*dummy_educ4+costs[25,1]*age+costs[26,1]*age_2 + lambda)  ///
if sector == 3
replace cost_entry5 = exp(costs[6+3,1]+costs[13+5,1]+costs[21,1]*gender+costs[22,1]*dummy_educ2+costs[23,1]*dummy_educ3+costs[24,1]*dummy_educ4+costs[25,1]*age+costs[26,1]*age_2 + lambda)  ///
if sector == 3
replace cost_entry6 = exp(costs[6+3,1]+costs[13+6,1]+costs[21,1]*gender+costs[22,1]*dummy_educ2+costs[23,1]*dummy_educ3+costs[24,1]*dummy_educ4+costs[25,1]*age+costs[26,1]*age_2 + lambda)  ///
if sector == 3
replace cost_entry7 = exp(costs[6+3,1]+costs[13+7,1]+costs[21,1]*gender+costs[22,1]*dummy_educ2+costs[23,1]*dummy_educ3+costs[24,1]*dummy_educ4+costs[25,1]*age+costs[26,1]*age_2 + lambda)  ///
if sector == 3


replace cost_entry1 = exp(costs[6+4,1]+costs[13+1,1]+costs[21,1]*gender+costs[22,1]*dummy_educ2+costs[23,1]*dummy_educ3+costs[24,1]*dummy_educ4+costs[25,1]*age+costs[26,1]*age_2 + lambda)  ///
if sector == 4
replace cost_entry2 = exp(costs[6+4,1]+costs[13+2,1]+costs[21,1]*gender+costs[22,1]*dummy_educ2+costs[23,1]*dummy_educ3+costs[24,1]*dummy_educ4+costs[25,1]*age+costs[26,1]*age_2 + lambda)  ///
if sector == 4
replace cost_entry3 = exp(costs[6+4,1]+costs[13+3,1]+costs[21,1]*gender+costs[22,1]*dummy_educ2+costs[23,1]*dummy_educ3+costs[24,1]*dummy_educ4+costs[25,1]*age+costs[26,1]*age_2 + lambda) ///
if sector == 4
replace cost_entry4 = . ///
if sector == 4
replace cost_entry5 = exp(costs[6+4,1]+costs[13+5,1]+costs[21,1]*gender+costs[22,1]*dummy_educ2+costs[23,1]*dummy_educ3+costs[24,1]*dummy_educ4+costs[25,1]*age+costs[26,1]*age_2 + lambda) ///
if sector == 4
replace cost_entry6 = exp(costs[6+4,1]+costs[13+6,1]+costs[21,1]*gender+costs[22,1]*dummy_educ2+costs[23,1]*dummy_educ3+costs[24,1]*dummy_educ4+costs[25,1]*age+costs[26,1]*age_2 + lambda) ///
if sector == 4
replace cost_entry7 = exp(costs[6+4,1]+costs[13+7,1]+costs[21,1]*gender+costs[22,1]*dummy_educ2+costs[23,1]*dummy_educ3+costs[24,1]*dummy_educ4+costs[25,1]*age+costs[26,1]*age_2 + lambda)  ///
if sector == 4


replace cost_entry1 = exp(costs[6+5,1]+costs[13+1,1]+costs[21,1]*gender+costs[22,1]*dummy_educ2+costs[23,1]*dummy_educ3+costs[24,1]*dummy_educ4+costs[25,1]*age+costs[26,1]*age_2 + lambda)  ///
if sector == 5
replace cost_entry2 = exp(costs[6+5,1]+costs[13+2,1]+costs[21,1]*gender+costs[22,1]*dummy_educ2+costs[23,1]*dummy_educ3+costs[24,1]*dummy_educ4+costs[25,1]*age+costs[26,1]*age_2 + lambda)  ///
if sector == 5
replace cost_entry3 = exp(costs[6+5,1]+costs[13+3,1]+costs[21,1]*gender+costs[22,1]*dummy_educ2+costs[23,1]*dummy_educ3+costs[24,1]*dummy_educ4+costs[25,1]*age+costs[26,1]*age_2 + lambda)  ///
if sector == 5
replace cost_entry4 = exp(costs[6+5,1]+costs[13+4,1]+costs[21,1]*gender+costs[22,1]*dummy_educ2+costs[23,1]*dummy_educ3+costs[24,1]*dummy_educ4+costs[25,1]*age+costs[26,1]*age_2 + lambda)  ///
if sector == 5
replace cost_entry5 = . ///
if sector == 5
replace cost_entry6 = exp(costs[6+5,1]+costs[13+6,1]+costs[21,1]*gender+costs[22,1]*dummy_educ2+costs[23,1]*dummy_educ3+costs[24,1]*dummy_educ4+costs[25,1]*age+costs[26,1]*age_2 + lambda)  ///
if sector == 5
replace cost_entry7 = exp(costs[6+5,1]+costs[13+7,1]+costs[21,1]*gender+costs[22,1]*dummy_educ2+costs[23,1]*dummy_educ3+costs[24,1]*dummy_educ4+costs[25,1]*age+costs[26,1]*age_2 + lambda)  ///
if sector == 5


replace cost_entry1 = exp(costs[6+6,1]+costs[13+1,1]+costs[21,1]*gender+costs[22,1]*dummy_educ2+costs[23,1]*dummy_educ3+costs[24,1]*dummy_educ4+costs[25,1]*age+costs[26,1]*age_2 + lambda)  ///
if sector == 6
replace cost_entry2 = exp(costs[6+6,1]+costs[13+2,1]+costs[21,1]*gender+costs[22,1]*dummy_educ2+costs[23,1]*dummy_educ3+costs[24,1]*dummy_educ4+costs[25,1]*age+costs[26,1]*age_2 + lambda)  ///
if sector == 6
replace cost_entry3 = exp(costs[6+6,1]+costs[13+3,1]+costs[21,1]*gender+costs[22,1]*dummy_educ2+costs[23,1]*dummy_educ3+costs[24,1]*dummy_educ4+costs[25,1]*age+costs[26,1]*age_2 + lambda)  ///
if sector == 6
replace cost_entry4 = exp(costs[6+6,1]+costs[13+4,1]+costs[21,1]*gender+costs[22,1]*dummy_educ2+costs[23,1]*dummy_educ3+costs[24,1]*dummy_educ4+costs[25,1]*age+costs[26,1]*age_2 + lambda)  ///
if sector == 6
replace cost_entry5 = exp(costs[6+6,1]+costs[13+5,1]+costs[21,1]*gender+costs[22,1]*dummy_educ2+costs[23,1]*dummy_educ3+costs[24,1]*dummy_educ4+costs[25,1]*age+costs[26,1]*age_2 + lambda)  ///
if sector == 6
replace cost_entry6 = . ///
if sector == 6
replace cost_entry7 = exp(costs[6+6,1]+costs[13+7,1]+costs[21,1]*gender+costs[22,1]*dummy_educ2+costs[23,1]*dummy_educ3+costs[24,1]*dummy_educ4+costs[25,1]*age+costs[26,1]*age_2 + lambda)  ///
if sector == 6


replace cost_entry1 = exp(costs[6+7,1]+costs[13+1,1]+costs[21,1]*gender+costs[22,1]*dummy_educ2+costs[23,1]*dummy_educ3+costs[24,1]*dummy_educ4+costs[25,1]*age+costs[26,1]*age_2 + lambda)  ///
if sector == 7
replace cost_entry2 = exp(costs[6+7,1]+costs[13+2,1]+costs[21,1]*gender+costs[22,1]*dummy_educ2+costs[23,1]*dummy_educ3+costs[24,1]*dummy_educ4+costs[25,1]*age+costs[26,1]*age_2 + lambda)  ///
if sector == 7
replace cost_entry3 = exp(costs[6+7,1]+costs[13+3,1]+costs[21,1]*gender+costs[22,1]*dummy_educ2+costs[23,1]*dummy_educ3+costs[24,1]*dummy_educ4+costs[25,1]*age+costs[26,1]*age_2 + lambda)  ///
if sector == 7
replace cost_entry4 = exp(costs[6+7,1]+costs[13+4,1]+costs[21,1]*gender+costs[22,1]*dummy_educ2+costs[23,1]*dummy_educ3+costs[24,1]*dummy_educ4+costs[25,1]*age+costs[26,1]*age_2 + lambda)  ///
if sector == 7
replace cost_entry5 = exp(costs[6+7,1]+costs[13+5,1]+costs[21,1]*gender+costs[22,1]*dummy_educ2+costs[23,1]*dummy_educ3+costs[24,1]*dummy_educ4+costs[25,1]*age+costs[26,1]*age_2 + lambda)  ///
if sector == 7
replace cost_entry6 = exp(costs[6+7,1]+costs[13+6,1]+costs[21,1]*gender+costs[22,1]*dummy_educ2+costs[23,1]*dummy_educ3+costs[24,1]*dummy_educ4+costs[25,1]*age+costs[26,1]*age_2 + lambda)  ///
if sector == 7
replace cost_entry7 = . ///
if sector == 7

**************************************************
**************************************************

gen Wcost_entry1 = cost_entry1 / WageHat
gen Wcost_entry2 = cost_entry2 / WageHat
gen Wcost_entry3 = cost_entry3 / WageHat
gen Wcost_entry4 = cost_entry4 / WageHat
gen Wcost_entry5 = cost_entry5 / WageHat
gen Wcost_entry6 = cost_entry6 / WageHat
gen Wcost_entry7 = cost_entry7 / WageHat

**************************************************
**************************************************

matrix Wcost_entry_formal = J(7,3,0)

forvalues i = 1(1)7{
	sum Wcost_entry`i' if sector ~= `i' & sector ~= 0, detail
	matrix Wcost_entry_formal[`i',1] = r(p10)
	matrix Wcost_entry_formal[`i',2] = r(p50)
	matrix Wcost_entry_formal[`i',3] = r(mean)
}	

**************************************************
**************************************************

matrix Wcost_entry_residual = J(7,3,0)

forvalues i = 1(1)7{
	sum Wcost_entry`i' if sector ~= `i' & sector == 0, detail
	matrix Wcost_entry_residual[`i',1] = r(p10)
	matrix Wcost_entry_residual[`i',2] = r(p50)
	matrix Wcost_entry_residual[`i',3] = r(mean)
}	

**************************************************
**************************************************

gen Vcost_entry1 = cost_entry1 / PV
gen Vcost_entry2 = cost_entry2 / PV
gen Vcost_entry3 = cost_entry3 / PV
gen Vcost_entry4 = cost_entry4 / PV
gen Vcost_entry5 = cost_entry5 / PV
gen Vcost_entry6 = cost_entry6 / PV
gen Vcost_entry7 = cost_entry7 / PV

**************************************************
**************************************************

matrix Vcost_entry_formal = J(7,3,0)

forvalues i = 1(1)7{
	sum Vcost_entry`i' if sector ~= `i' & sector ~= 0, detail
	matrix Vcost_entry_formal[`i',1] = r(p10)
	matrix Vcost_entry_formal[`i',2] = r(p50)
	matrix Vcost_entry_formal[`i',3] = r(mean)
}	

**************************************************
**************************************************

matrix Vcost_entry_residual = J(7,3,0)

forvalues i = 1(1)7{
	sum Vcost_entry`i' if sector ~= `i' & sector == 0, detail
	matrix Vcost_entry_residual[`i',1] = r(p10)
	matrix Vcost_entry_residual[`i',2] = r(p50)
	matrix Vcost_entry_residual[`i',3] = r(mean)
}	

**************************************************
**************************************************


preserve

cd "C:\Users\Dix-Carneiro\Documents\Trade and Labor Dynamics\Revision_Econometrica\CostsOfMobility\TablesGraphs"

clear
svmat Wcost_entry_formal

rename Wcost_entry_formal1 p10
rename Wcost_entry_formal2 p50
rename Wcost_entry_formal3 mean

outsheet _all using Wcost_entry_formal.csv, comma replace

clear
svmat Wcost_entry_residual

rename Wcost_entry_residual1 p10
rename Wcost_entry_residual2 p50
rename Wcost_entry_residual3 mean

outsheet _all using Wcost_entry_residual.csv, comma replace

clear
svmat Vcost_entry_formal

rename Vcost_entry_formal1 p10
rename Vcost_entry_formal2 p50
rename Vcost_entry_formal3 mean

outsheet _all using Vcost_entry_formal.csv, comma replace

clear
svmat Vcost_entry_residual

rename Vcost_entry_residual1 p10
rename Vcost_entry_residual2 p50
rename Vcost_entry_residual3 mean

outsheet _all using Vcost_entry_residual.csv, comma replace


restore


****************************************************************************************************************************************************************************************************
****************************************************************************************************************************************************************************************************
cd "C:\Users\Dix-Carneiro\Documents\Trade and Labor Dynamics\Revision_Econometrica\CostsOfMobility\TablesGraphs"

kdensity Wcost_entry1 if sector ~= 1 & sector ~= 0 & Wcost_entry1 <= 6, bwidth(0.2) n(50) lcolor(red) ytitle("") xtitle("Cost of Entry Agriculture/Mining From Formal Sector") xlabel(0 1 2 3 4 5 6)
graph save density_cost_entry1_formal, replace
kdensity Wcost_entry2 if sector ~= 2 & sector ~= 0 & Wcost_entry2 <= 6, bwidth(0.2) n(50) lcolor(red) ytitle("") xtitle("Cost of Entry Low Tech Manuf From Formal Sector") xlabel(0 1 2 3 4 5 6)
graph save density_cost_entry2_formal, replace
kdensity Wcost_entry3 if sector ~= 3 & sector ~= 0 & Wcost_entry3 <= 6, bwidth(0.2) n(50) lcolor(red) ytitle("") xtitle("Cost of Entry High tech Manuf From Formal Sector") xlabel(0 1 2 3 4 5 6)
graph save density_cost_entry3_formal, replace
kdensity Wcost_entry4 if sector ~= 4 & sector ~= 0 & Wcost_entry4 <= 6, bwidth(0.2) n(50) lcolor(red) ytitle("") xtitle("Cost of Entry Construction From Formal Sector") xlabel(0 1 2 3 4 5 6)
graph save density_cost_entry4_formal, replace
kdensity Wcost_entry5 if sector ~= 5 & sector ~= 0 & Wcost_entry5 <= 6, bwidth(0.2) n(50) lcolor(red) ytitle("") xtitle("Cost of Entry Trade From Formal Sector") xlabel(0 1 2 3 4 5 6)
graph save density_cost_entry5_formal, replace
kdensity Wcost_entry6 if sector ~= 6 & sector ~= 0 & Wcost_entry6 <= 6, bwidth(0.2) n(50) lcolor(red) ytitle("") xtitle("Cost of Entry Trans/Util From Formal Sector") xlabel(0 1 2 3 4 5 6)
graph save density_cost_entry6_formal, replace
kdensity Wcost_entry7 if sector ~= 7 & sector ~= 0 & Wcost_entry7 <= 6, bwidth(0.2) n(50) lcolor(red) ytitle("") xtitle("Cost of Entry Service From Formal Sector") xlabel(0 1 2 3 4 5 6)
graph save density_cost_entry7_formal, replace 	 	 



kdensity Wcost_entry1 if sector == 0 & Wcost_entry1 <= 25, bwidth(1.0) n(100) lcolor(red) ytitle("") xtitle("Cost of Entry Agriculture/Mining From Residual Sector") xlabel(0 5 10 15 20 25)
graph save density_cost_entry1_residual, replace
kdensity Wcost_entry2 if sector == 0 & Wcost_entry2 <= 25, bwidth(1.0) n(100) lcolor(red) ytitle("") xtitle("Cost of Entry Low Tech Manuf From Residual Sector") xlabel(0 5 10 15 20 25)
graph save density_cost_entry2_residual, replace
kdensity Wcost_entry3 if sector == 0 & Wcost_entry3 <= 25, bwidth(1.0) n(100) lcolor(red) ytitle("") xtitle("Cost of Entry High tech Manuf From Residual Sector") xlabel(0 5 10 15 20 25)
graph save density_cost_entry3_residual, replace
kdensity Wcost_entry4 if sector == 0 & Wcost_entry4 <= 25, bwidth(1.0) n(100) lcolor(red) ytitle("") xtitle("Cost of Entry Construction From Residual Sector") xlabel(0 5 10 15 20 25)
graph save density_cost_entry4_residual, replace
kdensity Wcost_entry5 if sector == 0 & Wcost_entry5 <= 25, bwidth(1.0) n(100) lcolor(red) ytitle("") xtitle("Cost of Entry Trade From Residual Sector") xlabel(0 5 10 15 20 25)
graph save density_cost_entry5_residual, replace
kdensity Wcost_entry6 if sector == 0 & Wcost_entry6 <= 25, bwidth(1.0) n(100) lcolor(red) ytitle("") xtitle("Cost of Entry Trans/Util From Residual Sector") xlabel(0 5 10 15 20 25)
graph save density_cost_entry6_residual, replace
kdensity Wcost_entry7 if sector == 0 & Wcost_entry7 <= 25, bwidth(1.0) n(100) lcolor(red) ytitle("") xtitle("Cost of Entry Service From Residual Sector") xlabel(0 5 10 15 20 25)
graph save density_cost_entry7_residual, replace

graph combine density_cost_entry1_formal.gph density_cost_entry2_formal.gph ///
density_cost_entry3_formal.gph density_cost_entry4_formal.gph density_cost_entry5_formal.gph ///
density_cost_entry6_formal.gph density_cost_entry7_formal.gph, saving(CostDensity_Formal, replace)

graph combine density_cost_entry1_residual.gph density_cost_entry2_residual.gph ///
density_cost_entry3_residual.gph density_cost_entry4_residual.gph density_cost_entry5_residual.gph ///
density_cost_entry6_residual.gph density_cost_entry7_residual.gph, saving(CostDensity_Residual, replace)

****************************************************************************************************************************************************************************************************
****************************************************************************************************************************************************************************************************

keep Wcost_entry* Vcost_entry* gender dummy_educ2 dummy_educ3 dummy_educ4 age age_2 exper1-exper7 sector
gen id = _n
sort id
save aux1, replace
keep Wcost_entry* Vcost_entry* id
reshape long Wcost_entry Vcost_entry, i(id) j(dest)
sort id
merge id using aux1
duplicates drop

gen log_Wcost = log(Wcost_entry)
gen log_Vcost = log(Vcost_entry)

sort sector dest 
egen gp = group(sector dest)

tab sector, gen(dummy_org)
tab dest, gen(dummy_dest)

reg log_Wcost dummy_org* dummy_dest* gender dummy_educ2 dummy_educ3 dummy_educ4 age age_2 exper1-exper7, noc 
outreg2 gender dummy_educ2 dummy_educ3 dummy_educ4 age age_2 exper1-exper7 using cost_reg, excel bracket bdec(5) noaster ctitle(Wcost) replace

reg log_Wcost dummy_org* dummy_dest* gender dummy_educ2 dummy_educ3 dummy_educ4 age age_2, noc 
outreg2 gender dummy_educ2 dummy_educ3 dummy_educ4 age age_2 using cost_reg, excel bracket bdec(5) noaster ctitle(Wcost) append

reg log_Vcost dummy_org* dummy_dest* gender dummy_educ2 dummy_educ3 dummy_educ4 age age_2 exper1-exper7, noc
outreg2 gender dummy_educ2 dummy_educ3 dummy_educ4 age age_2 exper1-exper7 using cost_reg, excel bracket bdec(5) noaster ctitle(Vcost) append

reg log_Vcost dummy_org* dummy_dest* gender dummy_educ2 dummy_educ3 dummy_educ4 age age_2, noc
outreg2 gender dummy_educ2 dummy_educ3 dummy_educ4 age age_2 using cost_reg, excel bracket bdec(5) noaster ctitle(Vcost) append



*********************************************************************************************************************
*********************************************************************************************************************


*****************************
* Costs incurred by switchers
*****************************

clear

cd "C:\Users\Dix-Carneiro\Documents\Trade and Labor Dynamics\Revision_Econometrica\CostsOfMobility"

insheet using "Data_set.csv"

rename v1  id
rename v2  coh 
rename v3  year
rename v4  educ
rename v5  gender
rename v6  choice
rename v7  choice1
rename v8  choice2
rename v9  choice3
rename v10 choice4
rename v11 choice5
rename v12 choice6
rename v13 choice7
rename v14 choice8
rename v15 choice9
rename v16 exper1
rename v17 exper2
rename v18 exper3
rename v19 exper4
rename v20 exper5
rename v21 exper6
rename v22 exper7
rename v23 wage
rename v24 value
rename v25 cost1
rename v26 cost2
rename v27 wgt
drop   v28

drop choice2-choice9

tab year, gen(year)
tab educ, gen(educ)

gen age = year - coh - 25
gen age_2 = age^2

replace gender = gender - 1

gen log_w_hat = wage_coef[1,1]*gender + wage_coef[1,2]*educ2 + wage_coef[1,3]*educ3 + ///
                wage_coef[1,4]*educ4 + wage_coef[1,5]*age + wage_coef[1,6]*age_2 + ///
				wage_coef[1,7]*exper1 + wage_coef[1,8]*exper2 + wage_coef[1,9]*exper3 + ///
				wage_coef[1,10]*exper4 + wage_coef[1,11]*exper5 + wage_coef[1,12]*exper6 + ///
				wage_coef[1,13]*exper7 + wage_coef[1,14]
				
gen w_hat = exp(log_w_hat)*exp(0.5*wage_rmse[1,1]^2)				

replace cost1 = cost1 / w_hat
replace cost2 = cost2 / w_hat

***************************************************
***************************************************

matrix cost1_switchers_formal = J(7,3,0)

forvalues i = 1(1)7{
	sum cost1 if choice == `i' & choice1 ~= `i' & choice1 ~= 0 [w =wgt], detail
	matrix cost1_switchers_formal[`i',1] = r(p50)
	matrix cost1_switchers_formal[`i',2] = r(p75)
	matrix cost1_switchers_formal[`i',3] = r(p90)
}

***************************************************
***************************************************

matrix cost1_switchers_residual = J(7,3,0)

forvalues i = 1(1)7{
	sum cost1 if choice == `i' & choice1 ~= `i' & choice1 == 0 [w =wgt], detail
	matrix cost1_switchers_residual[`i',1] = r(p50)
	matrix cost1_switchers_residual[`i',2] = r(p75)
	matrix cost1_switchers_residual[`i',3] = r(p90)
}	

***************************************************
***************************************************

cd "C:\Users\Dix-Carneiro\Documents\Trade and Labor Dynamics\Revision_Econometrica\CostsOfMobility\TablesGraphs"

clear
svmat cost1_switchers_formal

rename cost1_switchers_formal1 p50
rename cost1_switchers_formal2 p75
rename cost1_switchers_formal3 p90

outsheet _all using cost1_switchers_formal.csv, comma replace

clear
svmat cost1_switchers_residual

rename cost1_switchers_residual1 p50
rename cost1_switchers_residual2 p75
rename cost1_switchers_residual3 p90

outsheet _all using cost1_switchers_residual.csv, comma replace

erase rsk
