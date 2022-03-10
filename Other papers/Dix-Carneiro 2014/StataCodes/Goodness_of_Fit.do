clear
clear matrix

set matsize 5000


cd "C:\Users\Dix-Carneiro\Documents\Trade and Labor Dynamics\Revision_Econometrica\GoodnessOfFit"

********************************************************************************
* We first read the coefficients of the auxiliary model applied to the real data 
* and those of the auxiliary model applied to the simulated data
********************************************************************************

clear

insheet using "CompareCoefficients.csv"

rename v1 coef_type
rename v5 coefData
rename v6 coefSim
rename v7 sigmaData

gen weight = 1/(sigmaData^2)

cd "C:\Users\Dix-Carneiro\Documents\Trade and Labor Dynamics\Revision_Econometrica\GoodnessOfFit\TablesGraphs"

reg coefData coefSim [w=weight]
outreg2 coefSim using IndInfReg, excel bracket bdec(7) noaster replace

scatter coefData coefSim, ytitle("Coeficients Data") ///
xtitle("Coeficients Simulation") legend(off) msymbol(circle_hollow) color(blue) || ///
lfit coefData coefSim [w=weight], saving(IndirectInference1, replace)

gen coefData_st = coefData / sigmaData
gen coefSim_st = coefSim / sigmaData
gen ones = 1 / sigmaData

scatter coefData_st coefSim_st, ytitle("Standardized Coeficients Data") ///
xtitle("Standardized Coeficients Simulation") legend(off) msymbol(circle_hollow) color(blue) || ///
lfit coefData_st coefSim_st, saving(IndirectInference2, replace)




clear

cd "C:\Users\Dix-Carneiro\Documents\Trade and Labor Dynamics\Revision_Econometrica\GoodnessOfFit"

********************************************************************************
* We first read the coefficients of the auxiliary model applied to the real data 
* and those of the auxiliary model applied to the simulated data
********************************************************************************

clear

insheet using "CompareCoefficients.csv"

preserve

keep if v1 == "beta"

keep v5 v6

mkmat v5 v6 , matrix(beta_aux)

restore

matrix betaData = J(24,7,0)
matrix betaSim = J(24,7,0)

forvalues i = 1(1)24 {
	matrix betaData[`i',1] = beta_aux[`i',1]
	matrix betaData[`i',2] = beta_aux[24+`i',1]
	matrix betaData[`i',3] = beta_aux[48+`i',1]
	matrix betaData[`i',4] = beta_aux[72+`i',1]
	matrix betaData[`i',5] = beta_aux[96+`i',1]
	matrix betaData[`i',6] = beta_aux[120+`i',1]
	matrix betaData[`i',7] = beta_aux[144+`i',1]
	
	matrix betaSim[`i',1]  = beta_aux[`i',2]
	matrix betaSim[`i',2]  = beta_aux[24+`i',2]
	matrix betaSim[`i',3]  = beta_aux[48+`i',2]
	matrix betaSim[`i',4]  = beta_aux[72+`i',2]	
	matrix betaSim[`i',5]  = beta_aux[96+`i',2]	
	matrix betaSim[`i',6]  = beta_aux[120+`i',2]	
	matrix betaSim[`i',7]  = beta_aux[144+`i',2]	
}




preserve

keep if v1 == "gamma"

keep v5 v6

mkmat v5 v6 , matrix(gamma_aux)

restore


matrix gammaData = J(24,8,0)
matrix gammaSim  = J(24,8,0)

forvalues i = 1(1)24 {
	matrix gammaData[`i',1] = gamma_aux[`i',1]
	matrix gammaData[`i',2] = gamma_aux[24+`i',1]
	matrix gammaData[`i',3] = gamma_aux[48+`i',1]
	matrix gammaData[`i',4] = gamma_aux[72+`i',1]
	matrix gammaData[`i',5] = gamma_aux[96+`i',1]
	matrix gammaData[`i',6] = gamma_aux[120+`i',1]
	matrix gammaData[`i',7] = gamma_aux[144+`i',1]
	matrix gammaData[`i',8] = gamma_aux[168+`i',1]

	matrix gammaSim[`i',1]  = gamma_aux[`i',2]
	matrix gammaSim[`i',2]  = gamma_aux[24+`i',2]	
	matrix gammaSim[`i',3]  = gamma_aux[48+`i',2]
	matrix gammaSim[`i',4]  = gamma_aux[72+`i',2]
	matrix gammaSim[`i',5]  = gamma_aux[96+`i',2]
	matrix gammaSim[`i',6]  = gamma_aux[120+`i',2]
	matrix gammaSim[`i',7]  = gamma_aux[144+`i',2]
	matrix gammaSim[`i',8]  = gamma_aux[168+`i',2]
}



preserve

keep if v1 == "phi"

keep v5 v6

mkmat v5 v6 , matrix(phi_aux)

restore

matrix phiData = J(24,64,0)
matrix phiSim  = J(24,64,0)

forvalues s1 = 1(1)8{
	forvalues s2 = 1(1)8{
		forvalues i = 1(1)24{
		
			matrix phiData[`i',(`s1'-1)*8+`s2']  = phi_aux[((`s1'-1)*8+`s2'-1)*24+`i',1]
			matrix phiSim[`i',(`s1'-1)*8+`s2']  = phi_aux[((`s1'-1)*8+`s2'-1)*24+`i',2]
		
		}
	}
}




preserve

keep if v1 == "xsi1998"

keep v5 v6

mkmat v5 v6 , matrix(xsi1998_aux)

restore


matrix xsi1998Data = J(21,8,0)
matrix xsi1998Sim  = J(21,8,0)

forvalues i = 1(1)21 {
	matrix xsi1998Data[`i',1] = xsi1998_aux[`i',1]
	matrix xsi1998Data[`i',2] = xsi1998_aux[21+`i',1]
	matrix xsi1998Data[`i',3] = xsi1998_aux[42+`i',1]
	matrix xsi1998Data[`i',4] = xsi1998_aux[63+`i',1]
	matrix xsi1998Data[`i',5] = xsi1998_aux[84+`i',1]
	matrix xsi1998Data[`i',6] = xsi1998_aux[105+`i',1]
	matrix xsi1998Data[`i',7] = xsi1998_aux[126+`i',1]
	matrix xsi1998Data[`i',8] = xsi1998_aux[147+`i',1]

	matrix xsi1998Sim[`i',1]  = xsi1998_aux[`i',2]
	matrix xsi1998Sim[`i',2]  = xsi1998_aux[21+`i',2]	
	matrix xsi1998Sim[`i',3]  = xsi1998_aux[42+`i',2]
	matrix xsi1998Sim[`i',4]  = xsi1998_aux[63+`i',2]
	matrix xsi1998Sim[`i',5]  = xsi1998_aux[84+`i',2]
	matrix xsi1998Sim[`i',6]  = xsi1998_aux[105+`i',2]
	matrix xsi1998Sim[`i',7]  = xsi1998_aux[126+`i',2]
	matrix xsi1998Sim[`i',8]  = xsi1998_aux[147+`i',2]
}



preserve

keep if v1 == "xsi2000"

keep v5 v6

mkmat v5 v6 , matrix(xsi2000_aux)

restore


matrix xsi2000Data = J(21,8,0)
matrix xsi2000Sim  = J(21,8,0)

forvalues i = 1(1)21 {
	matrix xsi2000Data[`i',1] = xsi2000_aux[`i',1]
	matrix xsi2000Data[`i',2] = xsi2000_aux[21+`i',1]
	matrix xsi2000Data[`i',3] = xsi2000_aux[42+`i',1]
	matrix xsi2000Data[`i',4] = xsi2000_aux[63+`i',1]
	matrix xsi2000Data[`i',5] = xsi2000_aux[84+`i',1]
	matrix xsi2000Data[`i',6] = xsi2000_aux[105+`i',1]
	matrix xsi2000Data[`i',7] = xsi2000_aux[126+`i',1]
	matrix xsi2000Data[`i',8] = xsi2000_aux[147+`i',1]

	matrix xsi2000Sim[`i',1]  = xsi2000_aux[`i',2]
	matrix xsi2000Sim[`i',2]  = xsi2000_aux[21+`i',2]	
	matrix xsi2000Sim[`i',3]  = xsi2000_aux[42+`i',2]
	matrix xsi2000Sim[`i',4]  = xsi2000_aux[63+`i',2]
	matrix xsi2000Sim[`i',5]  = xsi2000_aux[84+`i',2]
	matrix xsi2000Sim[`i',6]  = xsi2000_aux[105+`i',2]
	matrix xsi2000Sim[`i',7]  = xsi2000_aux[126+`i',2]
	matrix xsi2000Sim[`i',8]  = xsi2000_aux[147+`i',2]
}


preserve

keep if v1 == "xsi2005"

keep v5 v6

mkmat v5 v6 , matrix(xsi2005_aux)

restore


matrix xsi2005Data = J(21,8,0)
matrix xsi2005Sim  = J(21,8,0)

forvalues i = 1(1)21 {
	matrix xsi2005Data[`i',1] = xsi2005_aux[`i',1]
	matrix xsi2005Data[`i',2] = xsi2005_aux[21+`i',1]
	matrix xsi2005Data[`i',3] = xsi2005_aux[42+`i',1]
	matrix xsi2005Data[`i',4] = xsi2005_aux[63+`i',1]
	matrix xsi2005Data[`i',5] = xsi2005_aux[84+`i',1]
	matrix xsi2005Data[`i',6] = xsi2005_aux[105+`i',1]
	matrix xsi2005Data[`i',7] = xsi2005_aux[126+`i',1]
	matrix xsi2005Data[`i',8] = xsi2005_aux[147+`i',1]

	matrix xsi2005Sim[`i',1]  = xsi2005_aux[`i',2]
	matrix xsi2005Sim[`i',2]  = xsi2005_aux[21+`i',2]	
	matrix xsi2005Sim[`i',3]  = xsi2005_aux[42+`i',2]
	matrix xsi2005Sim[`i',4]  = xsi2005_aux[63+`i',2]
	matrix xsi2005Sim[`i',5]  = xsi2005_aux[84+`i',2]
	matrix xsi2005Sim[`i',6]  = xsi2005_aux[105+`i',2]
	matrix xsi2005Sim[`i',7]  = xsi2005_aux[126+`i',2]
	matrix xsi2005Sim[`i',8]  = xsi2005_aux[147+`i',2]
}




preserve

keep if v1 == "eta"

keep v5 v6

mkmat v5 v6 , matrix(eta_aux)

restore


matrix etaData = J(21,8,0)
matrix etaSim  = J(21,8,0)

forvalues i = 1(1)21 {
	matrix etaData[`i',1] = eta_aux[`i',1]
	matrix etaData[`i',2] = eta_aux[21+`i',1]
	matrix etaData[`i',3] = eta_aux[42+`i',1]
	matrix etaData[`i',4] = eta_aux[63+`i',1]
	matrix etaData[`i',5] = eta_aux[84+`i',1]
	matrix etaData[`i',6] = eta_aux[105+`i',1]
	matrix etaData[`i',7] = eta_aux[126+`i',1]
	matrix etaData[`i',8] = eta_aux[147+`i',1]

	matrix etaSim[`i',1]  = eta_aux[`i',2]
	matrix etaSim[`i',2]  = eta_aux[21+`i',2]	
	matrix etaSim[`i',3]  = eta_aux[42+`i',2]
	matrix etaSim[`i',4]  = eta_aux[63+`i',2]
	matrix etaSim[`i',5]  = eta_aux[84+`i',2]
	matrix etaSim[`i',6]  = eta_aux[105+`i',2]
	matrix etaSim[`i',7]  = eta_aux[126+`i',2]
	matrix etaSim[`i',8]  = eta_aux[147+`i',2]
}



preserve

keep if v1 == "rho"

keep v5 v6

mkmat v5 v6 , matrix(rho_aux)

restore


matrix rhoData = J(24,8,0)
matrix rhoSim  = J(24,8,0)

forvalues i = 1(1)24 {
	matrix rhoData[`i',1] = rho_aux[`i',1]
	matrix rhoData[`i',2] = rho_aux[24+`i',1]
	matrix rhoData[`i',3] = rho_aux[48+`i',1]
	matrix rhoData[`i',4] = rho_aux[72+`i',1]
	matrix rhoData[`i',5] = rho_aux[96+`i',1]
	matrix rhoData[`i',6] = rho_aux[120+`i',1]
	matrix rhoData[`i',7] = rho_aux[144+`i',1]
	matrix rhoData[`i',8] = rho_aux[168+`i',1]

	matrix rhoSim[`i',1]  = rho_aux[`i',2]
	matrix rhoSim[`i',2]  = rho_aux[24+`i',2]	
	matrix rhoSim[`i',3]  = rho_aux[48+`i',2]
	matrix rhoSim[`i',4]  = rho_aux[72+`i',2]
	matrix rhoSim[`i',5]  = rho_aux[96+`i',2]
	matrix rhoSim[`i',6]  = rho_aux[120+`i',2]
	matrix rhoSim[`i',7]  = rho_aux[144+`i',2]
	matrix rhoSim[`i',8]  = rho_aux[168+`i',2]
}






****************************************************************************************
* We generate Wage Moments in the real data and those predicted by the model (Simulated)
* Then, we plot the Wage Moments against the Simulated Moments
****************************************************************************************

cd "C:\Users\Dix-Carneiro\Documents\Trade and Labor Dynamics\Revision_Econometrica\GoodnessOfFit"

clear

use New_Painel

keep pis ano sector

sort pis ano

reshape wide sector, i(pis) j(ano)

sort pis

gen sector0_1998 = (sector1998 == 0) 
gen sector1_1998 = (sector1998 == 1) 
gen sector2_1998 = (sector1998 == 2) 
gen sector3_1998 = (sector1998 == 3) 
gen sector4_1998 = (sector1998 == 4) 
gen sector5_1998 = (sector1998 == 5) 
gen sector6_1998 = (sector1998 == 6) 
gen sector7_1998 = (sector1998 == 7) 

gen sector0_2000 = (sector2000 == 0) 
gen sector1_2000 = (sector2000 == 1) 
gen sector2_2000 = (sector2000 == 2) 
gen sector3_2000 = (sector2000 == 3) 
gen sector4_2000 = (sector2000 == 4) 
gen sector5_2000 = (sector2000 == 5) 
gen sector6_2000 = (sector2000 == 6) 
gen sector7_2000 = (sector2000 == 7) 

gen sector0_2005 = (sector2005 == 0) 
gen sector1_2005 = (sector2005 == 1) 
gen sector2_2005 = (sector2005 == 2) 
gen sector3_2005 = (sector2005 == 3) 
gen sector4_2005 = (sector2005 == 4) 
gen sector5_2005 = (sector2005 == 5) 
gen sector6_2005 = (sector2005 == 6) 
gen sector7_2005 = (sector2005 == 7) 


gen freq0 = 0
gen freq1 = 0
gen freq2 = 0
gen freq3 = 0
gen freq4 = 0
gen freq5 = 0
gen freq6 = 0
gen freq7 = 0

forvalues i = 1995(1)2005{
	replace freq0 = freq0 + 1 if sector`i' == 0 
	replace freq1 = freq1 + 1 if sector`i' == 1
	replace freq2 = freq2 + 1 if sector`i' == 2
	replace freq3 = freq3 + 1 if sector`i' == 3
	replace freq4 = freq4 + 1 if sector`i' == 4
	replace freq5 = freq5 + 1 if sector`i' == 5
	replace freq6 = freq6 + 1 if sector`i' == 6
	replace freq7 = freq7 + 1 if sector`i' == 7
}

save choices_wide, replace


clear 

use New_Painel

gen dummy_educ2 = (educ_new == 2)
gen dummy_educ3 = (educ_new == 3)
gen dummy_educ4 = (educ_new == 4)

gen lag_sector0 = (lag_sector == 0)
gen lag_sector1 = (lag_sector == 1)
gen lag_sector2 = (lag_sector == 2)
gen lag_sector3 = (lag_sector == 3)
gen lag_sector4 = (lag_sector == 4)
gen lag_sector5 = (lag_sector == 5)
gen lag_sector6 = (lag_sector == 6)
gen lag_sector7 = (lag_sector == 7)

gen age = idade_new - 25
gen age_2 = (idade_new - 25)^2

preserve

forvalues i = 1(1)7 {

	gen Wage`i'Data =   betaData[1,`i']*year95        + ///
						betaData[2,`i']*year96        + ///
						betaData[3,`i']*year97        + ///
						betaData[4,`i']*year98        + ///
						betaData[5,`i']*year99        + ///
						betaData[6,`i']*year00        + ///
						betaData[7,`i']*year01        + ///
						betaData[8,`i']*year02        + ///
						betaData[9,`i']*year03        + /// 
						betaData[10,`i']*year04       + ///
						betaData[11,`i']*year05       + ///
						betaData[12,`i']*gender       + ///
						betaData[13,`i']*dummy_educ2  + ///
						betaData[14,`i']*dummy_educ3  + ///
						betaData[15,`i']*dummy_educ4  + ///
						betaData[16,`i']*age          + ///
						betaData[17,`i']*age_2        + ///
						betaData[18,`i']*exper1       + ///
						betaData[19,`i']*exper2       + ///
						betaData[20,`i']*exper3       + ///
						betaData[21,`i']*exper4       + ///
						betaData[22,`i']*exper5       + ///
						betaData[23,`i']*exper6       + ///
						betaData[24,`i']*exper7 


	gen Wage`i'Sim =    betaSim[1,`i']*year95          + ///
			     	    betaSim[2,`i']*year96          + ///
						betaSim[3,`i']*year97          + ///
						betaSim[4,`i']*year98          + ///
						betaSim[5,`i']*year99          + ///
						betaSim[6,`i']*year00          + ///
						betaSim[7,`i']*year01          + ///
						betaSim[8,`i']*year02          + ///
						betaSim[9,`i']*year03          + ///
						betaSim[10,`i']*year04         + ///
						betaSim[11,`i']*year05         + ///
						betaSim[12,`i']*gender         + ///
						betaSim[13,`i']*dummy_educ2    + ///
						betaSim[14,`i']*dummy_educ3    + ///
						betaSim[15,`i']*dummy_educ4    + ///
						betaSim[16,`i']*age            + ///
						betaSim[17,`i']*age_2          + ///
						betaSim[18,`i']*exper1         + ///
						betaSim[19,`i']*exper2         + ///
						betaSim[20,`i']*exper3         + ///
						betaSim[21,`i']*exper4         + ///
						betaSim[22,`i']*exper5         + ///
						betaSim[23,`i']*exper6         + ///
						betaSim[24,`i']*exper7

}

cd "C:\Users\Dix-Carneiro\Documents\Trade and Labor Dynamics\Revision_Econometrica\GoodnessOfFit\TablesGraphs"

scatter Wage1Data Wage1Sim if sector == 1, title("Agriculture/Mining") ///
ytitle("Best Linear Predictor Data") xtitle("Best Linear Predictor Model") msymbol(p) color(blue) legend(off) || ///
lfit Wage1Data Wage1Data if sector == 1, lcolor(red) lwidth(medium) legend(off) saving(wage1, replace) 

scatter Wage2Data Wage2Sim if sector == 2, title("LT Manuf") ///
ytitle("Best Linear Predictor Data") xtitle("Best Linear Predictor Model") msymbol(p) color(blue) legend(off) || ///
lfit Wage2Data Wage2Data if sector == 2, lcolor(red) lwidth(medium) legend(off) saving(wage2, replace) 

scatter Wage3Data Wage3Sim if sector == 3, title("HT Manuf") ///
ytitle("Best Linear Predictor Data") xtitle("Best Linear Predictor Model") msymbol(p) color(blue) legend(off) || ///
lfit Wage3Data Wage3Data if sector == 3, lcolor(red) lwidth(medium) legend(off) saving(wage3, replace) 

scatter Wage4Data Wage4Sim if sector == 4, title("Construction") ///
ytitle("Best Linear Predictor Data") xtitle("Best Linear Predictor Model") msymbol(p) color(blue) legend(off) || ///
lfit Wage4Data Wage4Data if sector == 4, lcolor(red) lwidth(medium) legend(off) saving(wage4, replace) 

scatter Wage5Data Wage5Sim if sector == 5, title("Trade") ///
ytitle("Best Linear Predictor Data") xtitle("Best Linear Predictor Model") msymbol(p) color(blue) legend(off) || ///
lfit Wage5Data Wage5Data if sector == 5, lcolor(red) lwidth(medium) legend(off) saving(wage5, replace) 

scatter Wage6Data Wage6Sim if sector == 6, title("Trans/Util") ///
ytitle("Best Linear Predictor Data") xtitle("Best Linear Predictor Model") msymbol(p) color(blue) legend(off) || ///
lfit Wage6Data Wage6Data if sector == 6, lcolor(red) lwidth(medium) legend(off) saving(wage6, replace) 

scatter Wage7Data Wage7Sim if sector == 7, title("Service") ///
ytitle("Best Linear Predictor Data") xtitle("Best Linear Predictor Model") msymbol(p) color(blue) legend(off) || ///
lfit Wage7Data Wage7Data if sector == 7, lcolor(red) lwidth(medium) legend(off) saving(wage7, replace) 


graph combine wage1.gph wage2.gph wage3.gph wage4.gph wage5.gph wage6.gph wage7.gph, saving(WageDataSim, replace)



reg Wage1Data Wage1Sim if sector == 1
outreg2 Wage1Sim using WageBLPDataSim, excel bracket bdec(7) noaster replace
reg Wage2Data Wage2Sim if sector == 2
outreg2 Wage2Sim using WageBLPDataSim, excel bracket bdec(7) noaster append
reg Wage3Data Wage3Sim if sector == 3
outreg2 Wage3Sim using WageBLPDataSim, excel bracket bdec(7) noaster append
reg Wage4Data Wage4Sim if sector == 4 
outreg2 Wage4Sim using WageBLPDataSim, excel bracket bdec(7) noaster append
reg Wage5Data Wage5Sim if sector == 5 
outreg2 Wage5Sim using WageBLPDataSim, excel bracket bdec(7) noaster append
reg Wage6Data Wage6Sim if sector == 6 
outreg2 Wage6Sim using WageBLPDataSim, excel bracket bdec(7) noaster append
reg Wage7Data Wage7Sim if sector == 7
outreg2 Wage7Sim using WageBLPDataSim, excel bracket bdec(7) noaster append

restore

********************
* Employment Moments
********************

preserve

forvalues i = 0(1)7 {

	gen Prob`i'Data =   gammaData[1,1+`i']*year95       + ///
						gammaData[2,1+`i']*year96       + ///
						gammaData[3,1+`i']*year97       + ///
						gammaData[4,1+`i']*year98       + ///
						gammaData[5,1+`i']*year99       + ///
						gammaData[6,1+`i']*year00       + ///
						gammaData[7,1+`i']*year01       + ///
						gammaData[8,1+`i']*year02       + ///
						gammaData[9,1+`i']*year03       + ///
						gammaData[10,1+`i']*year04      + ///
						gammaData[11,1+`i']*year05      + ///
						gammaData[12,1+`i']*gender      + ///
						gammaData[13,1+`i']*dummy_educ2 + ///
						gammaData[14,1+`i']*dummy_educ3 + ///
						gammaData[15,1+`i']*dummy_educ4 + ///
						gammaData[16,1+`i']*age         + ///
						gammaData[17,1+`i']*age_2       + ///
						gammaData[18,1+`i']*exper1      + ///
						gammaData[19,1+`i']*exper2      + ///
						gammaData[20,1+`i']*exper3      + ///
						gammaData[21,1+`i']*exper4      + ///
						gammaData[22,1+`i']*exper5      + ///
						gammaData[23,1+`i']*exper6      + ///
						gammaData[24,1+`i']*exper7

	gen Prob`i'Sim =    gammaSim[1,1+`i']*year95       + ///
						gammaSim[2,1+`i']*year96       + ///
						gammaSim[3,1+`i']*year97       + ///
						gammaSim[4,1+`i']*year98       + ///
						gammaSim[5,1+`i']*year99       + ///
						gammaSim[6,1+`i']*year00       + ///
						gammaSim[7,1+`i']*year01       + ///
						gammaSim[8,1+`i']*year02       + ///
						gammaSim[9,1+`i']*year03       + ///
						gammaSim[10,1+`i']*year04      + ///
						gammaSim[11,1+`i']*year05      + ///
						gammaSim[12,1+`i']*gender      + ///
						gammaSim[13,1+`i']*dummy_educ2 + ///
						gammaSim[14,1+`i']*dummy_educ3 + ///
						gammaSim[15,1+`i']*dummy_educ4 + ///
						gammaSim[16,1+`i']*age         + ///
						gammaSim[17,1+`i']*age_2       + ///
						gammaSim[18,1+`i']*exper1      + ///
						gammaSim[19,1+`i']*exper2      + ///
						gammaSim[20,1+`i']*exper3      + ///
						gammaSim[21,1+`i']*exper4      + ///
						gammaSim[22,1+`i']*exper5      + ///
						gammaSim[23,1+`i']*exper6      + ///
						gammaSim[24,1+`i']*exper7      

}

scatter Prob0Data Prob0Sim, title("Residual Sector") ytitle("Best Linear Predictor Data") ///
xtitle("Best Linear Predictor Model") legend(off) msymbol(p) color(blue) || ///
lfit Prob0Data Prob0Data, range(-0.2 1.2) lcolor(red) lwidth(medium) saving(Prob0, replace)

scatter Prob1Data Prob1Sim, title("Agriculture/Mining") ytitle("Best Linear Predictor Data") ///
xtitle("Best Linear Predictor Model") legend(off) msymbol(p) color(blue) || ///
lfit Prob1Data Prob1Data, range(-0.2 1.2) lcolor(red) lwidth(medium) saving(Prob1, replace)

scatter Prob2Data Prob2Sim, title("LT Manuf") ytitle("Best Linear Predictor Data") ///
xtitle("Best Linear Predictor Model") legend(off) msymbol(p) color(blue) || ///
lfit Prob2Data Prob2Data, range(-0.2 1.2) lcolor(red) lwidth(medium) saving(Prob2, replace)

scatter Prob3Data Prob3Sim, title("HT Manuf") ytitle("Best Linear Predictor Data") ///
xtitle("Best Linear Predictor Model") legend(off) msymbol(p) color(blue) || ///
lfit Prob3Data Prob3Data, range(-0.2 1.2) lcolor(red) lwidth(medium) saving(Prob3, replace)

scatter Prob4Data Prob4Sim, title("Construction") ytitle("Best Linear Predictor Data") ///
xtitle("Best Linear Predictor Model") legend(off) msymbol(p) color(blue) || ///
lfit Prob4Data Prob4Data, range(-0.2 1.2) lcolor(red) lwidth(medium) saving(Prob4, replace)

scatter Prob5Data Prob5Sim, title("Trade") ytitle("Best Linear Predictor Data") ///
xtitle("Best Linear Predictor Model") legend(off) msymbol(p) color(blue) || ///
lfit Prob5Data Prob5Data, range(-0.2 1.2) lcolor(red) lwidth(medium) saving(Prob5, replace)

scatter Prob6Data Prob6Sim, title("Trans/Util") ytitle("Best Linear Predictor Data") ///
xtitle("Best Linear Predictor Model") legend(off) msymbol(p) color(blue) || ///
lfit Prob6Data Prob6Data, range(-0.2 1.2) lcolor(red) lwidth(medium) saving(Prob6, replace)

scatter Prob7Data Prob7Sim, title("Service") ytitle("Best Linear Predictor Data") ///
xtitle("Best Linear Predictor Model") legend(off) msymbol(p) color(blue) || ///
lfit Prob7Data Prob7Data, range(-0.2 1.2) lcolor(red) lwidth(medium) saving(Prob7, replace) 	 	 	 

graph combine Prob0.gph Prob1.gph Prob2.gph Prob3.gph Prob4.gph Prob5.gph Prob6.gph Prob7.gph, saving(EmpDataSim, replace)

reg Prob0Data Prob0Sim
outreg2 Prob0Sim using ProbBLPDataSim, excel bracket bdec(7) noaster replace
reg Prob1Data Prob1Sim
outreg2 Prob1Sim using ProbBLPDataSim, excel bracket bdec(7) noaster append
reg Prob2Data Prob2Sim
outreg2 Prob2Sim using ProbBLPDataSim, excel bracket bdec(7) noaster append
reg Prob3Data Prob3Sim
outreg2 Prob3Sim using ProbBLPDataSim, excel bracket bdec(7) noaster append
reg Prob4Data Prob4Sim
outreg2 Prob4Sim using ProbBLPDataSim, excel bracket bdec(7) noaster append
reg Prob5Data Prob5Sim
outreg2 Prob5Sim using ProbBLPDataSim, excel bracket bdec(7) noaster append
reg Prob6Data Prob6Sim
outreg2 Prob6Sim using ProbBLPDataSim, excel bracket bdec(7) noaster append
reg Prob7Data Prob7Sim
outreg2 Prob7Sim using ProbBLPDataSim, excel bracket bdec(7) noaster append

restore


******************
* Transition Rates
******************

preserve

forvalues s1 = 0(1)7 {
	forvalues s2 = 0(1)7 {

		gen Tr`s1'`s2'Data = phiData[1,`s1'*8+`s2' + 1]*year95 		 + ///
			     		     phiData[2,`s1'*8+`s2' + 1]*year96 		 + ///
					         phiData[3,`s1'*8+`s2' + 1]*year97 		 + ///
							 phiData[4,`s1'*8+`s2' + 1]*year98 		 + ///
							 phiData[5,`s1'*8+`s2' + 1]*year99 		 + ///
							 phiData[6,`s1'*8+`s2' + 1]*year00 		 + ///
							 phiData[7,`s1'*8+`s2' + 1]*year01 		 + ///
							 phiData[8,`s1'*8+`s2' + 1]*year02 		 + ///
							 phiData[9,`s1'*8+`s2' + 1]*year03 		 + ///
							 phiData[10,`s1'*8+`s2' + 1]*year04 	 + ///
							 phiData[11,`s1'*8+`s2' + 1]*year05 	 + ///
							 phiData[12,`s1'*8+`s2' + 1]*gender 	 + ///
							 phiData[13,`s1'*8+`s2' + 1]*dummy_educ2 + ///
							 phiData[14,`s1'*8+`s2' + 1]*dummy_educ3 + ///
							 phiData[15,`s1'*8+`s2' + 1]*dummy_educ4 + ///
							 phiData[16,`s1'*8+`s2' + 1]*age         + ///
							 phiData[17,`s1'*8+`s2' + 1]*age_2 		 + ///
							 phiData[18,`s1'*8+`s2' + 1]*exper1 	 + ///
							 phiData[19,`s1'*8+`s2' + 1]*exper2 	 + ///
							 phiData[20,`s1'*8+`s2' + 1]*exper3 	 + ///
							 phiData[21,`s1'*8+`s2' + 1]*exper4 	 + ///
							 phiData[22,`s1'*8+`s2' + 1]*exper5 	 + ///
							 phiData[23,`s1'*8+`s2' + 1]*exper6 	 + ///
							 phiData[24,`s1'*8+`s2' + 1]*exper7
							 
							 

		gen Tr`s1'`s2'Sim =  phiSim[1,`s1'*8+`s2' + 1]*year95 	    + ///
							 phiSim[2,`s1'*8+`s2' + 1]*year96 		+ ///
							 phiSim[3,`s1'*8+`s2' + 1]*year97 		+ ///
							 phiSim[4,`s1'*8+`s2' + 1]*year98 		+ ///
							 phiSim[5,`s1'*8+`s2' + 1]*year99 		+ ///
							 phiSim[6,`s1'*8+`s2' + 1]*year00 		+ ///
							 phiSim[7,`s1'*8+`s2' + 1]*year01 		+ ///
							 phiSim[8,`s1'*8+`s2' + 1]*year02 		+ ///
							 phiSim[9,`s1'*8+`s2' + 1]*year03 		+ ///
							 phiSim[10,`s1'*8+`s2' + 1]*year04 		+ ///
							 phiSim[11,`s1'*8+`s2' + 1]*year05 		+ ///
							 phiSim[12,`s1'*8+`s2' + 1]*gender 		+ ///
							 phiSim[13,`s1'*8+`s2' + 1]*dummy_educ2 + ///
							 phiSim[14,`s1'*8+`s2' + 1]*dummy_educ3 + ///
							 phiSim[15,`s1'*8+`s2' + 1]*dummy_educ4 + ///
							 phiSim[16,`s1'*8+`s2' + 1]*age 		+ ///
							 phiSim[17,`s1'*8+`s2' + 1]*age_2 		+ ///
							 phiSim[18,`s1'*8+`s2' + 1]*exper1 		+ ///
							 phiSim[19,`s1'*8+`s2' + 1]*exper2 		+ ///
							 phiSim[20,`s1'*8+`s2' + 1]*exper3 		+ ///
							 phiSim[21,`s1'*8+`s2' + 1]*exper4 		+ ///
							 phiSim[22,`s1'*8+`s2' + 1]*exper5 		+ ///
							 phiSim[23,`s1'*8+`s2' + 1]*exper6 		+ ///
							 phiSim[24,`s1'*8+`s2' + 1]*exper7

		scatter Tr`s1'`s2'Data Tr`s1'`s2'Sim if lag_sector == `s1', title("Tr`s1'`s2'") ytitle("Best Linear Predictor Data") xtitle("Best Linear Predictor Model") ///
            legend(off)  msymbol(p) color(blue) || ///
		lfit Tr`s1'`s2'Data Tr`s1'`s2'Data if lag_sector == `s1', range(-0.2 1.2) lcolor(red) lwidth(medium) saving(Tr`s1'`s2', replace)

	}
}
			


graph combine Tr00.gph Tr01.gph Tr02.gph Tr03.gph Tr04.gph Tr05.gph Tr06.gph Tr07.gph ///
		      Tr10.gph Tr11.gph Tr12.gph Tr13.gph Tr14.gph Tr15.gph Tr16.gph Tr17.gph ///
		      Tr20.gph Tr21.gph Tr22.gph Tr23.gph Tr24.gph Tr25.gph Tr26.gph Tr27.gph ///
		      Tr30.gph Tr31.gph Tr32.gph Tr33.gph Tr34.gph Tr45.gph Tr36.gph Tr37.gph ///
              Tr40.gph Tr41.gph Tr42.gph Tr43.gph Tr44.gph Tr45.gph Tr46.gph Tr47.gph ///
			  Tr50.gph Tr51.gph Tr52.gph Tr53.gph Tr54.gph Tr55.gph Tr56.gph Tr57.gph ///
			  Tr60.gph Tr61.gph Tr62.gph Tr63.gph Tr64.gph Tr65.gph Tr66.gph Tr67.gph ///
			  Tr70.gph Tr71.gph Tr72.gph Tr73.gph Tr74.gph Tr75.gph Tr76.gph Tr77.gph, ///
              saving(TrDataSim, replace)



			  
reg Tr00Data Tr00Sim if lag_sector == 0
outreg2 Tr00Sim using TrDataSim, excel bracket bdec(7) noaster replace	  
reg Tr01Data Tr01Sim if lag_sector == 0
outreg2 Tr01Sim using TrDataSim, excel bracket bdec(7) noaster append			  
reg Tr02Data Tr02Sim if lag_sector == 0
outreg2 Tr02Sim using TrDataSim, excel bracket bdec(7) noaster append			  
reg Tr03Data Tr03Sim if lag_sector == 0
outreg2 Tr03Sim using TrDataSim, excel bracket bdec(7) noaster append			  
reg Tr04Data Tr04Sim if lag_sector == 0
outreg2 Tr04Sim using TrDataSim, excel bracket bdec(7) noaster append		
reg Tr05Data Tr05Sim if lag_sector == 0
outreg2 Tr05Sim using TrDataSim, excel bracket bdec(7) noaster append
reg Tr06Data Tr06Sim if lag_sector == 0
outreg2 Tr06Sim using TrDataSim, excel bracket bdec(7) noaster append
reg Tr07Data Tr07Sim if lag_sector == 0
outreg2 Tr07Sim using TrDataSim, excel bracket bdec(7) noaster append



reg Tr10Data Tr10Sim if lag_sector == 1
outreg2 Tr10Sim using TrDataSim, excel bracket bdec(7) noaster append	  
reg Tr11Data Tr11Sim if lag_sector == 1
outreg2 Tr11Sim using TrDataSim, excel bracket bdec(7) noaster append			  
reg Tr12Data Tr12Sim if lag_sector == 1
outreg2 Tr12Sim using TrDataSim, excel bracket bdec(7) noaster append			  
reg Tr13Data Tr13Sim if lag_sector == 1
outreg2 Tr13Sim using TrDataSim, excel bracket bdec(7) noaster append			  
reg Tr14Data Tr14Sim if lag_sector == 1
outreg2 Tr14Sim using TrDataSim, excel bracket bdec(7) noaster append
reg Tr15Data Tr15Sim if lag_sector == 1
outreg2 Tr15Sim using TrDataSim, excel bracket bdec(7) noaster append
reg Tr16Data Tr16Sim if lag_sector == 1
outreg2 Tr16Sim using TrDataSim, excel bracket bdec(7) noaster append
reg Tr17Data Tr17Sim if lag_sector == 1
outreg2 Tr17Sim using TrDataSim, excel bracket bdec(7) noaster append	  


reg Tr20Data Tr20Sim if lag_sector == 2
outreg2 Tr20Sim using TrDataSim, excel bracket bdec(7) noaster append	  
reg Tr21Data Tr21Sim if lag_sector == 2
outreg2 Tr21Sim using TrDataSim, excel bracket bdec(7) noaster append			  
reg Tr22Data Tr22Sim if lag_sector == 2
outreg2 Tr22Sim using TrDataSim, excel bracket bdec(7) noaster append			  
reg Tr23Data Tr23Sim if lag_sector == 2
outreg2 Tr23Sim using TrDataSim, excel bracket bdec(7) noaster append			  
reg Tr24Data Tr24Sim if lag_sector == 2
outreg2 Tr24Sim using TrDataSim, excel bracket bdec(7) noaster append
reg Tr25Data Tr25Sim if lag_sector == 2
outreg2 Tr25Sim using TrDataSim, excel bracket bdec(7) noaster append
reg Tr26Data Tr26Sim if lag_sector == 2
outreg2 Tr26Sim using TrDataSim, excel bracket bdec(7) noaster append
reg Tr27Data Tr27Sim if lag_sector == 2
outreg2 Tr27Sim using TrDataSim, excel bracket bdec(7) noaster append


reg Tr30Data Tr30Sim if lag_sector == 3
outreg2 Tr30Sim using TrDataSim, excel bracket bdec(7) noaster append	  
reg Tr31Data Tr31Sim if lag_sector == 3
outreg2 Tr31Sim using TrDataSim, excel bracket bdec(7) noaster append			  
reg Tr32Data Tr32Sim if lag_sector == 3
outreg2 Tr32Sim using TrDataSim, excel bracket bdec(7) noaster append			  
reg Tr33Data Tr33Sim if lag_sector == 3
outreg2 Tr33Sim using TrDataSim, excel bracket bdec(7) noaster append			  
reg Tr34Data Tr34Sim if lag_sector == 3
outreg2 Tr34Sim using TrDataSim, excel bracket bdec(7) noaster append
reg Tr35Data Tr35Sim if lag_sector == 3
outreg2 Tr35Sim using TrDataSim, excel bracket bdec(7) noaster append
reg Tr36Data Tr36Sim if lag_sector == 3
outreg2 Tr36Sim using TrDataSim, excel bracket bdec(7) noaster append
reg Tr37Data Tr37Sim if lag_sector == 3
outreg2 Tr37Sim using TrDataSim, excel bracket bdec(7) noaster append


reg Tr40Data Tr40Sim if lag_sector == 4
outreg2 Tr40Sim using TrDataSim, excel bracket bdec(7) noaster append	  
reg Tr41Data Tr41Sim if lag_sector == 4
outreg2 Tr41Sim using TrDataSim, excel bracket bdec(7) noaster append			  
reg Tr42Data Tr42Sim if lag_sector == 4
outreg2 Tr42Sim using TrDataSim, excel bracket bdec(7) noaster append			  
reg Tr43Data Tr43Sim if lag_sector == 4
outreg2 Tr43Sim using TrDataSim, excel bracket bdec(7) noaster append			  
reg Tr44Data Tr44Sim if lag_sector == 4
outreg2 Tr44Sim using TrDataSim, excel bracket bdec(7) noaster append
reg Tr45Data Tr45Sim if lag_sector == 4
outreg2 Tr45Sim using TrDataSim, excel bracket bdec(7) noaster append
reg Tr46Data Tr46Sim if lag_sector == 4
outreg2 Tr46Sim using TrDataSim, excel bracket bdec(7) noaster append
reg Tr47Data Tr47Sim if lag_sector == 4
outreg2 Tr47Sim using TrDataSim, excel bracket bdec(7) noaster append


reg Tr50Data Tr50Sim if lag_sector == 5
outreg2 Tr50Sim using TrDataSim, excel bracket bdec(7) noaster append	  
reg Tr51Data Tr51Sim if lag_sector == 5
outreg2 Tr51Sim using TrDataSim, excel bracket bdec(7) noaster append			  
reg Tr52Data Tr52Sim if lag_sector == 5
outreg2 Tr52Sim using TrDataSim, excel bracket bdec(7) noaster append			  
reg Tr53Data Tr53Sim if lag_sector == 5
outreg2 Tr53Sim using TrDataSim, excel bracket bdec(7) noaster append			  
reg Tr54Data Tr54Sim if lag_sector == 5
outreg2 Tr54Sim using TrDataSim, excel bracket bdec(7) noaster append
reg Tr55Data Tr55Sim if lag_sector == 5
outreg2 Tr55Sim using TrDataSim, excel bracket bdec(7) noaster append
reg Tr56Data Tr56Sim if lag_sector == 5
outreg2 Tr56Sim using TrDataSim, excel bracket bdec(7) noaster append
reg Tr57Data Tr57Sim if lag_sector == 5
outreg2 Tr57Sim using TrDataSim, excel bracket bdec(7) noaster append


reg Tr60Data Tr60Sim if lag_sector == 6
outreg2 Tr60Sim using TrDataSim, excel bracket bdec(7) noaster append	  
reg Tr61Data Tr61Sim if lag_sector == 6
outreg2 Tr61Sim using TrDataSim, excel bracket bdec(7) noaster append			  
reg Tr62Data Tr52Sim if lag_sector == 6
outreg2 Tr62Sim using TrDataSim, excel bracket bdec(7) noaster append			  
reg Tr63Data Tr63Sim if lag_sector == 6
outreg2 Tr63Sim using TrDataSim, excel bracket bdec(7) noaster append			  
reg Tr64Data Tr64Sim if lag_sector == 6
outreg2 Tr64Sim using TrDataSim, excel bracket bdec(7) noaster append
reg Tr65Data Tr65Sim if lag_sector == 6
outreg2 Tr65Sim using TrDataSim, excel bracket bdec(7) noaster append
reg Tr66Data Tr66Sim if lag_sector == 6
outreg2 Tr66Sim using TrDataSim, excel bracket bdec(7) noaster append
reg Tr67Data Tr67Sim if lag_sector == 6
outreg2 Tr67Sim using TrDataSim, excel bracket bdec(7) noaster append


reg Tr70Data Tr70Sim if lag_sector == 7
outreg2 Tr70Sim using TrDataSim, excel bracket bdec(7) noaster append	  
reg Tr71Data Tr71Sim if lag_sector == 7
outreg2 Tr61Sim using TrDataSim, excel bracket bdec(7) noaster append			  
reg Tr72Data Tr72Sim if lag_sector == 7
outreg2 Tr72Sim using TrDataSim, excel bracket bdec(7) noaster append			  
reg Tr73Data Tr73Sim if lag_sector == 7
outreg2 Tr73Sim using TrDataSim, excel bracket bdec(7) noaster append			  
reg Tr74Data Tr74Sim if lag_sector == 7
outreg2 Tr64Sim using TrDataSim, excel bracket bdec(7) noaster append
reg Tr75Data Tr75Sim if lag_sector == 7
outreg2 Tr75Sim using TrDataSim, excel bracket bdec(7) noaster append
reg Tr76Data Tr76Sim if lag_sector == 7
outreg2 Tr76Sim using TrDataSim, excel bracket bdec(7) noaster append
reg Tr77Data Tr77Sim if lag_sector == 7
outreg2 Tr77Sim using TrDataSim, excel bracket bdec(7) noaster append

restore

cd "C:\Users\Dix-Carneiro\Documents\Trade and Labor Dynamics\Revision_Econometrica\GoodnessOfFit"

********************
* Persistence 1998
********************

sort pis

merge pis using choices_wide

cd "C:\Users\Dix-Carneiro\Documents\Trade and Labor Dynamics\Revision_Econometrica\GoodnessOfFit\TablesGraphs"

preserve

keep if idade_new >= 25 & idade_new <= 60-3 & ano == 1995

forvalues i = 0(1)7 {

	gen Pers1998_`i'Data =  xsi1998Data[1,1+`i']*lag_sector0        + ///
							xsi1998Data[2,1+`i']*lag_sector1        + ///
							xsi1998Data[3,1+`i']*lag_sector2        + ///
							xsi1998Data[4,1+`i']*lag_sector3        + ///
							xsi1998Data[5,1+`i']*lag_sector4        + ///
							xsi1998Data[6,1+`i']*lag_sector5        + ///
							xsi1998Data[7,1+`i']*lag_sector6        + ///
							xsi1998Data[8,1+`i']*lag_sector7        + ///
							xsi1998Data[9,1+`i']*gender      		+ ///
							xsi1998Data[10,1+`i']*dummy_educ2 		+ ///
							xsi1998Data[11,1+`i']*dummy_educ3 		+ ///
							xsi1998Data[12,1+`i']*dummy_educ4 		+ ///
							xsi1998Data[13,1+`i']*age         		+ ///
							xsi1998Data[14,1+`i']*age_2       		+ ///
							xsi1998Data[15,1+`i']*exper1      		+ ///
							xsi1998Data[16,1+`i']*exper2      		+ ///
							xsi1998Data[17,1+`i']*exper3      		+ ///
							xsi1998Data[18,1+`i']*exper4      		+ ///
							xsi1998Data[19,1+`i']*exper5      		+ ///
							xsi1998Data[20,1+`i']*exper6      		+ ///
							xsi1998Data[21,1+`i']*exper7

	gen Pers1998_`i'Sim =   xsi1998Sim[1,1+`i']*lag_sector0         + ///
							xsi1998Sim[2,1+`i']*lag_sector1         + ///
							xsi1998Sim[3,1+`i']*lag_sector2         + ///
							xsi1998Sim[4,1+`i']*lag_sector3         + ///
							xsi1998Sim[5,1+`i']*lag_sector4         + ///
							xsi1998Sim[6,1+`i']*lag_sector5         + ///
							xsi1998Sim[7,1+`i']*lag_sector6         + ///
							xsi1998Sim[8,1+`i']*lag_sector7         + ///
							xsi1998Sim[9,1+`i']*gender       		+ ///
							xsi1998Sim[10,1+`i']*dummy_educ2 		+ ///
							xsi1998Sim[11,1+`i']*dummy_educ3 		+ ///
							xsi1998Sim[12,1+`i']*dummy_educ4 		+ ///
							xsi1998Sim[13,1+`i']*age         		+ ///
							xsi1998Sim[14,1+`i']*age_2       		+ ///
							xsi1998Sim[15,1+`i']*exper1      		+ ///
							xsi1998Sim[16,1+`i']*exper2      		+ ///
							xsi1998Sim[17,1+`i']*exper3      		+ ///
							xsi1998Sim[18,1+`i']*exper4      		+ ///
							xsi1998Sim[19,1+`i']*exper5      		+ ///
							xsi1998Sim[20,1+`i']*exper6      		+ ///
							xsi1998Sim[21,1+`i']*exper7      

}

scatter Pers1998_0Data Pers1998_0Sim, title("Residual Sector") ytitle("Best Linear Predictor Data") ///
xtitle("Best Linear Predictor Model") legend(off) msymbol(p) color(blue) || ///
lfit Pers1998_0Data Pers1998_0Data, range(-0.2 1.2) lcolor(red) lwidth(medium) saving(Pers1998_0, replace)

scatter Pers1998_1Data Pers1998_1Sim, title("Agriculture/Mining") ytitle("Best Linear Predictor Data") ///
xtitle("Best Linear Predictor Model") legend(off) msymbol(p) color(blue) || ///
lfit Pers1998_1Data Pers1998_1Data, range(-0.2 1.2) lcolor(red) lwidth(medium) saving(Pers1998_1, replace)

scatter Pers1998_2Data Pers1998_2Sim, title("LT Manuf") ytitle("Best Linear Predictor Data") ///
xtitle("Best Linear Predictor Model") legend(off) msymbol(p) color(blue) || ///
lfit Pers1998_2Data Pers1998_2Data, range(-0.2 1.2) lcolor(red) lwidth(medium) saving(Pers1998_2, replace)

scatter Pers1998_3Data Pers1998_3Sim, title("HT Manuf") ytitle("Best Linear Predictor Data") ///
xtitle("Best Linear Predictor Model") legend(off) msymbol(p) color(blue) || ///
lfit Pers1998_3Data Pers1998_3Data, range(-0.2 1.2) lcolor(red) lwidth(medium) saving(Pers1998_3, replace)

scatter Pers1998_4Data Pers1998_4Sim, title("Construction") ytitle("Best Linear Predictor Data") ///
xtitle("Best Linear Predictor Model") legend(off) msymbol(p) color(blue) || ///
lfit Pers1998_4Data Pers1998_4Data, range(-0.2 1.2) lcolor(red) lwidth(medium) saving(Pers1998_4, replace)

scatter Pers1998_5Data Pers1998_5Sim, title("Trade") ytitle("Best Linear Predictor Data") ///
xtitle("Best Linear Predictor Model") legend(off) msymbol(p) color(blue) || ///
lfit Pers1998_5Data Pers1998_5Data, range(-0.2 1.2) lcolor(red) lwidth(medium) saving(Pers1998_5, replace)

scatter Pers1998_6Data Pers1998_6Sim, title("Trans/Util") ytitle("Best Linear Predictor Data") ///
xtitle("Best Linear Predictor Model") legend(off) msymbol(p) color(blue) || ///
lfit Pers1998_6Data Pers1998_6Data, range(-0.2 1.2) lcolor(red) lwidth(medium) saving(Pers1998_6, replace)

scatter Pers1998_7Data Pers1998_7Sim, title("Service") ytitle("Best Linear Predictor Data") ///
xtitle("Best Linear Predictor Model") legend(off) msymbol(p) color(blue) || ///
lfit Pers1998_7Data Pers1998_7Data, range(-0.2 1.2) lcolor(red) lwidth(medium) saving(Pers1998_7, replace) 

graph combine Pers1998_0.gph Pers1998_1.gph Pers1998_2.gph Pers1998_3.gph Pers1998_4.gph Pers1998_5.gph Pers1998_6.gph Pers1998_7.gph, saving(Pers1998DataSim, replace)

restore


********************
* Persistence 2000
********************

preserve

keep if idade_new >= 25 & idade_new <= 60-5 & ano == 1995

forvalues i = 0(1)7 {

	gen Pers2000_`i'Data =  xsi2000Data[1,1+`i']*lag_sector0        + ///
							xsi2000Data[2,1+`i']*lag_sector1        + ///
							xsi2000Data[3,1+`i']*lag_sector2        + ///
							xsi2000Data[4,1+`i']*lag_sector3        + ///
							xsi2000Data[5,1+`i']*lag_sector4        + ///
							xsi2000Data[6,1+`i']*lag_sector5        + ///
							xsi2000Data[7,1+`i']*lag_sector6        + ///
							xsi2000Data[8,1+`i']*lag_sector7        + ///
							xsi2000Data[9,1+`i']*gender      		+ ///
							xsi2000Data[10,1+`i']*dummy_educ2 		+ ///
							xsi2000Data[11,1+`i']*dummy_educ3 		+ ///
							xsi2000Data[12,1+`i']*dummy_educ4 		+ ///
							xsi2000Data[13,1+`i']*age         		+ ///
							xsi2000Data[14,1+`i']*age_2       		+ ///
							xsi2000Data[15,1+`i']*exper1      		+ ///
							xsi2000Data[16,1+`i']*exper2      		+ ///
							xsi2000Data[17,1+`i']*exper3      		+ ///
							xsi2000Data[18,1+`i']*exper4      		+ ///
							xsi2000Data[19,1+`i']*exper5      		+ ///
							xsi2000Data[20,1+`i']*exper6      		+ ///
							xsi2000Data[21,1+`i']*exper7

	gen Pers2000_`i'Sim =   xsi2000Sim[1,1+`i']*lag_sector0         + ///
							xsi2000Sim[2,1+`i']*lag_sector1         + ///
							xsi2000Sim[3,1+`i']*lag_sector2         + ///
							xsi2000Sim[4,1+`i']*lag_sector3         + ///
							xsi2000Sim[5,1+`i']*lag_sector4         + ///
							xsi2000Sim[6,1+`i']*lag_sector5         + ///
							xsi2000Sim[7,1+`i']*lag_sector6         + ///
							xsi2000Sim[8,1+`i']*lag_sector7         + ///
							xsi2000Sim[9,1+`i']*gender       		+ ///
							xsi2000Sim[10,1+`i']*dummy_educ2 		+ ///
							xsi2000Sim[11,1+`i']*dummy_educ3 		+ ///
							xsi2000Sim[12,1+`i']*dummy_educ4 		+ ///
							xsi2000Sim[13,1+`i']*age         		+ ///
							xsi2000Sim[14,1+`i']*age_2       		+ ///
							xsi2000Sim[15,1+`i']*exper1      		+ ///
							xsi2000Sim[16,1+`i']*exper2      		+ ///
							xsi2000Sim[17,1+`i']*exper3      		+ ///
							xsi2000Sim[18,1+`i']*exper4      		+ ///
							xsi2000Sim[19,1+`i']*exper5      		+ ///
							xsi2000Sim[20,1+`i']*exper6      		+ ///
							xsi2000Sim[21,1+`i']*exper7      

}

scatter Pers2000_0Data Pers2000_0Sim, title("Residual Sector") ytitle("Best Linear Predictor Data") ///
xtitle("Best Linear Predictor Model") legend(off) msymbol(p) color(blue) || ///
lfit Pers2000_0Data Pers2000_0Data, range(-0.2 1.2) lcolor(red) lwidth(medium) saving(Pers2000_0, replace)

scatter Pers2000_1Data Pers2000_1Sim, title("Agriculture/Mining") ytitle("Best Linear Predictor Data") ///
xtitle("Best Linear Predictor Model") legend(off) msymbol(p) color(blue) || ///
lfit Pers2000_1Data Pers2000_1Data, range(-0.2 1.2) lcolor(red) lwidth(medium) saving(Pers2000_1, replace)

scatter Pers2000_2Data Pers2000_2Sim, title("LT Manuf") ytitle("Best Linear Predictor Data") ///
xtitle("Best Linear Predictor Model") legend(off) msymbol(p) color(blue) || ///
lfit Pers2000_2Data Pers2000_2Data, range(-0.2 1.2) lcolor(red) lwidth(medium) saving(Pers2000_2, replace)

scatter Pers2000_3Data Pers2000_3Sim, title("HT Manuf") ytitle("Best Linear Predictor Data") ///
xtitle("Best Linear Predictor Model") legend(off) msymbol(p) color(blue) || ///
lfit Pers2000_3Data Pers2000_3Data, range(-0.2 1.2) lcolor(red) lwidth(medium) saving(Pers2000_3, replace)

scatter Pers2000_4Data Pers2000_4Sim, title("Construction") ytitle("Best Linear Predictor Data") ///
xtitle("Best Linear Predictor Model") legend(off) msymbol(p) color(blue) || ///
lfit Pers2000_4Data Pers2000_4Data, range(-0.2 1.2) lcolor(red) lwidth(medium) saving(Pers2000_4, replace)

scatter Pers2000_5Data Pers2000_5Sim, title("Trade") ytitle("Best Linear Predictor Data") ///
xtitle("Best Linear Predictor Model") legend(off) msymbol(p) color(blue) || ///
lfit Pers2000_5Data Pers2000_5Data, range(-0.2 1.2) lcolor(red) lwidth(medium) saving(Pers2000_5, replace)

scatter Pers2000_6Data Pers2000_6Sim, title("Trans/Util") ytitle("Best Linear Predictor Data") ///
xtitle("Best Linear Predictor Model") legend(off) msymbol(p) color(blue) || ///
lfit Pers2000_6Data Pers2000_6Data, range(-0.2 1.2) lcolor(red) lwidth(medium) saving(Pers2000_6, replace)

scatter Pers2000_7Data Pers2000_7Sim, title("Service") ytitle("Best Linear Predictor Data") ///
xtitle("Best Linear Predictor Model") legend(off) msymbol(p) color(blue) || ///
lfit Pers2000_7Data Pers2000_7Data, range(-0.2 1.2) lcolor(red) lwidth(medium) saving(Pers2000_7, replace) 

graph combine Pers2000_0.gph Pers2000_1.gph Pers2000_2.gph Pers2000_3.gph Pers2000_4.gph Pers2000_5.gph Pers2000_6.gph Pers2000_7.gph, saving(Pers2000DataSim, replace)

restore			  



********************
* Persistence 2005
********************

preserve

keep if idade_new >= 25 & idade_new <= 60-10 & ano == 1995

forvalues i = 0(1)7 {

	gen Pers2005_`i'Data =  xsi2005Data[1,1+`i']*lag_sector0        + ///
							xsi2005Data[2,1+`i']*lag_sector1        + ///
							xsi2005Data[3,1+`i']*lag_sector2        + ///
							xsi2005Data[4,1+`i']*lag_sector3        + ///
							xsi2005Data[5,1+`i']*lag_sector4        + ///
							xsi2005Data[6,1+`i']*lag_sector5        + ///
							xsi2005Data[7,1+`i']*lag_sector6        + ///
							xsi2005Data[8,1+`i']*lag_sector7        + ///
							xsi2005Data[9,1+`i']*gender      		+ ///
							xsi2005Data[10,1+`i']*dummy_educ2 		+ ///
							xsi2005Data[11,1+`i']*dummy_educ3 		+ ///
							xsi2005Data[12,1+`i']*dummy_educ4 		+ ///
							xsi2005Data[13,1+`i']*age         		+ ///
							xsi2005Data[14,1+`i']*age_2       		+ ///
							xsi2005Data[15,1+`i']*exper1      		+ ///
							xsi2005Data[16,1+`i']*exper2      		+ ///
							xsi2005Data[17,1+`i']*exper3      		+ ///
							xsi2005Data[18,1+`i']*exper4      		+ ///
							xsi2005Data[19,1+`i']*exper5      		+ ///
							xsi2005Data[20,1+`i']*exper6      		+ ///
							xsi2005Data[21,1+`i']*exper7

	gen Pers2005_`i'Sim =   xsi2005Sim[1,1+`i']*lag_sector0         + ///
							xsi2005Sim[2,1+`i']*lag_sector1         + ///
							xsi2005Sim[3,1+`i']*lag_sector2         + ///
							xsi2005Sim[4,1+`i']*lag_sector3         + ///
							xsi2005Sim[5,1+`i']*lag_sector4         + ///
							xsi2005Sim[6,1+`i']*lag_sector5         + ///
							xsi2005Sim[7,1+`i']*lag_sector6         + ///
							xsi2005Sim[8,1+`i']*lag_sector7         + ///
							xsi2005Sim[9,1+`i']*gender       		+ ///
							xsi2005Sim[10,1+`i']*dummy_educ2 		+ ///
							xsi2005Sim[11,1+`i']*dummy_educ3 		+ ///
							xsi2005Sim[12,1+`i']*dummy_educ4 		+ ///
							xsi2005Sim[13,1+`i']*age         		+ ///
							xsi2005Sim[14,1+`i']*age_2       		+ ///
							xsi2005Sim[15,1+`i']*exper1      		+ ///
							xsi2005Sim[16,1+`i']*exper2      		+ ///
							xsi2005Sim[17,1+`i']*exper3      		+ ///
							xsi2005Sim[18,1+`i']*exper4      		+ ///
							xsi2005Sim[19,1+`i']*exper5      		+ ///
							xsi2005Sim[20,1+`i']*exper6      		+ ///
							xsi2005Sim[21,1+`i']*exper7      

}

scatter Pers2005_0Data Pers2005_0Sim, title("Residual Sector") ytitle("Best Linear Predictor Data") ///
xtitle("Best Linear Predictor Model") legend(off) msymbol(p) color(blue) || ///
lfit Pers2005_0Data Pers2005_0Data, range(-0.2 1.2) lcolor(red) lwidth(medium) saving(Pers2005_0, replace)

scatter Pers2005_1Data Pers2005_1Sim, title("Agriculture/Mining") ytitle("Best Linear Predictor Data") ///
xtitle("Best Linear Predictor Model") legend(off) msymbol(p) color(blue) || ///
lfit Pers2005_1Data Pers2005_1Data, range(-0.2 1.2) lcolor(red) lwidth(medium) saving(Pers2005_1, replace)

scatter Pers2005_2Data Pers2005_2Sim, title("LT Manuf") ytitle("Best Linear Predictor Data") ///
xtitle("Best Linear Predictor Model") legend(off) msymbol(p) color(blue) || ///
lfit Pers2005_2Data Pers2005_2Data, range(-0.2 1.2) lcolor(red) lwidth(medium) saving(Pers2005_2, replace)

scatter Pers2005_3Data Pers2005_3Sim, title("HT Manuf") ytitle("Best Linear Predictor Data") ///
xtitle("Best Linear Predictor Model") legend(off) msymbol(p) color(blue) || ///
lfit Pers2005_3Data Pers2005_3Data, range(-0.2 1.2) lcolor(red) lwidth(medium) saving(Pers2005_3, replace)

scatter Pers2005_4Data Pers2005_4Sim, title("Construction") ytitle("Best Linear Predictor Data") ///
xtitle("Best Linear Predictor Model") legend(off) msymbol(p) color(blue) || ///
lfit Pers2005_4Data Pers2005_4Data, range(-0.2 1.2) lcolor(red) lwidth(medium) saving(Pers2005_4, replace)

scatter Pers2005_5Data Pers2005_5Sim, title("Trade") ytitle("Best Linear Predictor Data") ///
xtitle("Best Linear Predictor Model") legend(off) msymbol(p) color(blue) || ///
lfit Pers2005_5Data Pers2005_5Data, range(-0.2 1.2) lcolor(red) lwidth(medium) saving(Pers2005_5, replace)

scatter Pers2005_6Data Pers2005_6Sim, title("Trans/Util") ytitle("Best Linear Predictor Data") ///
xtitle("Best Linear Predictor Model") legend(off) msymbol(p) color(blue) || ///
lfit Pers2005_6Data Pers2005_6Data, range(-0.2 1.2) lcolor(red) lwidth(medium) saving(Pers2005_6, replace)

scatter Pers2005_7Data Pers2005_7Sim, title("Service") ytitle("Best Linear Predictor Data") ///
xtitle("Best Linear Predictor Model") legend(off) msymbol(p) color(blue) || ///
lfit Pers2005_7Data Pers2005_7Data, range(-0.2 1.2) lcolor(red) lwidth(medium) saving(Pers2005_7, replace) 

graph combine Pers2005_0.gph Pers2005_1.gph Pers2005_2.gph Pers2005_3.gph Pers2005_4.gph Pers2005_5.gph Pers2005_6.gph Pers2005_7.gph, saving(Pers2005DataSim, replace)

restore		



***********
* Frequency
***********

preserve

keep if idade_new >= 25 & idade_new <= 60-10 & ano == 1995

forvalues i = 0(1)7 {

	gen Freq`i'Data =       etaData[1,1+`i']*lag_sector0        + ///
							etaData[2,1+`i']*lag_sector1        + ///
							etaData[3,1+`i']*lag_sector2        + ///
							etaData[4,1+`i']*lag_sector3        + ///
							etaData[5,1+`i']*lag_sector4        + ///
							etaData[6,1+`i']*lag_sector5        + ///
							etaData[7,1+`i']*lag_sector6        + ///
							etaData[8,1+`i']*lag_sector7        + ///
							etaData[9,1+`i']*gender      		+ ///
							etaData[10,1+`i']*dummy_educ2 		+ ///
							etaData[11,1+`i']*dummy_educ3 		+ ///
							etaData[12,1+`i']*dummy_educ4 		+ ///
							etaData[13,1+`i']*age         		+ ///
							etaData[14,1+`i']*age_2       		+ ///
							etaData[15,1+`i']*exper1      		+ ///
							etaData[16,1+`i']*exper2      		+ ///
							etaData[17,1+`i']*exper3      		+ ///
							etaData[18,1+`i']*exper4      		+ ///
							etaData[19,1+`i']*exper5      		+ ///
							etaData[20,1+`i']*exper6      		+ ///
							etaData[21,1+`i']*exper7

	gen Freq`i'Sim =        etaSim[1,1+`i']*lag_sector0         + ///
							etaSim[2,1+`i']*lag_sector1         + ///
							etaSim[3,1+`i']*lag_sector2         + ///
							etaSim[4,1+`i']*lag_sector3         + ///
							etaSim[5,1+`i']*lag_sector4         + ///
							etaSim[6,1+`i']*lag_sector5         + ///
							etaSim[7,1+`i']*lag_sector6         + ///
							etaSim[8,1+`i']*lag_sector7         + ///
							etaSim[9,1+`i']*gender       		+ ///
							etaSim[10,1+`i']*dummy_educ2 		+ ///
							etaSim[11,1+`i']*dummy_educ3 		+ ///
							etaSim[12,1+`i']*dummy_educ4 		+ ///
							etaSim[13,1+`i']*age         		+ ///
							etaSim[14,1+`i']*age_2       		+ ///
							etaSim[15,1+`i']*exper1      		+ ///
							etaSim[16,1+`i']*exper2      		+ ///
							etaSim[17,1+`i']*exper3      		+ ///
							etaSim[18,1+`i']*exper4      		+ ///
							etaSim[19,1+`i']*exper5      		+ ///
							etaSim[20,1+`i']*exper6      		+ ///
							etaSim[21,1+`i']*exper7      

}

scatter Freq0Data Freq0Sim, title("Residual Sector") ytitle("Best Linear Predictor Data") ///
xtitle("Best Linear Predictor Model") legend(off) msymbol(p) color(blue) || ///
lfit Freq0Data Freq0Data, range(-0.5 8) lcolor(red) lwidth(medium) saving(Freq0, replace)

scatter Freq1Data Freq1Sim, title("Agriculture/Mining") ytitle("Best Linear Predictor Data") ///
xtitle("Best Linear Predictor Model") legend(off) msymbol(p) color(blue) || ///
lfit Freq1Data Freq1Data, range(-0.5 8) lcolor(red) lwidth(medium) saving(Freq1, replace)

scatter Freq2Data Freq2Sim, title("LT Manuf") ytitle("Best Linear Predictor Data") ///
xtitle("Best Linear Predictor Model") legend(off) msymbol(p) color(blue) || ///
lfit Freq2Data Freq2Data, range(-0.5 8) lcolor(red) lwidth(medium) saving(Freq2, replace)

scatter Freq3Data Freq3Sim, title("HT Manuf") ytitle("Best Linear Predictor Data") ///
xtitle("Best Linear Predictor Model") legend(off) msymbol(p) color(blue) || ///
lfit Freq3Data Freq3Data, range(-0.5 8) lcolor(red) lwidth(medium) saving(Freq3, replace)

scatter Freq4Data Freq4Sim, title("Construction") ytitle("Best Linear Predictor Data") ///
xtitle("Best Linear Predictor Model") legend(off) msymbol(p) color(blue) || ///
lfit Freq4Data Freq4Data, range(-0.5 8) lcolor(red) lwidth(medium) saving(Freq4, replace)

scatter Freq5Data Freq5Sim, title("Trade") ytitle("Best Linear Predictor Data") ///
xtitle("Best Linear Predictor Model") legend(off) msymbol(p) color(blue) || ///
lfit Freq5Data Freq5Data, range(-0.5 8) lcolor(red) lwidth(medium) saving(Freq5, replace)

scatter Freq6Data Freq6Sim, title("Trans/Util") ytitle("Best Linear Predictor Data") ///
xtitle("Best Linear Predictor Model") legend(off) msymbol(p) color(blue) || ///
lfit Freq6Data Freq6Data, range(-0.5 8) lcolor(red) lwidth(medium) saving(Freq6, replace)

scatter Freq7Data Freq7Sim, title("Service") ytitle("Best Linear Predictor Data") ///
xtitle("Best Linear Predictor Model") legend(off) msymbol(p) color(blue) || ///
lfit Freq7Data Freq7Data, range(-0.5 8) lcolor(red) lwidth(medium) saving(Freq7, replace) 

graph combine Freq0.gph Freq1.gph Freq2.gph Freq3.gph Freq4.gph Freq5.gph Freq6.gph Freq7.gph, saving(FreqDataSim, replace)

restore


********
* Return
********

preserve

forvalues i = 0(1)7 {

	gen Return`i'Data =     rhoData[1,1+`i']*year95       + ///
							rhoData[2,1+`i']*year96       + ///
							rhoData[3,1+`i']*year97       + ///
							rhoData[4,1+`i']*year98       + ///
							rhoData[5,1+`i']*year99       + ///
							rhoData[6,1+`i']*year00       + ///
							rhoData[7,1+`i']*year01       + ///
							rhoData[8,1+`i']*year02       + ///
							rhoData[9,1+`i']*year03       + ///
							rhoData[10,1+`i']*year04      + ///
							rhoData[11,1+`i']*year05      + ///
							rhoData[12,1+`i']*gender      + ///
							rhoData[13,1+`i']*dummy_educ2 + ///
							rhoData[14,1+`i']*dummy_educ3 + ///
							rhoData[15,1+`i']*dummy_educ4 + ///
							rhoData[16,1+`i']*age         + ///
							rhoData[17,1+`i']*age_2       + ///
							rhoData[18,1+`i']*exper1      + ///
							rhoData[19,1+`i']*exper2      + ///
							rhoData[20,1+`i']*exper3      + ///
							rhoData[21,1+`i']*exper4      + ///
							rhoData[22,1+`i']*exper5      + ///
							rhoData[23,1+`i']*exper6      + ///
							rhoData[24,1+`i']*exper7

	gen Return`i'Sim =      rhoSim[1,1+`i']*year95       + ///
							rhoSim[2,1+`i']*year96       + ///
							rhoSim[3,1+`i']*year97       + ///
							rhoSim[4,1+`i']*year98       + ///
							rhoSim[5,1+`i']*year99       + ///
							rhoSim[6,1+`i']*year00       + ///
							rhoSim[7,1+`i']*year01       + ///
							rhoSim[8,1+`i']*year02       + ///
							rhoSim[9,1+`i']*year03       + ///
							rhoSim[10,1+`i']*year04      + ///
							rhoSim[11,1+`i']*year05      + ///
							rhoSim[12,1+`i']*gender      + ///
							rhoSim[13,1+`i']*dummy_educ2 + ///
							rhoSim[14,1+`i']*dummy_educ3 + ///
							rhoSim[15,1+`i']*dummy_educ4 + ///
							rhoSim[16,1+`i']*age         + ///
							rhoSim[17,1+`i']*age_2       + ///
							rhoSim[18,1+`i']*exper1      + ///
							rhoSim[19,1+`i']*exper2      + ///
							rhoSim[20,1+`i']*exper3      + ///
							rhoSim[21,1+`i']*exper4      + ///
							rhoSim[22,1+`i']*exper5      + ///
							rhoSim[23,1+`i']*exper6      + ///
							rhoSim[24,1+`i']*exper7   
							
		scatter Return`i'Data Return`i'Sim if lag2_sector == `i' & lag_sector ~= `i' & lag_sector ~= . , title("Return`i'") ytitle("Best Linear Predictor Data") xtitle("Best Linear Predictor Model") ///
        legend(off)  msymbol(p) color(blue) || ///
		lfit Return`i'Data Return`i'Data if lag2_sector == `i' & lag_sector ~= `i' & lag_sector ~= . , range(-0.2 1.2) lcolor(red) lwidth(medium) saving(Return`i', replace)

}		

graph combine Return0.gph Return1.gph Return2.gph Return3.gph Return4.gph Return5.gph Return6.gph Return7.gph, saving(ReturnDataSim, replace)
			  
****************************  
* Simple Moments in the Data			  
****************************  

clear
clear matrix
cd "C:\Users\Dix-Carneiro\Documents\Trade and Labor Dynamics\Revision_Econometrica\GoodnessOfFit"

use new_painel

sort pis

merge pis using choices_wide
			  
gen wage = exp(log_w)

matrix WageMeans = J(7,2,0)
		  
forvalues i = 1(1)7 {
	sum wage if sector == `i'
	matrix WageMeans[`i',1] = r(mean)
	sum log_w if sector == `i' 
	matrix WageMeans[`i',2] = r(mean)
}



tab sector, gen(emp)

matrix EmpMeans = J(8,1,0)

forvalues i = 1(1)8 {
	sum emp`i'
	matrix EmpMeans[`i',1] = 100.0*r(mean)
}


matrix Transition = J(8,8,0)

forvalues s1 = 1(1)8 {
	forvalues i = 1(1)8 {
		sum emp`i' if lag_sector == `s1'-1 
		matrix Transition[`s1',`i'] = 100.0*r(mean)
	}
}

matrix Return = J(8,1,0)

forvalues i = 1(1)8 {
	sum emp`i' if lag2_sector == `i'-1 & lag_sector ~= `i'-1 & lag_sector ~= .
	matrix Return[`i',1] = 100.0*r(mean)
}


matrix Freq = J(8,1,0)

forvalues i = 0(1)7 {
	sum freq`i' if idade_new >= 25 & idade_new <= 50 & ano == 1995
	matrix Freq[`i'+1,1] = r(mean)
}


clear
cd "C:\Users\Dix-Carneiro\Documents\Trade and Labor Dynamics\Revision_Econometrica\GoodnessOfFit\TablesGraphs"

svmat Return

outsheet Return* using ReturnData.csv, comma replace

svmat Freq

outsheet Freq* using FreqData.csv, comma replace

svmat WageMeans 

outsheet WageMeans* using WageMeansData.csv, comma replace

svmat EmpMeans 

outsheet EmpMeans* using EmpMeansData.csv, comma replace

svmat Transition

outsheet Transition* using TransitionData.csv, comma replace



**************************
* Simple Simulated Moments			  
**************************

clear
clear matrix
cd "C:\Users\Dix-Carneiro\Documents\Trade and Labor Dynamics\Revision_Econometrica\GoodnessOfFit"

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


keep id coh year choice

egen gp = group(id coh)

reshape wide choice, i(gp) j(year)

sort gp

save choices_wide, replace


clear

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

egen gp = group(id coh)

sort gp

merge gp using choices_wide

gen freq0 = 0
gen freq1 = 0
gen freq2 = 0
gen freq3 = 0
gen freq4 = 0
gen freq5 = 0
gen freq6 = 0
gen freq7 = 0

forvalues i = 1995(1)2005{
	replace freq0 = freq0 + 1 if choice == 0 
	replace freq1 = freq1 + 1 if choice == 1
	replace freq2 = freq2 + 1 if choice == 2
	replace freq3 = freq3 + 1 if choice == 3
	replace freq4 = freq4 + 1 if choice == 4
	replace freq5 = freq5 + 1 if choice == 5
	replace freq6 = freq6 + 1 if choice == 6
	replace freq7 = freq7 + 1 if choice == 7
}

matrix FreqSim = J(8,1,0)

forvalues i = 0(1)7 {
	sum freq`i' if year == 1995 & year - coh >= 25 & year - coh <= 50 [iw = wgt]
	matrix FreqSim[`i'+1,1] = r(mean)
}


*****************************************************************************************
*****************************************************************************************

tab choice [iw = wgt], gen(emp)

matrix ReturnSim = J(8,1,0)

forvalues i = 1(1)8 {
	sum emp`i' if choice2 == `i'-1 & choice1 ~= `i'-1 & choice1 ~= . [iw = wgt]
	matrix ReturnSim[`i',1] = 100.0*r(mean)
}

*****************************************************************************************
*****************************************************************************************

gen log_w = log(wage)


gen age = year - coh


keep if age >= 25 & age <= 60


replace age = year - coh - 25

matrix WageMeans = J(7,2,0)

forvalues i = 1(1)7 {
	sum wage if choice == `i' [iw = wgt]
	matrix WageMeans[`i',1] = r(mean)
	sum log_w if choice == `i' [iw = wgt]
	matrix WageMeans[`i',2] = r(mean)
}

*****************************************************************************************
*****************************************************************************************

matrix EmpMeans = J(8,1,0)

forvalues i = 1(1)8 {
	sum emp`i' [iw = wgt]
	matrix EmpMeans[`i',1] = 100.0*r(mean)
}

*****************************************************************************************
*****************************************************************************************

matrix Transition = J(8,8,0)

forvalues s1 = 1(1)8 {
	forvalues i = 1(1)8 {
		sum emp`i' if choice1 == `s1'-1 [iw = wgt]
		matrix Transition[`s1',`i'] = 100.0*r(mean)
	}
}

clear
cd "C:\Users\Dix-Carneiro\Documents\Trade and Labor Dynamics\Revision_Econometrica\GoodnessOfFit\TablesGraphs"

svmat ReturnSim

outsheet ReturnSim* using ReturnSim.csv, comma replace

svmat FreqSim

outsheet FreqSim* using FreqSim.csv, comma replace

svmat WageMeans 

outsheet WageMeans* using WageMeansSim.csv, comma replace

svmat EmpMeans 

outsheet EmpMeans* using EmpMeansSim.csv, comma replace

svmat Transition

outsheet Transition* using TransitionSim.csv, comma replace
