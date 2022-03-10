/*
For paper version sent to AER
*/
/*
global doFiles "/Users/sharontraiberman/Dropbox/DanishLEEDProject/structuralEstimation_Final/Estimation/doFiles"
global dataPath "/Users/sharontraiberman/Dropbox/DanishLEEDProject/structuralEstimation_Final/Estimation/inputs"
global outPath "/Users/sharontraiberman/Dropbox/DanishLEEDProject/structuralEstimation_Final/Estimation/outputs"
global matlabPath "/Users/sharontraiberman/Dropbox/DanishLEEDProject/structuralEstimation_Final/Estimation/matlabFiles"

global figurePath "/Users/sharontraiberman/Dropbox/DanishLEEDProject/structuralEstimation_Final/Estimation/figures"
global tablePath "/Users/sharontraiberman/Dropbox/DanishLEEDProject/structuralEstimation_Final/Estimation/tables"
global workerFiles /Users/sharontraiberman/Dropbox/DanishLEEDProject/structuralEstimation_Final/Calibration/workerFiles

global calibrationPathWorkers /Users/sharontraiberman/Dropbox/DanishLEEDProject/structuralEstimation_Final/Calibration/workerFiles
*/

global doFiles "/Users/sharont/Dropbox/DanishLEEDProject/structuralEstimation_Final/Estimation/doFiles"
global dataPath "/Users/sharont/Dropbox/DanishLEEDProject/structuralEstimation_Final/Estimation/inputs"
global outPath "/Users/sharont/Dropbox/DanishLEEDProject/structuralEstimation_Final/Estimation/outputs"
global matlabPath "/Users/sharont/Dropbox/DanishLEEDProject/structuralEstimation_Final/Estimation/matlabFiles"

* ~/Dropbox/DanishLEEDProject/structuralEstimation_Final/Estimation/matlabFiles/

global figurePath "/Users/sharont/Dropbox/DanishLEEDProject/structuralEstimation_Final/Estimation/figures"
global tablePath "/Users/sharont/Dropbox/DanishLEEDProject/structuralEstimation_Final/Estimation/tables"
global workerFiles /Users/sharont/Dropbox/DanishLEEDProject/structuralEstimation_Final/Calibration/workerFiles

global calibrationPathWorkers /Users/sharont/Dropbox/DanishLEEDProject/structuralEstimation_Final/Calibration/workerFiles
global statsDKoutput /Users/sharont/Dropbox/DanishLEEDProject/structuralEstimation_Final/Estimation/statdk_tables


cd $doFiles

/*
/*******
Make EDF
*******/
*% State = (age, type, tenure, occ)

forval year = 1996(1)2007{
	insheet using $dataPath/CENSORED_EDF_`year'.csv, c clear
	reshape long q, i(year skills employed newcode tenure smoothedage) j(type)
	rename q freq

	egen group = group(skills type)
	drop skills type

	gen age = smoothedage*3+23
	save temp, replace
	gen correct = 1
	append using temp
	replace age = age+1 if correct==.
	replace correct = 1
	append using temp
	replace age = age+2 if correct==.
	replace correct = 1

	* 59 got grouped into its own group so need to correct
	drop if age>=60

	replace freq = freq/3 if age<59

	* Matlab Model Assumes Age>25 [Can we fix this]?
	drop if age<25

	egen total = total(freq)
	gen double share = freq/total
	drop total freq smoothedage correct

	order year age group tenure newcode employed
	save tempEDF, replace

	* Make a grid for potential 0s...
	clear
	set ob 38
	gen newcode = _n
	save temp, replace

	gen age = 25
	forval j = 26(1)59{
		append using temp
		replace age = `j' if age==.
	}
	save temp, replace

	gen tenure = 0
	forval j = 1(1)5{
		append using temp
		replace tenure = `j' if tenure==.
	}
	save temp, replace

	gen group = 1
	forval j = 2(1)6{
		append using temp
		replace group = `j' if group==.
	}

	save temp, replace
	gen employed = 1
	append using temp
	replace employed = 0 if employed==.

	merge 1:1 age group tenure newcode employed using tempEDF
	drop _merge
	replace year = 1996 
	replace share = 0 if share==.

	* Population Normalization
	replace share = share*1000

	* Normalize tenure to be end-of-period
	replace tenure = tenure+1

	sort employed year age group tenure newcode
	order employed year age group tenure newcode

	
	outsheet year age group tenure newcode share if employed==1 using $calibrationPathWorkers/workerEDF`year'.csv, replace
	outsheet year age group tenure newcode share if employed==0 using $calibrationPathWorkers/UworkerEDF`year'.csv, replace
}

*/


/***************************************
	UNPACK PARAMETERS FOR CALIBRATION
***************************************/

/************
	COSTS
*************/
import delimited $statsDKoutput/XQ_BOOTS.txt, clear
* Original estimate
keep v1
gen param = _n

* RHO
matrix RHO = 1/v1[1]

* PRODUCTIVITY FUNCTION
matrix FTYPE = J(1,6,0)
forval j = 1(1)5{
	local ind = `j'+1
	matrix FTYPE[1,`j'] = v[`ind']
}

matrix FAGE = v[7]
matrix FAGE2 = v[8]

* GAMMA_CONS
matrix GAMMA0 = v[9]

* GAMMA_VEC
matrix GAMMA_VEC = J(1,22,.)
forval j = 1(1)22{
	local ind = `j'+9
	matrix GAMMA_VEC[1,`j'] = v[`ind']
}

save $outPath/xQ, replace
	
** CREATING OUTPUT FOR MATLAB ** 
* 1. RHO
clear
svmat RHO
outsheet using $matlabPath/RHOQ.csv, c replace nol non noq

* 2. GAMMA
matrix GAMMAQ = [GAMMA0, 0, GAMMA_VEC]'
clear
svmat GAMMAQ
outsheet using $matlabPath/GAMMAQ.csv, c replace nol non noq

* 3. F(OMEGA)
matrix SEARCHQ = [FAGE, FAGE2, FTYPE]'
clear
svmat SEARCHQ
outsheet using $matlabPath/SEARCHQ.csv, c replace nol non noq

* 4. U(OMEGA) val32-val38
use $outPath/xQ, clear
keep if param>=32 & param<=38
drop param
set obs 8
replace v1 = 0 if v1==.
outsheet using $matlabPath/NONEMPQ.csv, c replace nol non noq

* 5. ETA val39 - val76
use $outPath/xQ, clear
keep if param>=39
drop param
outsheet using $matlabPath/ETAQ.csv, c replace nol non noq



* 6. SKILL PRICES

import delimited $statsDKoutput/BETAYEARS.txt, clear
keep if sample==0
drop sample
reshape long betay, i(newcode) j(period)
save $outPath\BETAYEARS0, replace





** CREATING TABLES ** 
use $outPath/occDistMatrix, clear

sort newCodeA newCodeB 
rename upDiff1-downDiff10 char#,addnumber(1)
gen char21 = broadSectorA!= broadSectorB
gen char22 = disco2A!= disco2B

gen costs = GAMMA0[1,1]
forval j = 1/22{
	local k = `j'+2
	replace costs = costs + char`j'*GAMMA_VEC[1,`j']
}

keep newCodeA newCodeB costs
save $outPath/costMatrix, replace



use $outPath/costMatrix, clear

rename newCodeA newcode
merge m:1 newcode using $workerFiles/newcodenames
drop _merge
rename newcode newCodeA
rename (disco2 broadSector) (disco2A broadSectorA)

rename newCodeB newcode
merge m:1 newcode using $workerFiles/newcodenames
drop _merge
rename newcode newCodeB
rename (disco2 broadSector) (disco2B broadSectorB)

drop if newCodeA==newCodeB

replace costs = exp(costs)*RHO[1,1]



preserve
collapse (mean) costs, by(broadSectorA)
gen broadSectorB = 5
save meanExitCosts, replace
restore

preserve
collapse (mean) costs, by(broadSectorB)
gen broadSectorA = 5
save meanEntryCosts, replace
restore

preserve
drop if broadSectorA==broadSectorB
collapse (mean) costs, by(broadSectorA)
gen broadSectorB = 6
save meanExitCostsOffDiagonal, replace
restore

preserve
drop if broadSectorA==broadSectorB
collapse (mean) costs, by(broadSectorB)
gen broadSectorA = 6
save meanEntryCostsOffDiagonal, replace
restore

append using meanExitCosts
append using meanEntryCosts

append using meanExitCostsOffDiagonal
append using meanEntryCostsOffDiagonal


collapse (mean) costs, by(broadSectorA broadSectorB)
reshape wide costs, i(broadSectorA) j(broadSectorB)


format costs* %6.2f
tostring costs*, replace u force

forval j = 1(1)6{
	replace costs`j' = "" if costs`j'=="."
}

cap label drop costMatName
label define costMatName 1 "Manufacturing",modify
label define costMatName 2 "Services",modify
label define costMatName 3 "FIRE",modify
label define costMatName 4 "H \& E",modify
label define costMatName 5 "Avg. Entry (w/ Diag.)",modify
label define costMatName 6 "Avg. Entry (w/o Diag.)",modify

label values broadSectorA costMatName 

decode broadSectorA, gen(broadSectorStr)
replace broadSectorStr = trim(broadSectorStr)

gen amp = "&"
gen lb = "\\"

order broadSectorStr costs*

outsheet broadSectorStr amp costs1 amp costs2 amp costs3 amp costs4 amp costs5 amp costs6 lb using $tablePath/costMatrix.tex, replace non nol noq



**********************
** Parameter Tables **
**********************
* Wage Parameters
import delimited $statsDKoutput/BETAMAT.txt, clear
cap rename (v1-v11) (sample newcode beta1 beta2 beta3 beta4 beta5 beta6 beta7 beta8 sigma)
rename sigma beta9

* Adjusting Age Parameters (just for reading):
* Coef on a to Coef on (a-23)
* Multiply coef on a^2*1000 for readability
replace beta2 = beta2+23*2*beta3
replace beta3 = beta3*1000

forval j = 1(1)9{
	egen stdErr`j' = sd(beta`j'),by(newcode)
}

keep if sample==0
drop sample

format beta* stdErr* %6.3f
tostring beta*, replace u force
tostring stdErr*, replace u force

forval j = 1(1)9{
	replace beta`j' = "" if beta`j'=="0.000"
	replace stdErr`j' = "" if beta`j'==""
	
	*replace beta`j' = "$" + beta`j' + "$"
	replace stdErr`j' = "("+stdErr`j'+")" if stdErr`j'!=""
}


rename beta# v#1
rename stdErr# v#2

reshape long v1 v2 v3 v4 v5 v6 v7 v8 v9, i(newcode) j(row)

rename v* beta*

merge m:1 newcode using $workerFiles/newcodenames
drop _merge

cap label drop discoNames
label define discoNames 12 "Managers",add
label define discoNames 21 "Science Prof.",modify
label define discoNames 22 "Health Prof.",modify
label define discoNames 23 "Teachers",modify
label define discoNames 24 "Other Prof.",modify
label define discoNames 31 "Science Assc. Prof.",modify
label define discoNames 32 "Health Assc. Prof.",modify
label define discoNames 33 "Teaching Assc. Prof.",modify
label define discoNames 34 "Other Assc. Prof.",modify
label define discoNames 41 "Clerks",modify
label define discoNames 42 "Customer Service",modify
label define discoNames 51 "Personal Services",modify
label define discoNames 52 "Retail Occs.",modify
label define discoNames 61 "Agriculture ",modify
label define discoNames 71 "Building Trades",modify
label define discoNames 72 "Metal Trades",modify
label define discoNames 74 "Other Crafts",modify
label define discoNames 81 "Plant Operator",modify
label define discoNames 82 "Machine Operator",modify
label define discoNames 83 "Drivers",modify
label define discoNames 91 "Elementary Occs.",modify
label define discoNames 93 "Laborers",modify

label values disco2 discoNames 

decode disco2, gen(discoStr)
replace discoStr = "" if row==2

sort broadSector disco2 row

gen amp = "&"
gen lb = "\\"

outsheet discoStr amp beta1 amp beta2 amp beta3 amp beta4 amp beta5 amp beta6 amp beta7 amp beta8 amp beta9 lb using $tablePath/wageParameters_MAN.tex if broadSector==1, noq replace non 
outsheet discoStr amp beta1 amp beta2 amp beta3 amp beta4 amp beta5 amp beta6 amp beta7 amp beta8 amp beta9 lb using $tablePath/wageParameters_SERV.tex if broadSector==2, noq replace non 
outsheet discoStr amp beta1 amp beta2 amp beta3 amp beta4 amp beta5 amp beta6 amp beta7 amp beta8 amp beta9 lb using $tablePath/wageParameters_FIRE.tex if broadSector==3, noq replace non 
outsheet discoStr amp beta1 amp beta2 amp beta3 amp beta4 amp beta5 amp beta6 amp beta7 amp beta8 amp beta9 lb using $tablePath/wageParameters_HE.tex if broadSector==4, noq replace non 



import delimited $statsDKoutput/BETAYEARS.txt, clear
*reporting the time mean & putting in levels (not logs)
egen meanPrice = rowmean(betay1-betay12)
replace meanPrice = exp(meanPrice)

keep newcode sample meanPrice


egen stdErr = sd(meanPrice),by(newcode)
keep if sample==0
drop sample

format meanPrice stdErr %6.3f
tostring meanPrice, replace u force
*replace meanPrice = "$" + meanPrice + "$"

tostring stdErr, replace u force
replace stdErr = "("+stdErr+")"


rename (meanPrice stdErr) (v1 v2)
reshape long v, i(newcode) j(row)



merge m:1 newcode using $workerFiles/newcodenames
drop _merge

drop newcode
reshape wide v , i(disco2 row) j(broadSector)

cap label drop discoNames
label define discoNames 12 "Managers",add
label define discoNames 21 "Science Professional",modify
label define discoNames 22 "Health Professional",modify
label define discoNames 23 "Teachers",modify
label define discoNames 24 "Other Professional",modify
label define discoNames 31 "Science Assc. Professional",modify
label define discoNames 32 "Health Assc. Professional",modify
label define discoNames 33 "Teaching Assc. Professional",modify
label define discoNames 34 "Other Assc. Professional",modify
label define discoNames 41 "Clerks",modify
label define discoNames 42 "Customer Service",modify
label define discoNames 51 "Personal Workers",modify
label define discoNames 52 "Retail Workers",modify
label define discoNames 61 "Agriculture ",modify
label define discoNames 71 "Building Trades",modify
label define discoNames 72 "Metal Trades",modify
label define discoNames 74 "Other Crafts",modify
label define discoNames 81 "Plant Operator",modify
label define discoNames 82 "Machine Operator",modify
label define discoNames 83 "Drivers",modify
label define discoNames 91 "Elementary Occupations",modify
label define discoNames 93 "Laborers",modify
label values disco2 discoNames 

decode disco2, gen(discoStr)
replace discoStr = "" if row==2

gen missing = v1==""
replace missing = missing+1 if v2=="" & missing>0
replace missing = missing+1 if v3=="" & missing>1
replace missing = missing+1 if v4=="" & missing>2
sort missing disco2 row


gen amp = "&"
gen lb = "\\"

outsheet discoStr amp v1 amp v2 amp v3 amp v4 lb using $tablePath/skillPriceParameters.tex, noq replace non

***************************
**      CHI TESTS        **
***************************
matrix Ftests = [0, 0, 0, 0, 0 ,0, 0, 0]
import delimited $statsDKoutput/BETAMAT.txt, clear
cap rename (v1-v11) (sample newcode beta1 beta2 beta3 beta4 beta5 beta6 beta7 beta8 sigma)

* Need to skip 3 in the new version -- can't get enough decimal points
* instead, do this ON the server (still rejected)


forval j = 1(1)2{
	display `j'
	preserve
	keep newcode beta`j' sample
	rename beta`j' beta
	drop if beta==0
	
	egen betaMean = mean(beta), by(sample)
	gen thetaHat = beta-betaMean
	drop beta betaMean
	
	reshape wide thetaHat, i(sample) j(newcode)
	corr thetaHat*, cov
	
	matrix A = r(C)
	mean thetaHat* 
	*if sample==0
	matrix thetaHat = e(b)
	matrix Ftests[1,`j'] = thetaHat*inv(A)*thetaHat'
	restore
}

forval j = 4(1)8{
	display `j'
	preserve
	keep newcode beta`j' sample
	rename beta`j' beta
	drop if beta==0
	
	egen betaMean = mean(beta), by(sample)
	gen thetaHat = beta-betaMean
	drop beta betaMean
	
	reshape wide thetaHat, i(sample) j(newcode)
	corr thetaHat*, cov
	
	matrix A = r(C)
	mean thetaHat* 
	*if sample==0
	matrix thetaHat = e(b)
	matrix Ftests[1,`j'] = thetaHat*inv(A)*thetaHat'
	restore
}


***************************
** Elasticity Parameters **
***************************
import delimited $statsDKoutput/XQ_BOOTS.txt, clear
keep if _n==1
gen param = 1
reshape long v, i(param) j(boot)
replace boot = boot-1

egen stdErr = sd(v)
keep if boot==0
drop boot
drop param
gen spec = 1
save temp, replace


import delimited $statsDKoutput/XQ_BOOTS_TIME.txt, clear
keep if _n==1
gen param = 1
reshape long v, i(param) j(boot)
replace boot = boot-1

egen stdErrTime = sd(v)
keep if boot==0
drop boot
drop param
gen spec = 1
merge 1:1 spec using temp
drop _merge
order spec v stdErr stdErrTime

save temp, replace

import delimited $statsDKoutput/aggBoots.txt, clear
rename v2 boot
rename v1 v2
rename v3 v3

order boot v2 v3
egen stdErrTime2 = sd(v2)
egen stdErrTime3 = sd(v3)
keep if boot==0
reshape long stdErrTime v, i(boot) j(spec)
drop boot

append using temp
order spec v stdErr stdErrTime
save temp, replace
sort spec
duplicates drop

gen id = 1

format v stdErr* %6.3f
tostring v, replace u force

tostring stdErr, replace u force
replace stdErr = "("+stdErr+")"
replace stdErr = "" if stdErr=="(.)"

tostring stdErrTime, replace u force
replace stdErrTime = "["+stdErrTime+"]"

rename v bv
rename stdErr bstdErr
rename stdErrTime bstdErrTime
reshape wide bv bstdErr bstdErrTime, i(id) j(spec)


rename b*1 b1*
rename b*2 b2*
rename b*3 b3*
reshape long b1 b2 b3, i(id) j(stat) string

gen order = 1 if stat=="v"
replace order = 2 if stat=="stdErr"
replace order = 3 if stat==""
sort order

replace id = _n

replace stat = "" if stat!="v"
replace stat = "\$1/\rho$" if stat == "v"
gen amp = "&"
gen lb = "\\"

outsheet stat amp b1 amp b2 amp b3 lb using $tablePath/rhoAcrossSpecs.tex, noq replace non 


save rhoTemp, replace


***** Brings in Average C/rho *****
clear
set obs 0
gen cost = .
save meanGamma, replace

import delimited $statsDKoutput/XQ_BOOTS.txt, clear
gen param = _n
rename v* val#, addnumber(0)
reshape long val, i(param) j(boot)
reshape wide val, i(boot) j(param)

keep boot val9 - val31

rename val# gammaVal#, addnumber(0)

reshape long gammaVal, i(boot) j(char)

forval bootNum = 0(1)100{
	display "Starting Iteration `bootNum'"
	preserve
	keep if boot == `bootNum'
	save gammaTemp, replace

	use $outPath/occDistMatrix, clear

	sort newCodeA newCodeB 
	rename upDiff1-downDiff10 charVal#,addnumber(1)
	gen charVal21 = broadSectorA!= broadSectorB
	gen charVal22 = disco2A!= disco2B 
	gen charVal0 = 1
	
	reshape long charVal, i(newCodeA newCodeB) j(char)
	merge m:1 char using gammaTemp
	drop _merge
	gen cost = charVal*gammaVal
	collapse (sum) cost, by(newCodeA newCodeB boot)
	replace cost = exp(cost)
	
	collapse (mean) cost, by(boot)
	
	append using meanGamma
	save meanGamma, replace
	restore
}




clear
set obs 0
gen cost = .
save meanGamma_TIME, replace

import delimited $statsDKoutput/XQ_BOOTS_TIME.txt, clear
gen param = _n
rename v* val#, addnumber(0)
reshape long val, i(param) j(boot)
reshape wide val, i(boot) j(param)

keep boot val9 - val31

rename val# gammaVal#, addnumber(0)

reshape long gammaVal, i(boot) j(char)

forval bootNum = 0(1)100{
	display "Starting Iteration `bootNum'"
	preserve
	keep if boot == `bootNum'
	save gammaTemp, replace

	use $outPath/occDistMatrix, clear

	sort newCodeA newCodeB 
	rename upDiff1-downDiff10 charVal#,addnumber(1)
	gen charVal21 = broadSectorA!= broadSectorB
	gen charVal22 = disco2A!= disco2B 
	gen charVal0 = 1
	
	reshape long charVal, i(newCodeA newCodeB) j(char)
	merge m:1 char using gammaTemp
	drop _merge
	gen cost = charVal*gammaVal
	collapse (sum) cost, by(newCodeA newCodeB boot)
	replace cost = exp(cost)
	
	collapse (mean) cost, by(boot)
	
	append using meanGamma_TIME
	save meanGamma_TIME, replace
	restore
}


use meanGamma, clear
rename cost cost1
merge 1:1 boot using meanGamma_TIME
drop _merge
rename cost cost2

forval j = 1(1)2{
	egen stdErr`j' = sd(cost`j')
}
keep if boot == 0
drop boot
gen spec = 1
gen id = 1

drop cost2
rename cost1 cost
rename stdErr1 stdErr
rename stdErr2 stdErrTime

format cost stdErr* %6.3f
tostring cost, replace u force
tostring stdErr, replace u force
replace stdErr = "("+stdErr+")"

tostring stdErrTime, replace u force
replace stdErrTime = "["+stdErrTime+"]"


rename cost bcost
rename stdErr bstdErr
rename stdErrTime bstdErrTime
reshape wide bcost bstdErr bstdErrTime, i(id) j(spec)

rename b*1 b1*
reshape long b1, i(id) j(stat) string

* Temporary
gen b2 = ""
gen b3 = ""

replace stat = "" if stat!="cost"
replace stat = "Mean \$C/\rho$" if stat == "cost"
gen amp = "&"
gen lb = "\\"

*outsheet stat amp b1 amp b2 amp b3 lb using $tablePath/meanCostsAcrossSpecs.tex, noq replace non 



*********************
** Cost Parameters **
*********************
import delimited $statsDKoutput/XQ_BOOTS_TIME.txt, clear
gen param = _n
rename v* val#, addnumber(0)
reshape long val, i(param) j(boot)
reshape wide val, i(boot) j(param)

keep boot val9 - val31
rename val# val#, addnumber(0)
reshape long val, i(boot) j(Gamma)

egen stdErrTime = sd(val),by(Gamma)
keep if boot==0
save temp, replace

import delimited $statsDKoutput/XQ_BOOTS.txt, clear
gen param = _n
rename v* val#, addnumber(0)
reshape long val, i(param) j(boot)
reshape wide val, i(boot) j(param)

keep boot val9 - val31
rename val# val#, addnumber(0)
reshape long val, i(boot) j(Gamma)

egen stdErr = sd(val),by(Gamma)
keep if boot==0

merge 1:1 boot Gamma using temp
drop _merge

format val stdErr* %6.3f
tostring val, replace u force
*replace val = "$" + val + "$"
tostring stdErr, replace u force
replace stdErr = "("+stdErr+")"

tostring stdErrTime, replace u force
replace stdErrTime = "["+stdErrTime+"]"

rename (val stdErr stdErrTime) (v1 v2 v3)
reshape long v, i(boot Gamma) j(row)
drop boot
gen down = 1 if Gamma>=11 & Gamma<=20
replace down = 0 if down==.
replace Gamma = Gamma-10 if down==1
reshape wide v, i(Gamma row) j(down)

replace Gamma = -2 if Gamma==0
replace Gamma = -1 if Gamma==21
replace Gamma = 0 if Gamma==22

sort Gamma row
replace Gamma = . if row==2|row==3
tostring Gamma, replace

gen gammaString = "Task " + Gamma
replace gammaString = "Constant" if Gamma == "-2"
replace gammaString = "Occ. Dummy" if Gamma=="-1"
replace gammaString = "Sec. Dummy" if Gamma=="0"
replace gammaString = "" if row==2|row==3

drop Gamma row

gen amp = "&"
gen lb = "\\"

outsheet gammaStr amp v0 amp v1 lb using $tablePath/costParameters.tex, noq replace non 





**************************************
** Mobility Productivity Parameters **
**************************************
import delimited $statsDKoutput/XQ_BOOTS_TIME.txt, clear
gen param = _n
rename v* val#, addnumber(0)
reshape long val, i(param) j(boot)
reshape wide val, i(boot) j(param)

keep boot val2 - val8

* rearranging to put age out front
rename val# val#, addnumber(3)
rename val8 val1
rename val9 val2

* Make squared term readable
replace val2 = val2*1000

reshape long val, i(boot) j(Gamma)

egen stdErrTime = sd(val),by(Gamma)
keep if boot==0
drop boot val
save temp, replace

import delimited $statsDKoutput/XQ_BOOTS.txt, clear
gen param = _n
rename v* val#, addnumber(0)
reshape long val, i(param) j(boot)
reshape wide val, i(boot) j(param)

keep boot val2 - val8

* rearranging to put age out front
rename val# val#, addnumber(3)
rename val8 val1
rename val9 val2

* Make squared term readable
replace val2 = val2*1000

reshape long val, i(boot) j(Gamma)

egen stdErr = sd(val),by(Gamma)
keep if boot==0
drop boot
merge 1:1 Gamma using temp
drop _merge

* adjusting Age^2 up to make it readable
*replace val = val*1000 if Gamma==2
*åreplace stdErr = stdErr*1000 if Gamma==2

format val stdErr* %6.3f
tostring val, replace u force
*replace val = "$" + val + "$"
tostring stdErr, replace u force
replace stdErr = "("+stdErr+")"

tostring stdErrTime, replace u force
replace stdErrTime = "["+stdErrTime+"]"

rename (val stdErr stdErrTime) (v1 v2 v3)
reshape long v, i(Gamma) j(row)

cap label drop costProdParam
label define costProdParam 1 "Age", add 
label define costProdParam 2 "$\text{Age}^2(\times 1000)$", modify 
label define costProdParam 3 "Type 1", modify
label define costProdParam 4 "Type 2", modify
label define costProdParam 5 "Type 3", modify
label define costProdParam 6 "Type 4", modify
label define costProdParam 7 "Type 5", modify
label values Gamma costProdParam 

decode Gamma, gen(gammaStr)
replace gammaStr = "" if row==2|row==3

gen amp = "&"
gen lb = "\\"

gen famp = "\hspace{10mm} &"

gen bamp = "& \hspace{10mm}"

outsheet famp gammaStr amp v bamp lb using $tablePath/costProductivityParameters.tex, noq replace non


********************
** Eta Parameters **
********************
import delimited $statsDKoutput/XQ_BOOTS_TIME.txt, clear
gen param = _n
rename v* val#, addnumber(0)
reshape long val, i(param) j(boot)

* Adjusts for variance
egen rho = max(val*(param==1)),by(boot)
replace val = val/rho
drop rho

reshape wide val, i(boot) j(param)

keep boot val39 - val76
rename val# eta#, addnumber(1)
reshape long eta, i(boot) j(newcode)

egen stdErrTime = sd(eta),by(newcode)
keep if boot==0
drop boot eta

save $tablePath/etaVec, replace



import delimited $statsDKoutput/XQ_BOOTS.txt, clear
gen param = _n
rename v* val#, addnumber(0)
reshape long val, i(param) j(boot)

* Adjusts for variance
egen rho = max(val*(param==1)),by(boot)
replace val = val/rho
drop rho

reshape wide val, i(boot) j(param)

keep boot val39 - val76
rename val# eta#, addnumber(1)
reshape long eta, i(boot) j(newcode)

egen stdErr = sd(eta),by(newcode)
keep if boot==0
drop boot

merge 1:1 newcode using $tablePath/etaVec
drop _merge

save $tablePath/etaVec, replace

format eta stdErr* %6.3f
tostring eta, replace u force
*replace eta = "$" + eta + "$"
tostring stdErr, replace u force
replace stdErr = "("+stdErr+")"

tostring stdErrTime, replace u force
replace stdErrTime = "["+stdErrTime+"]"


rename (eta stdErr stdErrTime) (v1 v2 v3)
reshape long v, i(newcode) j(row)


merge m:1 newcode using $workerFiles/newcodenames
drop _merge

drop newcode
reshape wide v , i(disco2 row) j(broadSector)

cap label drop discoNames
label define discoNames 12 "Managers",add
label define discoNames 21 "Science Professional",modify
label define discoNames 22 "Health Professional",modify
label define discoNames 23 "Teachers",modify
label define discoNames 24 "Other Professional",modify
label define discoNames 31 "Science Assc. Professional",modify
label define discoNames 32 "Health Assc. Professional",modify
label define discoNames 33 "Teaching Assc. Professional",modify
label define discoNames 34 "Other Assc. Professional",modify
label define discoNames 41 "Clerks",modify
label define discoNames 42 "Customer Service",modify
label define discoNames 51 "Personal Workers",modify
label define discoNames 52 "Retail Workers",modify
label define discoNames 61 "Agriculture ",modify
label define discoNames 71 "Building Trades",modify
label define discoNames 72 "Metal Trades",modify
label define discoNames 74 "Other Crafts",modify
label define discoNames 81 "Plant Operator",modify
label define discoNames 82 "Machine Operator",modify
label define discoNames 83 "Drivers",modify
label define discoNames 91 "Elementary Occupations",modify
label define discoNames 93 "Laborers",modify
label values disco2 discoNames 

decode disco2, gen(discoStr)
replace discoStr = "" if row==2|row==3

gen missing = v1==""
replace missing = missing+1 if v2=="" & missing>0
replace missing = missing+1 if v3=="" & missing>1
replace missing = missing+1 if v4=="" & missing>2
sort missing disco2 row


gen amp = "&"
gen lb = "\\"

outsheet discoStr amp v1 amp v2 amp v3 amp v4 lb using $tablePath/etaParameters.tex, noq replace non






*******************************
** Non-Employment Parameters **
*******************************
import delimited $statsDKoutput/XQ_BOOTS_TIME.txt, clear
gen param = _n
rename v* val#, addnumber(0)
reshape long val, i(param) j(boot)

* Adjusts for variance
egen rho = max(val*(param==1)),by(boot)
replace val = val/rho
drop rho

reshape wide val, i(boot) j(param)

keep boot val32 - val38

rename val# val#, addnumber(1)
replace val2 = val2*1000


reshape long val, i(boot) j(param)

egen stdErrTime = sd(val),by(param)
keep if boot==0
drop boot val
save temp, replace


import delimited $statsDKoutput/XQ_BOOTS.txt, clear
gen param = _n
rename v* val#, addnumber(0)
reshape long val, i(param) j(boot)

* Adjusts for variance
egen rho = max(val*(param==1)),by(boot)
replace val = val/rho
drop rho

reshape wide val, i(boot) j(param)

keep boot val32 - val38

rename val# val#, addnumber(1)
replace val2 = val2*1000

reshape long val, i(boot) j(param)


egen stdErr = sd(val),by(param)
keep if boot==0
drop boot
merge 1:1 param using temp
drop _merge


format val stdErr* %6.3f
tostring val, replace u force
*replace val = "$" + val + "$"
tostring stdErr, replace u force
replace stdErr = "("+stdErr+")"
tostring stdErrTime, replace u force
replace stdErrTime = "["+stdErrTime+"]"


rename (val stdErr stdErrTime) (v1 v2 v3)
reshape long v, i(param) j(row)


cap label drop nonempParam
label define nonempParam 1 "Age", add 
label define nonempParam 2 "$\text{Age}^2(\times 1000)$", modify 
label define nonempParam 3 "Type 1", modify
label define nonempParam 4 "Type 2", modify
label define nonempParam 5 "Type 3", modify
label define nonempParam 6 "Type 4", modify
label define nonempParam 7 "Type 5", modify
label values param nonempParam 


decode param, gen(paramStr)
replace paramStr = "" if row==2|row==3

gen amp = "&"
gen lb = "\\"

gen famp = "\hspace{10mm} &"

gen bamp = "& \hspace{10mm}"

outsheet famp paramStr amp v bamp lb using $tablePath/nonEmploymentParameters.tex, noq replace non
