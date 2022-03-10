/*
Make Revision Tables

This document takes in output files produced on the Statistics Denmark Servers 
and pretty prints them for input into Latex.

*/



global doFiles "/Users/sharont/Dropbox/DanishLEEDProject/structuralEstimation_Final/Estimation/doFiles"
global dataPath "/Users/sharont/Dropbox/DanishLEEDProject/structuralEstimation_Final/Estimation/inputs"
global outPath "/Users/sharont/Dropbox/DanishLEEDProject/structuralEstimation_Final/Estimation/outputs"
global matlabPath "/Users/sharont/Dropbox/DanishLEEDProject/structuralEstimation_Final/Estimation/matlabFiles"

global figurePath "/Users/sharont/Dropbox/DanishLEEDProject/structuralEstimation_Final/Estimation/figures"
global tablePath "/Users/sharont/Dropbox/DanishLEEDProject/structuralEstimation_Final/Estimation/tables"
global workerFiles /Users/sharont/Dropbox/DanishLEEDProject/structuralEstimation_Final/Calibration/workerFiles

global rev1inputs /Users/SharonT/Dropbox/DanishLEEDProject/structuralEstimation_Final/Estimation/statdk_rev1
global rev2inputs /Users/SharonT/Dropbox/DanishLEEDProject/structuralEstimation_Final/Estimation/statdk_rev2

cd $doFiles

 

* Set 1: Occupation and Sector Capital




**********************
** Parameter Tables **
**********************
* Wage Parameters
import delimited $rev1inputs/BETAMAT_rev1.txt, clear
rename sigma beta10

* Adjusting Age Parameters (just for reading):
* Coef on a to Coef on (a-23)
* Multiply coef on a^2*1000 for readability
*replace beta3 = beta2+23*2*beta3
*replace beta4 = beta3*1000

forval j = 1(1)10{
	egen stdErr`j' = sd(beta`j'),by(newcode)
}

keep if sample==0
drop sample

format beta* stdErr* %6.3f
tostring beta*, replace u force
tostring stdErr*, replace u force

forval j = 1(1)10{
	replace beta`j' = "" if beta`j'=="0.000"
	replace stdErr`j' = "" if beta`j'==""
	
	*replace beta`j' = "$" + beta`j' + "$"
	replace stdErr`j' = "("+stdErr`j'+")" if stdErr`j'!=""
}


rename beta# v#1
rename stdErr# v#2

reshape long v1 v2 v3 v4 v5 v6 v7 v8 v9 v10, i(newcode) j(row)

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

outsheet discoStr amp beta1 amp beta2 amp beta3 amp beta4 amp beta5 amp beta6 amp beta7 amp beta8 amp beta9 amp beta10 lb using $tablePath/wageParameters_MAN_rev1.tex if broadSector==1, noq replace non 
outsheet discoStr amp beta1 amp beta2 amp beta3 amp beta4 amp beta5 amp beta6 amp beta7 amp beta8 amp beta9 amp beta10 lb using $tablePath/wageParameters_SERV_rev1.tex if broadSector==2, noq replace non 
outsheet discoStr amp beta1 amp beta2 amp beta3 amp beta4 amp beta5 amp beta6 amp beta7 amp beta8 amp beta9 amp beta10 lb using $tablePath/wageParameters_FIRE_rev1.tex if broadSector==3, noq replace non 
outsheet discoStr amp beta1 amp beta2 amp beta3 amp beta4 amp beta5 amp beta6 amp beta7 amp beta8 amp beta9 amp beta10 lb using $tablePath/wageParameters_HE_rev1.tex if broadSector==4, noq replace non 



import delimited $rev1inputs/BETAYEARS_rev1.txt, clear
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

outsheet discoStr amp v1 amp v2 amp v3 amp v4 lb using $tablePath/skillPriceParameters_rev1.tex, noq replace non



import delimited $rev1inputs/xQ_BOOTS_rev1.txt, clear
save $rev1inputs/xQ_BOOTS_rev1, replace

import delimited $rev1inputs/xQ_BOOTS_rev1_nocons.txt, clear
save $rev1inputs/xQ_BOOTS_rev1_nocons, replace

*  Inverse Rho
use $rev1inputs/xQ_BOOTS_rev1, clear
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


use $rev1inputs/xQ_BOOTS_rev1_nocons, clear
keep if _n==1
gen param = 1
reshape long v, i(param) j(boot)
replace boot = boot-1

egen stdErr = sd(v)
keep if boot==0
drop boot
drop param
gen spec = 2

append using temp
save temp, replace


sort spec

gen id = 1

format v stdErr %6.3f
tostring v, replace u force
tostring stdErr, replace u force
replace stdErr = "("+stdErr+")"
rename v bv
rename stdErr bstdErr
reshape wide bv bstdErr, i(id) j(spec)


rename b*1 b1*
rename b*2 b2*
reshape long b1 b2, i(id) j(stat) string

gsort -stat
replace id = _n

replace stat = "" if stat!="v"
replace stat = "\$1/\rho$" if stat == "v"
gen amp = "&"
gen lb = "\\"

outsheet stat amp b1 amp b2 lb using $tablePath/rho_rev1.tex, noq replace non 



***** Brings in Average C/rho *****
clear
set obs 0
gen cost = .
save meanGamma_rev1, replace

clear
set obs 0
gen cost = .
save meanGamma_rev1_nocons, replace

use $rev1inputs/xQ_BOOTS_rev1, clear

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
	
	append using meanGamma_rev1
	save meanGamma_rev1, replace
	restore
}


use $rev1inputs/xQ_BOOTS_rev1_nocons, clear

gen param = _n
rename v* val#, addnumber(0)
reshape long val, i(param) j(boot)
reshape wide val, i(boot) j(param)

keep boot val9 - val30

rename val# gammaVal#, addnumber(1)

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
	
	reshape long charVal, i(newCodeA newCodeB) j(char)
	merge m:1 char using gammaTemp
	drop _merge
	gen cost = charVal*gammaVal
	collapse (sum) cost, by(newCodeA newCodeB boot)
	replace cost = exp(cost)
	
	collapse (mean) cost, by(boot)
	
	append using meanGamma_rev1_nocons
	save meanGamma_rev1_nocons, replace
	restore
}



use meanGamma_rev1, clear
rename cost cost1
merge 1:1 boot using meanGamma_rev1_nocons
drop _merge
rename cost cost2

forval j = 1(1)2{
	egen stdErr`j' = sd(cost`j')
}
keep if boot == 0

reshape long cost stdErr, i(boot) j(spec)
gen id = 2
drop boot

format cost stdErr %6.3f
tostring cost, replace u force
tostring stdErr, replace u force
replace stdErr = "("+stdErr+")"

rename cost bcost
rename stdErr bstdErr
reshape wide bcost bstdErr, i(id) j(spec)

rename b*1 b1*
rename b*2 b2*
reshape long b1 b2, i(id) j(stat) string


replace stat = "" if stat!="cost"
replace stat = "Mean \$C/\rho$" if stat == "cost"
gen amp = "&"
gen lb = "\\"

outsheet stat amp b1 amp b2 lb using $tablePath/meanCosts_rev1.tex, noq replace non 






*********************
** Cost Parameters **
*********************

use $rev1inputs/xQ_BOOTS_rev1, clear

gen param = _n
rename v* val#, addnumber(0)
reshape long val, i(param) j(boot)
reshape wide val, i(boot) j(param)

keep boot val9 - val31

rename val# val#, addnumber(0)
reshape long val, i(boot) j(Gamma)

egen stdErr = sd(val),by(Gamma)
keep if boot==0


format val stdErr %6.3f
tostring val, replace u force
*replace val = "$" + val + "$"
tostring stdErr, replace u force
replace stdErr = "("+stdErr+")"

rename (val stdErr) (v1 v2)
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

rename (v0 v1) (v0_cons v1_cons)

save temp, replace


use $rev1inputs/xQ_BOOTS_rev1_nocons, clear

gen param = _n
rename v* val#, addnumber(0)
reshape long val, i(param) j(boot)
reshape wide val, i(boot) j(param)

keep boot val9 - val30

rename val# val#, addnumber(1)
reshape long val, i(boot) j(Gamma)

egen stdErr = sd(val),by(Gamma)
keep if boot==0


format val stdErr %6.3f
tostring val, replace u force
*replace val = "$" + val + "$"
tostring stdErr, replace u force
replace stdErr = "("+stdErr+")"

rename (val stdErr) (v1 v2)
reshape long v, i(boot Gamma) j(row)
drop boot
gen down = 1 if Gamma>=11 & Gamma<=20
replace down = 0 if down==.
replace Gamma = Gamma-10 if down==1
reshape wide v, i(Gamma row) j(down)

replace Gamma = -2 if Gamma==0
replace Gamma = -1 if Gamma==21
replace Gamma = 0 if Gamma==22


merge 1:1 Gamma row using temp
drop _merge

sort Gamma row
replace Gamma = . if row==2
tostring Gamma, replace

gen gammaString = "Task " + Gamma
replace gammaString = "Constant" if Gamma == "-2"
replace gammaString = "Occ. Dummy" if Gamma=="-1"
replace gammaString = "Sec. Dummy" if Gamma=="0"
replace gammaString = "" if row==2

drop Gamma row

gen amp = "&"
gen lb = "\\"

outsheet gammaStr amp v0_cons amp v1_cons amp v0 amp v1 lb using $tablePath/costParameters_rev1.tex, noq replace non 






**************************************
** Mobility Productivity Parameters **
**************************************

use $rev1inputs/xQ_BOOTS_rev1_nocons, clear
gen param = _n
rename v* val#, addnumber(0)
reshape long val, i(param) j(boot)
reshape wide val, i(boot) j(param)

keep boot val2 - val8

* rearranging to put age out front
rename val# val#, addnumber(3)
rename val8 val1
rename val9 val2
reshape long val, i(boot) j(Gamma)

egen stdErr = sd(val),by(Gamma)
keep if boot==0
drop boot

* adjusting Age^2 up to make it readable
replace val = val*1000 if Gamma==2
replace stdErr = stdErr*1000 if Gamma==2

format val stdErr %6.3f
tostring val, replace u force
*replace val = "$" + val + "$"
tostring stdErr, replace u force
replace stdErr = "("+stdErr+")"

rename (val stdErr) (v1 v2)
reshape long v, i(Gamma) j(row)
rename v v_nocons
save temp, replace


use $rev1inputs/xQ_BOOTS_rev1, clear
gen param = _n
rename v* val#, addnumber(0)
reshape long val, i(param) j(boot)
reshape wide val, i(boot) j(param)

keep boot val2 - val8

* rearranging to put age out front
rename val# val#, addnumber(3)
rename val8 val1
rename val9 val2
reshape long val, i(boot) j(Gamma)

egen stdErr = sd(val),by(Gamma)
keep if boot==0
drop boot

* adjusting Age^2 up to make it readable
replace val = val*1000 if Gamma==2
replace stdErr = stdErr*1000 if Gamma==2

format val stdErr %6.3f
tostring val, replace u force
*replace val = "$" + val + "$"
tostring stdErr, replace u force
replace stdErr = "("+stdErr+")"

rename (val stdErr) (v1 v2)
reshape long v, i(Gamma) j(row)
merge 1:1 Gamma row using temp
drop _merge

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
replace gammaStr = "" if row==2

gen amp = "&"
gen lb = "\\"

gen famp = "\hspace{6mm} &"

gen bamp = "& \hspace{6mm}"

outsheet famp gammaStr amp v amp v_nocons bamp lb using $tablePath/costProductivity_rev1.tex, noq replace non





*******************************
** Non-Employment Parameters **
*******************************
use $rev1inputs/xQ_BOOTS_rev1_nocons, clear
gen param = _n
rename v* val#, addnumber(0)
reshape long val, i(param) j(boot)

* Adjusts for variance
egen rho = max(val*(param==1)),by(boot)
replace val = val/rho
drop rho

reshape wide val, i(boot) j(param)

keep boot val31 - val37

rename val# val#, addnumber(1)
reshape long val, i(boot) j(param)


* For presentation purposes -- shift age range from [0,35] back to [25,60]
* This will make the numbers commensurate with wage parameters.
* Note: Linear rescaling is irrelevant for actual counterfactuals/presentation
gen temp = val if param==1
egen bAge = max(temp), by(boot)
drop temp

gen temp = val if param==2
egen bAge2 = max(temp), by(boot)
drop temp

replace val = bAge - 50*bAge2 if param==1
replace val = val - 25*bAge + 625*bAge2 if param>3
drop bAge bAge2


egen stdErr = sd(val),by(param)
keep if boot==0
drop boot


format val stdErr %6.3f
tostring val, replace u force
*replace val = "$" + val + "$"
tostring stdErr, replace u force
replace stdErr = "("+stdErr+")"


rename (val stdErr) (v1 v2)
reshape long v, i(param) j(row)
rename v v_nocons
save temp, replace


use $rev1inputs/xQ_BOOTS_rev1, clear
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
reshape long val, i(boot) j(param)

* For presentation purposes -- shift age range from [0,35] back to [25,60]
* This will make the numbers commensurate with wage parameters.
* Note: Linear rescaling is irrelevant for actual counterfactuals/presentation
gen temp = val if param==1
egen bAge = max(temp), by(boot)
drop temp

gen temp = val if param==2
egen bAge2 = max(temp), by(boot)
drop temp

replace val = bAge - 50*bAge2 if param==1
replace val = val - 25*bAge + 625*bAge2 if param>3
drop bAge bAge2


egen stdErr = sd(val),by(param)
keep if boot==0
drop boot


format val stdErr %6.3f
tostring val, replace u force
*replace val = "$" + val + "$"
tostring stdErr, replace u force
replace stdErr = "("+stdErr+")"


rename (val stdErr) (v1 v2)
reshape long v, i(param) j(row)
merge 1:1 param row using temp
drop _merge

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
replace paramStr = "" if row==2

gen amp = "&"
gen lb = "\\"

gen famp = "\hspace{6mm} &"

gen bamp = "& \hspace{6mm}"

outsheet famp paramStr amp v amp v_nocons bamp lb using $tablePath/nonEmployment_rev1.tex, noq replace non








********************
** Eta Parameters **
********************
use $rev1inputs/xQ_BOOTS_rev1_nocons, clear
gen param = _n
rename v* val#, addnumber(0)
reshape long val, i(param) j(boot)

* Adjusts for variance
egen rho = max(val*(param==1)),by(boot)
replace val = val/rho
drop rho

reshape wide val, i(boot) j(param)

keep boot val38 - val75
rename val# eta#, addnumber(1)
reshape long eta, i(boot) j(newcode)

egen stdErr = sd(eta),by(newcode)
keep if boot==0
drop boot

format eta stdErr %6.3f
tostring eta, replace u force
*replace eta = "$" + eta + "$"
tostring stdErr, replace u force
replace stdErr = "("+stdErr+")"


rename (eta stdErr) (v1 v2)
reshape long v, i(newcode) j(row)
rename v v_nocons
save temp, replace

use $rev1inputs/xQ_BOOTS_rev1, clear
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

save $tablePath/etaVecRev1, replace

format eta stdErr %6.3f
tostring eta, replace u force
*replace eta = "$" + eta + "$"
tostring stdErr, replace u force
replace stdErr = "("+stdErr+")"


rename (eta stdErr) (v1 v2)
reshape long v, i(newcode) j(row)
merge 1:1 newcode row using temp
drop _merge

merge m:1 newcode using $workerFiles/newcodenames
drop _merge

drop newcode
reshape wide v v_nocons, i(disco2 row) j(broadSector)

cap label drop discoNames
label define discoNames 12 "Managers",add
label define discoNames 21 "Science Pro",modify
label define discoNames 22 "Health Pro",modify
label define discoNames 23 "Teachers",modify
label define discoNames 24 "Other Pro",modify
label define discoNames 31 "Science Assc. Pro",modify
label define discoNames 32 "Health Assc. Pro",modify
label define discoNames 33 "Teaching Assc. Pro",modify
label define discoNames 34 "Other Assc. Pro",modify
label define discoNames 41 "Clerks",modify
label define discoNames 42 "Customer Service",modify
label define discoNames 51 "Personal Workers",modify
label define discoNames 52 "Retail Workers",modify
label define discoNames 61 "Agriculture ",modify
label define discoNames 71 "Building Trades",modify
label define discoNames 72 "Metal Trades",modify
label define discoNames 74 "Other Crafts",modify
label define discoNames 81 "Plant Op",modify
label define discoNames 82 "Machine Op",modify
label define discoNames 83 "Drivers",modify
label define discoNames 91 "Elementary Occs",modify
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

outsheet discoStr amp v1 amp v2 amp v3 amp v4 amp v_nocons1 amp v_nocons2 amp v_nocons3 amp v_nocons4 lb using $tablePath/etaParameters_rev1.tex, noq replace non







* Set 2: Alternative Estimator

import delimited $rev2inputs/xQ_BOOTS_rev2.txt, clear
save $rev2inputs/xQ_BOOTS_rev2, replace

* First, rho
use $rev2inputs/xQ_BOOTS_rev2, clear
keep if _n==1
gen param = 1
reshape long v, i(param) j(boot)
replace boot = boot-1

egen stdErr = sd(v)
keep if boot==0
drop boot
drop param

gen id = 1

format v stdErr %6.3f
tostring v, replace u force
tostring stdErr, replace u force
replace stdErr = "("+stdErr+")"
rename v bv
rename stdErr bstdErr
reshape long b, i(id) j(stat) string

gsort -stat
replace id = _n

replace stat = "" if stat!="v"
replace stat = "\$1/\rho$" if stat == "v"
gen amp = "&"
gen lb = "\\"

outsheet stat amp b lb using $tablePath/rho_rev2.tex, noq replace non 







***** Brings in Average C/rho *****
clear
set obs 0
gen cost = .
save meanGamma_rev2, replace

use $rev2inputs/xQ_BOOTS_rev2, clear
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
	
	append using meanGamma_rev2
	save meanGamma_rev2, replace
	restore
}





use meanGamma_rev2, clear
egen stdErr = sd(cost)
keep if boot == 0

gen id = 2
drop boot

format cost stdErr %6.3f
tostring cost, replace u force
tostring stdErr, replace u force
replace stdErr = "("+stdErr+")"

rename cost bcost
rename stdErr bstdErr

reshape long b, i(id) j(stat) string


replace stat = "" if stat!="cost"
replace stat = "Mean \$C/\rho$" if stat == "cost"
gen amp = "&"
gen lb = "\\"

outsheet stat amp b lb using $tablePath/meanCosts_rev2.tex, noq replace non 







*********************
** Cost Parameters **
*********************
use $rev2inputs/xQ_BOOTS_rev2, clear
gen param = _n
rename v* val#, addnumber(0)
reshape long val, i(param) j(boot)
reshape wide val, i(boot) j(param)

keep boot val9 - val31
rename val# val#, addnumber(0)
reshape long val, i(boot) j(Gamma)

egen stdErr = sd(val),by(Gamma)
keep if boot==0

format val stdErr %6.3f
tostring val, replace u force
*replace val = "$" + val + "$"
tostring stdErr, replace u force
replace stdErr = "("+stdErr+")"

rename (val stdErr) (v1 v2)
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
replace Gamma = . if row==2
tostring Gamma, replace

gen gammaString = "Task " + Gamma
replace gammaString = "Constant" if Gamma == "-2"
replace gammaString = "Occ. Dummy" if Gamma=="-1"
replace gammaString = "Sec. Dummy" if Gamma=="0"
replace gammaString = "" if row==2

drop Gamma row

gen amp = "&"
gen lb = "\\"

outsheet gammaStr amp v0 amp v1 lb using $tablePath/costParameters_rev2.tex, noq replace non 



**************************************
** Mobility Productivity Parameters **
**************************************
use $rev2inputs/xQ_BOOTS_rev2, clear
gen param = _n
rename v* val#, addnumber(0)
reshape long val, i(param) j(boot)
reshape wide val, i(boot) j(param)

keep boot val2 - val8

* rearranging to put age out front
rename val# val#, addnumber(3)
rename val8 val1
rename val9 val2
reshape long val, i(boot) j(Gamma)

egen stdErr = sd(val),by(Gamma)
keep if boot==0
drop boot

* adjusting Age^2 up to make it readable
replace val = val*1000 if Gamma==2
replace stdErr = stdErr*1000 if Gamma==2

format val stdErr %6.3f
tostring val, replace u force
*replace val = "$" + val + "$"
tostring stdErr, replace u force
replace stdErr = "("+stdErr+")"

rename (val stdErr) (v1 v2)
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
replace gammaStr = "" if row==2

gen amp = "&"
gen lb = "\\"

gen famp = "\hspace{6mm} &"

gen bamp = "& \hspace{6mm}"

outsheet famp gammaStr amp v bamp lb using $tablePath/costProductivity_rev2.tex, noq replace non






*******************************
** Non-Employment Parameters **
*******************************
use $rev2inputs/xQ_BOOTS_rev2, clear
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
reshape long val, i(boot) j(param)

* For presentation purposes -- shift age range from [0,35] back to [25,60]
* This will make the numbers commensurate with wage parameters.
* Note: Linear rescaling is irrelevant for actual counterfactuals/presentation
gen temp = val if param==1
egen bAge = max(temp), by(boot)
drop temp

gen temp = val if param==2
egen bAge2 = max(temp), by(boot)
drop temp

replace val = bAge - 50*bAge2 if param==1
replace val = val - 25*bAge + 625*bAge2 if param>3
drop bAge bAge2


egen stdErr = sd(val),by(param)
keep if boot==0
drop boot


format val stdErr %6.3f
tostring val, replace u force
*replace val = "$" + val + "$"
tostring stdErr, replace u force
replace stdErr = "("+stdErr+")"


rename (val stdErr) (v1 v2)
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
replace paramStr = "" if row==2

gen amp = "&"
gen lb = "\\"

gen famp = "\hspace{10mm} &"

gen bamp = "& \hspace{10mm}"

outsheet famp paramStr amp v bamp lb using $tablePath/nonEmploymentParameters_rev2.tex, noq replace non






********************
** Eta Parameters **
********************
use $rev2inputs/xQ_BOOTS_rev2, clear
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

save $tablePath/etaVec_rev2, replace

format eta stdErr %6.3f
tostring eta, replace u force
*replace eta = "$" + eta + "$"
tostring stdErr, replace u force
replace stdErr = "("+stdErr+")"


rename (eta stdErr) (v1 v2)
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

outsheet discoStr amp v1 amp v2 amp v3 amp v4 lb using $tablePath/etaParameters_rev2.tex, noq replace non




**************************************
** Mobility Productivity Parameters **
**************************************
use $rev2inputs/xQ_BOOTS_rev2, clear
gen param = _n
rename v* val#, addnumber(0)
reshape long val, i(param) j(boot)
reshape wide val, i(boot) j(param)

keep boot val2 - val8

* rearranging to put age out front
rename val# val#, addnumber(3)
rename val8 val1
rename val9 val2
reshape long val, i(boot) j(Gamma)

egen stdErr = sd(val),by(Gamma)
keep if boot==0
drop boot

* adjusting Age^2 up to make it readable
replace val = val*1000 if Gamma==2
replace stdErr = stdErr*1000 if Gamma==2

format val stdErr %6.3f
tostring val, replace u force
*replace val = "$" + val + "$"
tostring stdErr, replace u force
replace stdErr = "("+stdErr+")"

rename (val stdErr) (v1 v2)
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
replace gammaStr = "" if row==2

gen amp = "&"
gen lb = "\\"

gen famp = "\hspace{10mm} &"

gen bamp = "& \hspace{10mm}"

outsheet famp gammaStr amp v bamp lb using $tablePath/mobility_rev2.tex, noq replace non










use $tablePath/etaVec_rev2, clear
drop stdErr
rename eta eta_rev2
merge 1:1 newcode using $tablePath/etaVec
drop _merge
drop stdErr*


regress eta_rev2 eta
local b0 = int(_b[_cons]*100)/100
local b1 = int(_b[eta]*100)/100
local r2 = int(e(r2)*100)/100

graph tw (scatter eta_rev2 eta) (lfit eta_rev2 eta), text(-2.05 1.9 "y = `b0' + `b1'x,  R{superscript:2}=`r2'",orient(horizontal) placement(e) si(vsmall) justification(center) box fcolor(white) m(vsmall) ) scheme(s1mono) legend(off) xtitle("Original Estimate of {&eta}")  ytitle("New Estimate of {&eta}")
gr export $figurePath/etaComparison.pdf, replace


