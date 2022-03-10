global doFiles /Users/sharontraiberman/Dropbox/DanishLEEDProject/renewalEstimationIncome/Calibration/doFiles
global foreignPriceFiles /Users/sharontraiberman/Dropbox/DanishLEEDProject/renewalEstimationIncome/Calibration/foreignPriceIndices
global csvFiles /Users/sharontraiberman/Dropbox/DanishLEEDProject/renewalEstimationIncome/Calibration/csvFiles
global outFiles /Users/sharontraiberman/Dropbox/DanishLEEDProject/renewalEstimationIncome/Calibration/outFiles
global workerFiles /Users/sharontraiberman/Dropbox/DanishLEEDProject/renewalEstimationIncome/Calibration/workerFiles
global prodcomFiles /Users/sharontraiberman/Dropbox/DanishLEEDProject/renewalEstimationIncome/Calibration/prodcomFiles
global figurePath /Users/sharontraiberman/Dropbox/DanishLEEDProject/figures
global witsData /Users/sharontraiberman/Dropbox/DanishLEEDProject/renewalEstimationIncome/Calibration/witsFiles
global externalData /Users/sharontraiberman/Dropbox/DanishLEEDProject/renewalEstimationIncome/Calibration/externalData
global unctadData /Users/sharontraiberman/Dropbox/DanishLEEDProject/renewalEstimationIncome/dekleEatonKortum/unctad_data_new

cd $unctadData

use ../dekCountries, clear
drop if country == "china, hong kong sar"
drop if country == "china, macao sar"
drop if country == "indonesia (...2002)"
save ../dekCountries2ISO, replace
 
use $prodcomFiles/nace2hs, clear
gen year = 1995
forval j = 1996/2008{
append using $prodcomFiles/nace2hs
replace year = `j' if year==.
}
save nace2hs_year, replace

use $prodcomFiles/nace2hs, clear
gen hs4 = int(hs6/100)
gen hs2 = int(hs4/100)
merge m:1 hs6 using $witsData/elasticities/hs6Sigmas
drop if _merge==2
drop _merge
merge m:1 hs4 using $witsData/elasticities/hs4Sigmas, update
drop if _merge==2
drop _merge
merge m:1 hs2 using $witsData/elasticities/hs2Sigmas, update
drop if _merge==2
drop _merge
drop if sigma==.
collapse (min) minSigma = sigma (median) medSigma = sigma, by(nace2)

save nace2sigmas, replace

use $witsData/outFiles/denmarkImportFlows, clear

* Convering to DKK...
merge m:1 year using $externalData/usddk_forex,
keep if _merge==3
replace tradevalue = tradevalue * usd_dkk

* Dropping World and Other
drop if partneriso3=="WLD" | partneriso3=="UNS" | partneriso3=="SPE"

replace partneriso3="BEL" if partneriso3=="BLX" | partneriso3=="LUX"
replace partneriso3="ROU" if partneriso3=="ROM"
replace partneriso3="SRB" if partneriso3=="SER"

collapse (sum)  quantity tradevalue,by(year reporteriso3 partneriso3 hs6 qtytoken)

/*
rename partneriso3 alpha3
merge m:1 alpha3 using iso2region
drop if _merge==2
drop _merge

replace subregioncode = 30 if alpha3=="OAS"
drop if subregioncode == .

collapse (sum)  quantity tradevalue,by(year reporteriso3 subregioncode productcode qtytoken)
rename subregioncode cid
*/

encode partneriso3, gen(cid)

gen hs4 = int(hs6/100)
gen hs2 = int(hs4/100)

gen unitvalue = tradevalue/quantity
drop if unitvalue==. | unitvalue==0

* Deflates by the CPI -- as wages and domestic prices
gen cpi = 82.9 if year==1991
replace cpi = 80.9 if year ==1990
replace cpi = 84.6 if year ==1992
replace cpi = 85.7 if year ==1993
replace cpi = 87.4 if year ==1994
replace cpi = 89.2 if year ==1995
replace cpi = 91.1 if year ==1996
replace cpi = 93.1 if year ==1997
replace cpi = 94.8 if year ==1998
replace cpi = 97.2 if year ==1999
replace cpi = 100 if year ==2000
replace cpi = 102.4 if year ==2001
replace cpi = 104.8 if year ==2002
replace cpi = 107 if year ==2003
replace cpi = 108 if year ==2004
replace cpi = 110.2 if year ==2005
replace cpi = 112.3 if year ==2006
replace cpi = 114.2 if year ==2007
replace cpi = 118.1 if year ==2008
replace cpi = 119.7 if year ==2009
replace cpi = 122.4 if year ==2010
replace cpi = 125.8 if year ==2011

replace unitvalue = unitvalue*100/cpi

/*
Getting the right value of sigma for each good
*/

merge m:1 hs6 using $witsData/elasticities/hs6Sigmas
drop if _merge==2
drop _merge


merge m:1 hs4 using $witsData/elasticities/hs4Sigmas, update
drop if _merge==2
drop _merge

merge m:1 hs2 using $witsData/elasticities/hs2Sigmas, update
drop if _merge==2
drop _merge

egen varID = group(hs6 cid)
xtset varID year


save temp, replace

use temp, clear

collapse (sum) tradevalue, by(hs6 year)
rename tradevalue goodTotal

save goodTotals, replace


use temp, clear

rename partneriso3 iso3
merge m:1 iso3 using ../dekCountries2ISO
drop if _merge==2
replace country = "row" if _merge==1

replace country = "singapore" if country=="malaysia"|country=="philippines"

drop _merge

collapse (sum) tradevalue, by(hs6 sigma year country)

egen total = total(tradevalue), by(hs6 year)
gen share = tradevalue/total


merge m:1 year country using dkTradeCosts
keep if _merge==3
drop _merge

* DECIDE WHOM TO KEEP HERE!!!!!!
*keep if country=="china"
*keep if country=="china"|country=="poland"|country=="lithuania"|country=="turkey"


drop X*
egen tau96 = max(tau*(year==1996)), by(country)

gen tauHat = tau96/tau

gen pgHat = share*tauHat^(1-sigma)

collapse (sum) share pgHat, by(hs6 year sigma)
replace pgHat = (1-share + pgHat)^(1/(1-sigma))

*egen inShare = total(share), by(hs6 year)
*gen pgHat = (1-share + share*tauHat^(1-sigma))^(1/(1-sigma))

keep hs6 year pgHat

save $unctadData/priceChanges, replace


/*
use $prodcomFiles/nace2hs, clear

merge m:m hs6 using goodTotals
keep if _merge==3
drop _merge
*/

use $unctadData/goodTotals, clear
merge 1:m hs6 year using $unctadData/nace2hs_year
keep if _merge==3
drop _merge

merge m:1 hs6 year using $unctadData/priceChanges
drop if _merge==2
replace pgHat = 1 if pgHat==.

keep nace2 hs6 year pgHat goodTotal

egen total = total(goodTotal), by(nace2 year)

gen share = goodTotal/total

merge m:1 nace2 using $unctadData/nace2sigmas
keep if _merge==3
drop _merge

gen pFHat = pgHat^(1-minSigma)*share

collapse (sum) pFHat, by(nace2 year minSigma)
replace pFHat = pFHat^(1/(1-minSigma)) 

xtset nace2 year
tsfill, full
drop minSigma 

replace pFHat = 1 if year==1996
replace pFHat = L.pFHat if pFHat==.

save $unctadData/counterfactualChanges, replace


use  $outFiles/smoothedForeignPriceChanges, clear
merge 1:1 nace2 year using $unctadData/counterfactualChanges
drop if _merge==2
drop _merge 

replace pFHat = 1 if pFHat==.

xtset nace2 year
gen newDex = smoothedDex* pFHat/L.pFHat
replace newDex = smoothedDex if newDex==.

merge m:1 nace2 using $matlabFiles/nace2ind
drop _merge nace2

keep if tradable==1
drop tradable

replace year = year-1996

sort year indCode
order year indCode

outsheet year indCode newDex using $matlabFiles/DELPF_cf.csv, c replace nol non noq
outsheet year indCode pFHat using $matlabFiles/PRICERATIO_cf.csv, c replace nol non noq
