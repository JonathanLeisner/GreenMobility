/*********************
Counterfactuals Analysis
*********************/

global postFiles /Users/sharontraiberman/Dropbox/DanishLEEDProject/renewalEstimationIncome/Counterfactuals/ca_pf
global sandbox /Users/sharontraiberman/Dropbox/DanishLEEDProject/renewalEstimationIncome/Counterfactuals/ca_pf

cd $sandbox

/***************
  REAL WAGES
****************/


insheet using $sandbox/meanRealWageDiffs.csv, c clear
gen newcode = _n
order newcode

rename v* realWageDiff*
reshape long realWageDiff, i(newcode) j(period)
save $postFiles/wageDiffs, replace

merge m:1 newcode using $workerFiles/newcodenames
drop _merge

cap label drop codeNames
label define codeNames 1 "Managers",add
label define codeNames 2 "Science Professional",modify
label define codeNames 3 "Science Associate Professional",modify
label define codeNames 4 "Other Assc. Professional",modify
label define codeNames 5 "Clerks",modify
label define codeNames 6 "Agriculture ",modify
label define codeNames 7 "Building Trades",modify
label define codeNames 8 "Metal Trades",modify
label define codeNames 9 "Other Crafts",modify
label define codeNames 10 "Plant Operator",modify
label define codeNames 11 "Machine Operator",modify
label define codeNames 12 "Drivers",modify
label define codeNames 13 "Laborers",modify
label define codeNames 14 "Managers",modify
label define codeNames 15 "Science Professional",modify
label define codeNames 16 "Other Professional",modify
label define codeNames 17 "Science Associate Professional",modify
label define codeNames 18 "Other Assc. Professional",modify
label define codeNames 19 "Clerks",modify
label define codeNames 20 "Personal Workers",modify
label define codeNames 21 "Retail Workers",modify
label define codeNames 22 "Metal Trades",modify
label define codeNames 23 "Drivers",modify
label define codeNames 24 "Elementary Occupations",modify
label define codeNames 25 "Laborers",modify
label define codeNames 26 "Managers",modify
label define codeNames 27 "Science Professional",modify
label define codeNames 28 "Other Professional",modify
label define codeNames 29 "Science Associate Professional",modify
label define codeNames 30 "Other Assc. Professional",modify
label define codeNames 31 "Clerks",modify
label define codeNames 32 "Customer Service",modify
label define codeNames 33 "Health Professional",modify
label define codeNames 34 "Teachers",modify
label define codeNames 35 "Health Assc. Professional",modify
label define codeNames 36 "Teaching Assc. Professional",modify
label define codeNames 37 "Clerks",modify
label define codeNames 38 "Personal Workers",modify

label values newcode codeNames

gen year = period-9
xtset newcode year

su period
local max = r(max)

egen finalDiff = min(realWageDiff*(period==`max')), by(newcode)

gen percentOff =  (realWageDiff -  finalDiff)*100

*tsline percentOff if year>=-2, by(newcode, title("Convergence to New Steady State: Real Wages") note("Time 0 is normalized to first period after shocks (2006)""Difference is weighted by worker composition.", si(tiny))) xtitle("Year") ytitle("Percent Difference in Wage from Final SS")
*graph export $postFiles/convergenceTime.pdf, replace

replace realWageDiff = realWageDiff*100

*tsline realWageDiff if year>=-2, by(newcode, title("Differences Real Wages") note("Time 0 is normalized to first period after shocks (2006)""Difference is weighted by worker composition.", si(tiny))) xtitle("Year") ytitle("Real Wage Differences")
*graph export $postFiles/wageDiffs.pdf, replace

*xtline percentOff, overlay legend(off) ttitle(Period) ytitle(Percent Difference from Steady State)
*graph export $postFiles/convergenceTime_$flag.pdf, replace

*graph tw (kdensity realWageDiff if year==0, lc(gs0)) (kdensity realWageDiff if year==5, lc(gs8)) (kdensity realWageDiff if year==30, lc(gs13)), legend(lab(1 "Shocks End") lab(2 "+5 Periods") lab(3 "+25 Periods") si(small) row(1) symx(3) symy(2)) xtitle("w{subscript:Actual Prices} - w{subscript:Fixed Prices}") ytitle("Freq.") title("Density of Skill Price Differences") note("Differences are log points x 100""Unit of observation is occupation, no weights.",si(vsmall))



cap label drop codeNames
label define codeNames 1 "MANUFACTURING & Managers",add
label define codeNames 2 "&Science Professional",modify
label define codeNames 3 "&Science Assc. Professional",modify
label define codeNames 4 "&Other Assc. Professional",modify
label define codeNames 5 "&Clerks",modify
label define codeNames 6 "&Agriculture ",modify
label define codeNames 7 "&Building Trades",modify
label define codeNames 8 "&Metal Trades",modify
label define codeNames 9 "&Other Crafts",modify
label define codeNames 10 "&Plant Operator",modify
label define codeNames 11 "&Machine Operator",modify
label define codeNames 12 "&Drivers",modify
label define codeNames 13 "&Laborers",modify
label define codeNames 14 "SERVICES &Managers",modify
label define codeNames 15 "&Science Professional",modify
label define codeNames 16 "&Other Professional",modify
label define codeNames 17 "&Science Assc. Professional",modify
label define codeNames 18 "&Other Assc. Professional",modify
label define codeNames 19 "&Clerks",modify
label define codeNames 20 "&Personal Workers",modify
label define codeNames 21 "&Retail Workers",modify
label define codeNames 22 "&Metal Trades",modify
label define codeNames 23 "&Drivers",modify
label define codeNames 24 "&Elementary Occupations",modify
label define codeNames 25 "&Laborers",modify
label define codeNames 26 "FIRE &Managers",modify
label define codeNames 27 "&Science Professional",modify
label define codeNames 28 "&Other Professional",modify
label define codeNames 29 "&Science Assc. Professional",modify
label define codeNames 30 "&Other Assc. Professional",modify
label define codeNames 31 "&Clerks",modify
label define codeNames 32 "&Customer Service",modify
label define codeNames 33 "HEALTH \& EDUC. &Health Professional",modify
label define codeNames 34 "&Teachers",modify
label define codeNames 35 "&Health Assc. Professional",modify
label define codeNames 36 "&Teaching Assc. Professional",modify
label define codeNames 37 "&Clerks",modify
label define codeNames 38 "&Personal Workers",modify
label values newcode codeNames

keep newcode realWageDiff period broadSector
reshape wide realWageDiff, i(newcode) j(period)

format realWageDiff* %4.2f

estpost tabstat realWageDiff9 realWageDiff19 realWageDiff39 if broadSector==1, by(newcode) statistics(mean) nototal
esttab using $postFiles/realWageDiffsTable_man.tex, cells("realWageDiff9(fmt(2) l({})) realWageDiff39(fmt(2) l({})) realWageDiff26(fmt(2) l({}))") varlabels(`e(labels)')  replace  plain fragment noobs nonumber nomtitles

estpost tabstat realWageDiff9 realWageDiff19 realWageDiff39 if broadSector==2, by(newcode) statistics(mean) nototal
esttab using $postFiles/realWageDiffsTable_serv.tex, cells("realWageDiff9(fmt(2) l({})) realWageDiff39(fmt(2) l({})) realWageDiff26(fmt(2) l({}))") varlabels(`e(labels)')  replace  plain fragment noobs nonumber nomtitles

estpost tabstat realWageDiff9 realWageDiff19 realWageDiff39 if broadSector==3, by(newcode) statistics(mean) nototal
esttab using $postFiles/realWageDiffsTable_fire.tex, cells("realWageDiff9(fmt(2) l({})) realWageDiff39(fmt(2) l({})) realWageDiff26(fmt(2) l({}))") varlabels(`e(labels)')  replace  plain fragment noobs nonumber nomtitles

estpost tabstat realWageDiff9 realWageDiff19 realWageDiff39 if broadSector==4, by(newcode) statistics(mean) nototal
esttab using $postFiles/realWageDiffsTable_he.tex, cells("realWageDiff9(fmt(2) l({})) realWageDiff39(fmt(2) l({})) realWageDiff26(fmt(2) l({}))") varlabels(`e(labels)')  replace  plain fragment noobs nonumber nomtitles


save $postFiles/wageComparison, replace








insheet using $sandbox/skillPriceDiffs.csv, c clear
gen newcode = _n
order newcode

rename v* skillPriceDiff*

reshape long skillPriceDiff, i(newcode) j(period)

replace skillPriceDiff = skillPriceDiff*100

reshape wide skillPriceDiff, i(newcode) j(period)

merge m:1 newcode using $workerFiles/newcodenames
drop _merge

cap label drop codeNames
label define codeNames 1 "MANUFACTURING & Managers",add
label define codeNames 2 "&Science Professional",modify
label define codeNames 3 "&Science Assc. Professional",modify
label define codeNames 4 "&Other Assc. Professional",modify
label define codeNames 5 "&Clerks",modify
label define codeNames 6 "&Agriculture ",modify
label define codeNames 7 "&Building Trades",modify
label define codeNames 8 "&Metal Trades",modify
label define codeNames 9 "&Other Crafts",modify
label define codeNames 10 "&Plant Operator",modify
label define codeNames 11 "&Machine Operator",modify
label define codeNames 12 "&Drivers",modify
label define codeNames 13 "&Laborers",modify
label define codeNames 14 "SERVICES &Managers",modify
label define codeNames 15 "&Science Professional",modify
label define codeNames 16 "&Other Professional",modify
label define codeNames 17 "&Science Assc. Professional",modify
label define codeNames 18 "&Other Assc. Professional",modify
label define codeNames 19 "&Clerks",modify
label define codeNames 20 "&Personal Workers",modify
label define codeNames 21 "&Retail Workers",modify
label define codeNames 22 "&Metal Trades",modify
label define codeNames 23 "&Drivers",modify
label define codeNames 24 "&Elementary Occupations",modify
label define codeNames 25 "&Laborers",modify
label define codeNames 26 "FIRE &Managers",modify
label define codeNames 27 "&Science Professional",modify
label define codeNames 28 "&Other Professional",modify
label define codeNames 29 "&Science Assc. Professional",modify
label define codeNames 30 "&Other Assc. Professional",modify
label define codeNames 31 "&Clerks",modify
label define codeNames 32 "&Customer Service",modify
label define codeNames 33 "HEALTH \& EDUC. &Health Professional",modify
label define codeNames 34 "&Teachers",modify
label define codeNames 35 "&Health Assc. Professional",modify
label define codeNames 36 "&Teaching Assc. Professional",modify
label define codeNames 37 "&Clerks",modify
label define codeNames 38 "&Personal Workers",modify
label values newcode codeNames



format skillPriceDiff* %4.2f

estpost tabstat skillPriceDiff9 skillPriceDiff19 skillPriceDiff39 if broadSector==1, by(newcode) statistics(mean) nototal
esttab using $postFiles/skillPriceDiffsTable_man.tex, cells("skillPriceDiff9(fmt(2) l({})) skillPriceDiff19(fmt(2) l({})) skillPriceDiff39(fmt(2) l({}))") varlabels(`e(labels)')  replace  plain fragment noobs nonumber nomtitles

estpost tabstat skillPriceDiff9 skillPriceDiff19 skillPriceDiff39 if broadSector==2, by(newcode) statistics(mean) nototal
esttab using $postFiles/skillPriceDiffsTable_serv.tex, cells("skillPriceDiff9(fmt(2) l({})) skillPriceDiff19(fmt(2) l({})) skillPriceDiff39(fmt(2) l({}))") varlabels(`e(labels)')  replace  plain fragment noobs nonumber nomtitles

estpost tabstat skillPriceDiff9 skillPriceDiff19 skillPriceDiff39 if broadSector==3, by(newcode) statistics(mean) nototal
esttab using $postFiles/skillPriceDiffsTable_fire.tex, cells("skillPriceDiff9(fmt(2) l({})) skillPriceDiff19(fmt(2) l({})) skillPriceDiff39(fmt(2) l({}))") varlabels(`e(labels)')  replace  plain fragment noobs nonumber nomtitles

estpost tabstat skillPriceDiff9 skillPriceDiff19 skillPriceDiff39 if broadSector==4, by(newcode) statistics(mean) nototal
esttab using $postFiles/skillPriceDiffsTable_he.tex, cells("skillPriceDiff9(fmt(2) l({})) skillPriceDiff19(fmt(2) l({})) skillPriceDiff39(fmt(2) l({}))") varlabels(`e(labels)')  replace  plain fragment noobs nonumber nomtitles



save $postFiles/skillPriceDiff, replace




insheet using occDistPathcf_pf.csv, c clear
gen newcode = _n
reshape long v, i(newcode) j(period)
rename v occDistCF
save $postFiles/occDist, replace

insheet using occDistPathT_pf.csv, c clear
gen newcode = _n
reshape long v, i(newcode) j(period)
rename v occDistT
save $postFiles/occDistTrade, replace
merge 1:1 newcode period using $postFiles/occDist
drop _merge


su period
local max = r(max)
keep if period==`max'

rename occDistT ssOccDistT
rename occDistCF ssOccDist

cap label drop codeNames
label define codeNames 1 "Managers",add
label define codeNames 2 "Science Professional",modify
label define codeNames 3 "Science Associate Professional",modify
label define codeNames 4 "Other Assc. Professional",modify
label define codeNames 5 "Clerks",modify
label define codeNames 6 "Agriculture ",modify
label define codeNames 7 "Building Trades",modify
label define codeNames 8 "Metal Trades",modify
label define codeNames 9 "Other Crafts",modify
label define codeNames 10 "Plant Operator",modify
label define codeNames 11 "Machine Operator",modify
label define codeNames 12 "Drivers",modify
label define codeNames 13 "Laborers",modify
label define codeNames 14 "Managers",modify
label define codeNames 15 "Science Professional",modify
label define codeNames 16 "Other Professional",modify
label define codeNames 17 "Science Associate Professional",modify
label define codeNames 18 "Other Assc. Professional",modify
label define codeNames 19 "Clerks",modify
label define codeNames 20 "Personal Workers",modify
label define codeNames 21 "Retail Workers",modify
label define codeNames 22 "Metal Trades",modify
label define codeNames 23 "Drivers",modify
label define codeNames 24 "Elementary Occupations",modify
label define codeNames 25 "Laborers",modify
label define codeNames 26 "Managers",modify
label define codeNames 27 "Science Professional",modify
label define codeNames 28 "Other Professional",modify
label define codeNames 29 "Science Associate Professional",modify
label define codeNames 30 "Other Assc. Professional",modify
label define codeNames 31 "Clerks",modify
label define codeNames 32 "Customer Service",modify
label define codeNames 33 "Health Professional",modify
label define codeNames 34 "Teachers",modify
label define codeNames 35 "Health Assc. Professional",modify
label define codeNames 36 "Teaching Assc. Professional",modify
label define codeNames 37 "Clerks",modify
label define codeNames 38 "Personal Workers",modify
label values newcode codeNames

merge m:1 newcode using $workerFiles/newcodenames
drop _merge

gen delShare = (log(ssOccDistT) - log(ssOccDist))*100
separate delShare, by(broadSector)
rename delShare dShareFull


graph hbar delShare*, over(newcode, sort(dShareFull) desc label(labsize(vsmall))) ytitle("") title("Change in Share of Supply") nofill  legend(si(vsmall) row(1) lab(1 "Man") lab(2 "Services") lab(3 "FIRE") lab(4 "H & E") symy(1) symx(1))
graph export $postFiles/employmentChanges.pdf, replace

collapse (sum) ssOccDist*, by(broadSector)
gen sectoralGrowth = (log( ssOccDistT) - log(ssOccDist))*100

cap label drop sectorNames
label define sectorNames 1 "Manufacturing", add
label define sectorNames 2 "Services", modify
label define sectorNames 3 "FIRE", modify
label define sectorNames 4 "Health & Educ.", modify
label values broadSector sectorNames

graph hbar sectoralGrowth, over(broadSector) ytitle("") title("Change in Share of Workforce") nofill
graph export $postFiles/sectoralChanges.pdf, replace



******** TOTAL INCOME ********
insheet using $sandbox/totalIncomeDiff.csv, c clear

gen newcode = _n
reshape long v, i(newcode) j(period)
rename v diff

replace diff = diff*100

merge m:1 newcode using $workerFiles/newcodenames
drop _merge

keep newcode period diff disco2 broadSector
reshape wide diff, i(newcode disco2 broadSector) j(period)
sort newcode
save temp, replace

collapse (mean) diff*, by(broadSector)
gen newcode = broadSector + 38
append using temp
save temp, replace

keep if newcode<=38
collapse (sd) diff*, by(broadSector)
gen newcode = broadSector + 42
append using temp
save temp, replace

keep if newcode<=38
collapse (mean) diff*
gen newcode = 47
append using temp
save temp, replace

keep if newcode<=38
collapse (sd) diff*
gen newcode = 48
append using temp

sort broadSector newcode
order newcode broadSector

cap label drop codeNames
label define codeNames 1 "MANUFACTURING & Managers",add
label define codeNames 2 "&Science Professional",modify
label define codeNames 3 "&Science Assc. Professional",modify
label define codeNames 4 "&Other Assc. Professional",modify
label define codeNames 5 "&Clerks",modify
label define codeNames 6 "&Agriculture ",modify
label define codeNames 7 "&Building Trades",modify
label define codeNames 8 "&Metal Trades",modify
label define codeNames 9 "&Other Crafts",modify
label define codeNames 10 "&Plant Operator",modify
label define codeNames 11 "&Machine Operator",modify
label define codeNames 12 "&Drivers",modify
label define codeNames 13 "&Laborers",modify
label define codeNames 14 "SERVICES &Managers",modify
label define codeNames 15 "&Science Professional",modify
label define codeNames 16 "&Other Professional",modify
label define codeNames 17 "&Science Assc. Professional",modify
label define codeNames 18 "&Other Assc. Professional",modify
label define codeNames 19 "&Clerks",modify
label define codeNames 20 "&Personal Workers",modify
label define codeNames 21 "&Retail Workers",modify
label define codeNames 22 "&Metal Trades",modify
label define codeNames 23 "&Drivers",modify
label define codeNames 24 "&Elementary Occupations",modify
label define codeNames 25 "&Laborers",modify
label define codeNames 26 "FIRE &Managers",modify
label define codeNames 27 "&Science Professional",modify
label define codeNames 28 "&Other Professional",modify
label define codeNames 29 "&Science Assc. Professional",modify
label define codeNames 30 "&Other Assc. Professional",modify
label define codeNames 31 "&Clerks",modify
label define codeNames 32 "&Customer Service",modify
label define codeNames 33 "HEALTH \& EDUC. &Health Professional",modify
label define codeNames 34 "&Teachers",modify
label define codeNames 35 "&Health Assc. Professional",modify
label define codeNames 36 "&Teaching Assc. Professional",modify
label define codeNames 37 "&Clerks",modify
label define codeNames 38 "&Personal Workers",modify

label define codeNames 39 "&\bf Mean",modify
label define codeNames 40 "&\bf Mean",modify
label define codeNames 41 "&\bf Mean",modify
label define codeNames 42 "&\bf Mean",modify

label define codeNames 43 "&\bf SD",modify
label define codeNames 44 "&\bf SD",modify
label define codeNames 45 "&\bf SD",modify
label define codeNames 46 "&\bf SD",modify

label define codeNames 47 "&\bf GRAND MEAN",modify
label define codeNames 48 "&\bf GRAND SD",modify

label values newcode codeNames


* Choosing Rows
keep newcode diff2 diff9 diff19 diff29 diff40
gen amp = " & "
gen jump = " \\"

format diff* %4.2f

*gen hline = "\hline\\" if newcode==13| newcode==25 | newcode==32
gen hline = "\hline" if newcode>=43 & newcode<=46


format diff* %4.2f

outsheet newcode amp diff2 amp diff9 amp diff19 amp diff29 amp diff40 jump hline using $postFiles/wagesPlusCapital.tex, non noq  replace







******** REAL WAGES ********
insheet using $sandbox/meanRealWageDiffs.csv, c clear
gen newcode = _n
reshape long v, i(newcode) j(period)
rename v diff

replace diff = diff*100

merge m:1 newcode using $workerFiles/newcodenames
drop _merge

keep newcode period diff disco2 broadSector
reshape wide diff, i(newcode disco2 broadSector) j(period)
sort newcode
save temp, replace

collapse (mean) diff*, by(broadSector)
gen newcode = broadSector + 38
append using temp
save temp, replace

keep if newcode<=38
collapse (sd) diff*, by(broadSector)
gen newcode = broadSector + 42
append using temp
save temp, replace

keep if newcode<=38
collapse (mean) diff*
gen newcode = 47
append using temp
save temp, replace

keep if newcode<=38
collapse (sd) diff*
gen newcode = 48
append using temp

sort broadSector newcode
order newcode broadSector

cap label drop codeNames
label define codeNames 1 "MANUFACTURING & Managers",add
label define codeNames 2 "&Science Professional",modify
label define codeNames 3 "&Science Assc. Professional",modify
label define codeNames 4 "&Other Assc. Professional",modify
label define codeNames 5 "&Clerks",modify
label define codeNames 6 "&Agriculture ",modify
label define codeNames 7 "&Building Trades",modify
label define codeNames 8 "&Metal Trades",modify
label define codeNames 9 "&Other Crafts",modify
label define codeNames 10 "&Plant Operator",modify
label define codeNames 11 "&Machine Operator",modify
label define codeNames 12 "&Drivers",modify
label define codeNames 13 "&Laborers",modify
label define codeNames 14 "SERVICES &Managers",modify
label define codeNames 15 "&Science Professional",modify
label define codeNames 16 "&Other Professional",modify
label define codeNames 17 "&Science Assc. Professional",modify
label define codeNames 18 "&Other Assc. Professional",modify
label define codeNames 19 "&Clerks",modify
label define codeNames 20 "&Personal Workers",modify
label define codeNames 21 "&Retail Workers",modify
label define codeNames 22 "&Metal Trades",modify
label define codeNames 23 "&Drivers",modify
label define codeNames 24 "&Elementary Occupations",modify
label define codeNames 25 "&Laborers",modify
label define codeNames 26 "FIRE &Managers",modify
label define codeNames 27 "&Science Professional",modify
label define codeNames 28 "&Other Professional",modify
label define codeNames 29 "&Science Assc. Professional",modify
label define codeNames 30 "&Other Assc. Professional",modify
label define codeNames 31 "&Clerks",modify
label define codeNames 32 "&Customer Service",modify
label define codeNames 33 "HEALTH \& EDUC. &Health Professional",modify
label define codeNames 34 "&Teachers",modify
label define codeNames 35 "&Health Assc. Professional",modify
label define codeNames 36 "&Teaching Assc. Professional",modify
label define codeNames 37 "&Clerks",modify
label define codeNames 38 "&Personal Workers",modify

label define codeNames 39 "&\bf Mean",modify
label define codeNames 40 "&\bf Mean",modify
label define codeNames 41 "&\bf Mean",modify
label define codeNames 42 "&\bf Mean",modify

label define codeNames 43 "&\bf SD",modify
label define codeNames 44 "&\bf SD",modify
label define codeNames 45 "&\bf SD",modify
label define codeNames 46 "&\bf SD",modify

label define codeNames 47 "&\bf GRAND MEAN",modify
label define codeNames 48 "&\bf GRAND SD",modify

label values newcode codeNames


preserve

reshape long diff, i(newcode broadSector disco2) j(period)
keep if broadSector!=.
graph tw (kdensity diff if period==9, lc(gs0)) (kdensity diff if period==20, lc(gs8)) (kdensity diff if period==40, lc(gs13)), legend(lab(1 "0 Periods") lab(2 "10 Periods") lab(3 "30 Periods") si(small) row(1) symx(3) symy(2)) xtitle("w{subscript:Actual Prices} - w{subscript:Fixed Prices}") ytitle("Freq.") title("Density of Skill Price Differences") note("Differences are log points x 100""Personal Workers in H&E not shown for visibility.",si(vsmall))
graph export $postFiles/wageDensityProgression.pdf, replace
restore


preserve

reshape long diff, i(newcode broadSector disco2) j(period)
keep if broadSector==1
graph tw (kdensity diff if period==9, lc(gs0)) (kdensity diff if period==20, lc(gs8)) (kdensity diff if period==40, lc(gs13)), legend(lab(1 "0 Periods") lab(2 "10 Periods") lab(3 "30 Periods") si(small) row(1) symx(3) symy(2)) xtitle("w{subscript:Actual Prices} - w{subscript:Fixed Prices}") ytitle("Freq.") title("Density of Skill Price Differences") note("Differences are log points x 100",si(vsmall))
graph export $postFiles/wageDensityProgression_man.pdf, replace
restore

* Choosing Rows
keep newcode diff2 diff9 diff19 diff29 diff40
gen amp = " & "
gen jump = " \\"

format diff* %4.2f

*gen hline = "\hline\\" if newcode==13| newcode==25 | newcode==32
gen hline = "\hline" if newcode>=43 & newcode<=46


format diff* %4.2f

outsheet newcode amp diff2 amp diff9 amp diff19 amp diff29 amp diff40 jump hline using $postFiles/wages.tex, non noq  replace


