global figurePath "/Users/SharonT/Dropbox/DanishLEEDProject/structuralEstimation_Final/figures/"
global tablePath "/Users/sharont/Dropbox/DanishLEEDProject/structuralEstimation_Final/Estimation/tables"

insheet using historyCF_OFF.csv, c clear
drop v3
rename v# v#,addnumber(1)

rename v1 pid
rename v2 t
rename v3 age
rename v4 type
rename v5 tenure
rename v6 occCF
rename v7 employedCF
rename v8 switchCF
rename v9 wageCF
compress
save historyCF_OFF, replace

insheet using historyT_OFF.csv, c clear
drop v3
rename v# v#,addnumber(1)

rename v1 pid
rename v2 t
rename v3 age
rename v4 type
rename v5 tenure
rename v6 occT
rename v7 employedT
rename v8 switchT
rename v9 wageT
compress
save historyT, replace


use historyT, clear
merge 1:1 pid t using historyCF_OFF
drop _merge

save mergedHistories_OFF, replace





use mergedHistories_OFF, clear

rename occT newcode
merge m:1 newcode using eta
drop _merge
replace eta = 0 if employedT==0

rename eta etaT
rename newcode occT


rename occCF newcode
merge m:1 newcode using eta
replace eta = 0 if employedCF==0
drop _merge
rename eta etaCF
rename newcode occCF


drop if t <=1

gen lifetimeT = .96^(t-2)*wageT
gen lifetimeCF = .96^(t-2)*wageCF
gen lifetime2T = .96^(t-2)*(wageT+etaT)
gen lifetime2CF = .96^(t-2)*(wageCF+etaCF)

collapse (sum) lifetime* (min) age, by(pid type)

gen logDiff1 = (log(lifetimeT) - log(lifetimeCF))*100
gen logDiff2 = (log(lifetime2T) - log(lifetime2CF))*100
drop if logDiff1==.

gen dip1 = (logDiff1<0)*100
gen dip2 = (logDiff2<0)*100

cap label drop skillNames
label define skillNames 1 "Low [L]", modify
label define skillNames 2 "Low [H]", modify
label define skillNames 3 "Med [L]", modify
label define skillNames 4 "Med [H]", modify
label define skillNames 5 "High [L]", modify
label define skillNames 6 "High [H]", modify
label values type skillNames


eststo diff1Dist: estpost tabstat logDiff1, by(type) s(mean sd p5 p50 p95) columns(statistics)
eststo fracDip1: estpost tabstat dip1, by(type) s(mean) columns(statistics)


esttab diff1Dist fracDip1 using $tablePath/diff1Dist_OFF.tex,tex replace cells("mean(fmt(%4.2f)) sd(fmt(%4.2f) pattern(1 0)) p5(fmt(%4.2f) pattern(1 0)) p50(fmt(%4.2f) pattern(1 0)) p95(fmt(%4.2f) pattern(1 0))") varlabels(`e(labels)') fragment plain noobs non mlabels(none) collabels(none)


eststo diff2Dist: estpost tabstat logDiff2, by(type) s(mean sd p5 p50 p95) columns(statistics) 
eststo fracDip2: estpost tabstat dip2, by(type) s(mean) columns(statistics)

esttab diff2Dist fracDip2 using $tablePath/diff2Dist_OFF.tex,tex replace cells("mean(fmt(%4.2f)) sd(fmt(%4.2f) pattern(1 0)) p5(fmt(%4.2f) pattern(1 0)) p50(fmt(%4.2f) pattern(1 0)) p95(fmt(%4.2f) pattern(1 0))") varlabels(`e(labels)') fragment plain noobs non mlabels(none) collabels(none)




use mergedHistories_OFF, clear

rename occT newcode
merge m:1 newcode using eta
drop _merge
replace eta = 0 if employedT==0

rename eta etaT
rename newcode occT


rename occCF newcode
merge m:1 newcode using eta
replace eta = 0 if employedCF==0
drop _merge
rename eta etaCF
rename newcode occCF



drop if t <=1

gen lifetimeT = .96^(t-2)*wageT
gen lifetimeCF = .96^(t-2)*wageCF
gen lifetime2T = .96^(t-2)*(wageT+etaT)
gen lifetime2CF = .96^(t-2)*(wageCF+etaCF)

collapse (sum) lifetime* (min) age, by(pid type)


gen ageGroup = floor((age-1)/5)
collapse (mean) lifetimeT lifetimeCF, by(ageGroup type)
gen logDiff = log(lifetimeT)-log(lifetimeCF)

cap label drop typeNames
label define typeNames 1 "HS, Low AA", modify
label define typeNames 2 "HS, High AA", modify
label define typeNames 3 "Voc, Low AA", modify
label define typeNames 4 "Voc, High AA", modify
label define typeNames 5 "BA+, Low AA", modify
label define typeNames 6 "BA+, High AA", modify
label values type typeNames 



graph bar logDiff, over(age) by(type, note("Ages binned into 5 year groups.")) ytitle("Log Difference in Lifetime Earnings") scheme(s1mono)
graph export $figurePath/cfOutcomesByAge_OFF.eps, replace
