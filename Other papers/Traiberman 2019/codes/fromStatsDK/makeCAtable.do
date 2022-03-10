global doFiles "/Users/sharont/Dropbox/DanishLEEDProject/structuralEstimation_Final/Estimation/doFiles"
global dataPath "/Users/sharont/Dropbox/DanishLEEDProject/structuralEstimation_Final/Estimation/inputs"
global outPath "/Users/sharont/Dropbox/DanishLEEDProject/structuralEstimation_Final/Estimation/outputs"
global matlabPath "/Users/sharont/Dropbox/DanishLEEDProject/structuralEstimation_Final/Estimation/matlabFiles"

global figurePath "/Users/sharont/Dropbox/DanishLEEDProject/structuralEstimation_Final/Estimation/figures"
global tablePath "/Users/sharont/Dropbox/DanishLEEDProject/structuralEstimation_Final/Estimation/tables"
global workerFiles "/Users/sharont/Dropbox/DanishLEEDProject/structuralEstimation_Final/Calibration/workerFiles"


use $outPath/BETAMAT0, clear

keep newCode type betak

egen noHighEd = max((betak==0)*(type==5)), by(newCode)
egen noMedEd = max((betak==0)*(type==3)),by(newCode)
egen noLowEd = max((betak==0)*(type==1)),by(newCode)

replace betak = . if noHighEd & (type==5|type==6)
replace betak = . if noMedEd & (type==3|type==4)
replace betak = . if noLowEd & (type==1|type==2)

levelsof newCode if noHighEd, local(renorms)

gen thetak = exp(betak)
egen theta_1i = max(thetak*(newCode==1)),by(type)

gen caMeasure = thetak/theta_1i

keep newCode type caMeasure
reshape wide caMeasure, i(newCode) j(type)


foreach occ of local renorms{
	forval t = 1(1)4{
		replace caMeasure`t' = caMeasure`t'/caMeasure4 if newCode==`occ'
	}
}

rename newCode newcode
merge m:1 newcode using $workerFiles/newcodenames


cap label drop broadNames
label define broadNames 1 "Man.",add
label define broadNames 2 "Serv.",modify
label define broadNames 3 "FIRE",modify
label define broadNames 4 "H \& E",modify


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
label values broadSector broadNames

order broadSector disco2 caMeasure*

save $outPath/caMeasures, replace


use $outPath/caMeasures, clear
format caMeasure1-caMeasure5 %3.2f

tostring caMeasure*, replace force u
forval j = 1(1)6{	
	replace caMeasure`j' = "" if caMeasure`j' == "."
	replace caMeasure`j' = "1" if newcode==1
}

foreach occ of local renorms{
	replace caMeasure4 = "1" if newcode==`occ' & caMeasure4=="1.00"
}

sort broadSector disco2
by broadSector (disco2): gen first=_n==1

decode broadSector, gen(broadString)
decode disco2, gen(occString)

order broadString occString caMeasure*
replace broadString = "" if !first
sort broadSector disco2

gen amp = "&"
gen lb = "\\"
outsheet broadString amp occString amp caMeasure1 amp caMeasure2 amp caMeasure3 amp caMeasure4 amp caMeasure5 amp caMeasure6 lb using $tablePath/caTable.tex, noq replace non 


use $outPath/caMeasures, clear



use $outPath/BETAYEARS0, clear
collapse (mean) betay, by(newCode)
gsort -betay
gen rank = _n
rename newCode newcode
save tempRank, replace

use $outPath/BETAMAT0, clear

keep newCode type betak

egen noHighEd = max((betak==0)*(type==5)), by(newCode)
egen noMedEd = max((betak==0)*(type==3)),by(newCode)
egen noLowEd = max((betak==0)*(type==1)),by(newCode)

replace betak = . if noHighEd & (type==5|type==6)
replace betak = . if noMedEd & (type==3|type==4)
replace betak = . if noLowEd & (type==1|type==2)

levelsof newCode if noHighEd, local(renorms)

gen thetak = exp(betak)
egen theta_1i = max(thetak*(newCode==1)),by(type)

gen caMeasure = thetak/theta_1i



keep newCode type caMeasure betak
reshape wide caMeasure betak, i(newCode) j(type)


rename newCode newcode
merge m:1 newcode using $workerFiles/newcodenames
drop _merge

cap label drop broadNames
label define broadNames 1 "Man.",add
label define broadNames 2 "Serv.",modify
label define broadNames 3 "FIRE",modify
label define broadNames 4 "H \& E",modify


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
label values broadSector broadNames

/*
gsort -betak1 -betak2 -betak3 -betak4 -betak5
gen rank = _n
*/
merge 1:1 newcode using tempRank
drop _merge
*tsset rank

save caTempInfo, replace

corr betay caMeasure*

cap drop xPoint*
gsort -rank
forval j = 1(1)6{
	gen xPoint`j' = caMeasure`j'
	replace xPoint`j' = xPoint`j'[_n-1] if xPoint`j'[_n-1]!=.
}



tsset rank
forval j = 1(1)6{
replace xPoint`j' = L.xPoint`j' if xPoint`j'==.
*replace caMeasure`j' = L.caMeasure`j' if caMeasure`j'==.
}





forval j = 1(1)6{
	local xp`j'=xPoint`j'[38]
	display `xp`j''
}


graph tw (line caMeasure1 rank, lc(gs0) lw(.4))  ///
(line caMeasure2 rank, lc(gs0) lw(.4) lp(dash)) ///
(line caMeasure6 rank, lc(gs10) lw(.4) lp(dot)), ///
scheme(s1mono) legend(off) title("HS") ylabel(.5(.5)2,labs(small)) ///
xtitle("Mean Skill Price Rank of Occupation", si(small)) ///
ytitle("CA Index", si(small)) ///
text(`xp1' 38 "Type 1", si(vsmall) placement(north)) ///
text(`xp2' 38 "Type 2", si(vsmall) placement(north)) saving(graph1, replace)


graph tw (line caMeasure3 rank, lc(gs0) lw(.4))   ///
(line caMeasure4 rank, lc(gs0) lw(.4) lp(dash)) ///
(line caMeasure6 rank, lc(gs10) lw(.4) lp(dot)), ///
scheme(s1mono) legend(off) title("Voc.") ylabel(.5(.5)2,labs(small)) ///
xtitle("Mean Skill Price Rank of Occupation", si(small))  ///
ytitle("CA Index", si(small)) ///
text(`xp3' 38 "Type 3", si(vsmall) placement(north)) ///
text(`xp4' 38 "Type 4", si(vsmall) ) saving(graph2, replace)
*placement(south)

graph tw (line caMeasure5 rank, lc(gs0) lw(.4))  ///
(line caMeasure6 rank, lc(gs0) lw(.4) lp(dash)) ///
(line caMeasure6 rank, lc(gs10) lw(.4) lp(dot)), ///
scheme(s1mono) legend(off) title("BA+") ylabel(.5(.5)2,labs(small)) ///
xtitle("Mean Skill Price Rank of Occupation", si(small)) ///
ytitle("CA Index", si(small)) ///
text(`xp5' 38 "Type 5", si(vsmall) placement(south)) ///
text(`xp6' 38 "Type 6", si(vsmall) placement(south)) saving(graph3, replace)

gr combine graph1.gph graph2.gph graph3.gph, rows(2) scheme(s1mono) note("Skill prices averaged over time and ranked highest to lowest", si(tiny))
gr export $figurePath/caFigure.eps, replace


use caTempInfo, clear

corr betay caMeasure*
matrix bob = r(C)
clear
svmat bob

drop if _n==1
gen type = _n
keep bob1 type
order type bob1
save compCorrs, replace

use $outPath/BETAMAT0, clear

keep newCode type betak
gen thetak = exp(betak)
collapse (mean) thetak betak, by(type)

cap label drop typeNames
label define typeNames 1 "1 & HS, L",add
label define typeNames 2 "2 & HS, H",modify
label define typeNames 3 "3 & Voc, L",modify
label define typeNames 4 "4 & Voc, H",modify
label define typeNames 5 "5 & Col, L",modify
label define typeNames 6 "6 & Col, H",modify

label values type typeNames
decode type, gen(typeString)

* add income information... (from Stats DK-should merge properly)
gen relIncome = .
replace relIncome = .330 if type==1
replace relIncome = .823 if type==2
replace relIncome = .539 if type==3
replace relIncome = .921 if type==4
replace relIncome = .576 if type==5
replace relIncome = 1.000 if type==6

merge 1:1 type using compCorrs
drop _merge

*estpost tabstat thetak relIncome, by(type) nototal
*esttab using $tablePath/absoluteAdvantageTable.tex, cell(thetak(fmt(3)) relIncome(fmt(3))) noobs varlabels(`e(labels)') compress frag booktabs plain mlabels(none) collabels(none) replace
eststo clear

estpost tabstat thetak relIncome bob1, by(type) nototal statistics(mean)
*esttab using $tablePath/absoluteAdvantageTable.tex, main(mean) noobs varlabels(`e(labels)') compress frag booktabs plain mlabels(none) collabels(none) replace
esttab using $tablePath/absoluteAdvantageTable.tex, cell("thetak(fmt(3)) relIncome(fmt(3)) bob1(fmt(3))") noobs varlabels(`e(labels)') compress frag booktabs plain mlabels(none) collabels(none) replace





/*
graph tw ///
(line caMeasure1 rank, lc(gs10) lw(.2) lp(_-#)) (lowess caMeasure1 rank, lc(gs0) lw(.4) lp(_-#))  ///
(line caMeasure2 rank, lc(gs10) lw(.2) lp(--)) (lowess caMeasure2 rank, lc(gs0) lw(.4) lp(--)) ///
(line caMeasure3 rank, lc(gs10) lw(.2) lp(_--_#)) (lowess caMeasure3 rank, lc(gs0) lw(.4) lp(_--_#)) ///
(line caMeasure4 rank, lc(gs10) lw(.2) lp(.--.)) (lowess caMeasure4 rank, lc(gs0) lw(.4) lp(.--.)) ///
(line caMeasure5 rank, lc(gs10) lw(.2) lp(.-_.)) (lowess caMeasure5 rank, lc(gs0) lw(.4) lp(.-_.)) ///
(line caMeasure6 rank, lc(gs10) lw(.2) lp(solid)) (lowess caMeasure6 rank, lc(gs0) lw(.4) lp(solid)) , ///
scheme(s1mono) legend(order(2 4 6 8 10 12) row(2) lab(2 "Type 1") lab(4 "Type 2") lab(6 "Type 3") lab(8 "Type 4") lab(10 "Type 5") lab(12 "Type 6"))


tsset rank
forval j = 1(1)6{
	lowess caMeasure`j' rank, gen(sca`j')
}

gsort -rank
forval j = 1(1)6{
	gen xPoint`j' = sca`j'
	replace xPoint`j' = xPoint`j'[_n-1] if xPoint`j'[_n-1]!=.
}

	

/*
forval j = 1(1)6{
gen xPoint`j' = caMeasure`j'
replace xPoint`j' = xPoint`j'[_n-1] if xPoint`j'[_n-1]!=.
}
*/


tsset rank
forval j = 1(1)6{
replace xPoint`j' = L.xPoint`j' if xPoint`j'==.
}





forval j = 1(1)6{
	local xp`j'=xPoint`j'[38]
	display `xp`j''
}

* Order = 3 5 2 4 6 1
graph tw ///
(line caMeasure1 rank, lc(gs10) lw(.2) lp(dash)) (lowess caMeasure1 rank, lc(gs0) lw(.4) lp(dash))  ///
(line caMeasure2 rank, lc(gs10) lw(.2) lp(dash)) (lowess caMeasure2 rank, lc(gs0) lw(.4) lp(dash)) ///
(line caMeasure3 rank, lc(gs10) lw(.2) lp(__..)) (lowess caMeasure3 rank, lc(gs0) lw(.4) lp(__..)) ///
(line caMeasure4 rank, lc(gs10) lw(.2) lp(__..)) (lowess caMeasure4 rank, lc(gs0) lw(.4) lp(__..)) ///
(line caMeasure5 rank, lc(gs10) lw(.2)) (lowess caMeasure5 rank, lc(gs0) lw(.4)) ///
(line caMeasure6 rank, lc(gs10) lw(.2)) (lowess caMeasure6 rank, lc(gs0) lw(.4)) , ///
scheme(s1mono) legend(order(2 6 10)row(1) lab(2 "Short Ed") lab(6 "Med Ed") lab(10 "BA+"))  ///
xtitle("Occupation Rank") ///
ytitle("Comparative Advantage Index") ///
text(`xp1' 38 "Type 1", si(vsmall) placement(neast)) ///
text(`xp2' 38 "Type 2", si(vsmall) placement(east)) ///
text(`xp3' 38 "Type 3", si(vsmall) placement(seast)) ///
text(`xp4' 38 "Type 4", si(vsmall) placement(east)) ///
text(`xp5' 38 "Type 5", si(vsmall) placement(neast)) ///
text(`xp6' 38 "Type 6", si(vsmall) placement(east))




forval j = 1(1)6{
	local xp`j'=xPoint`j'[38]
	display `xp`j''
}


graph tw (line caMeasure1 rank, lc(gs10) lw(.2)) (lowess caMeasure1 rank, lc(gs0) lw(.4))  ///
(line caMeasure2 rank, lc(gs10) lw(.2)) (lowess caMeasure2 rank, lc(gs0) lw(.4)), ///
scheme(s1mono) legend(off) title("HS") ylabel(.5(.5)2) ///
xtitle("Occupation Rank") ///
ytitle("Comparative Advantage Index") ///
text(`xp1' 38 "Type 1", si(vsmall) placement(west)) ///
text(`xp2' 38 "Type 2", si(vsmall) placement(nwest)) saving(graph1, replace)


graph tw (line caMeasure3 rank, lc(gs10) lw(.2)) (lowess caMeasure3 rank, lc(gs0) lw(.4))  ///
(line caMeasure4 rank, lc(gs10) lw(.2)) (lowess caMeasure4 rank, lc(gs0) lw(.4)), ///
scheme(s1mono) legend(off) title("Voc.") ylabel(.5(.5)2) ///
xtitle("Occupation Rank") ///
ytitle("Comparative Advantage Index") ///
text(`xp3' 38 "Type 3", si(vsmall) placement(west)) ///
text(`xp4' 38 "Type 4", si(vsmall) placement(swest)) saving(graph2, replace)


graph tw (line caMeasure5 rank, lc(gs10) lw(.2)) (lowess caMeasure5 rank, lc(gs0) lw(.4))  ///
(line caMeasure6 rank, lc(gs10) lw(.2)) (lowess caMeasure6 rank, lc(gs0) lw(.4)), ///
scheme(s1mono) legend(off) title("BA+") ylabel(.5(.5)2) xlabel(#5) ///
xtitle("Occupation Rank") ///
ytitle("Comparative Advantage Index") ///
text(`xp5' 38 "Type 5", si(vsmall) placement(nwest)) ///
text(`xp6' 38 "Type 6", si(vsmall) placement(swest)) saving(graph3, replace)

gr combine graph1.gph graph2.gph graph3.gph, rows(2)

*legend(order(2 4 6 8 10 12) row(2) lab(2 "Type 1") lab(4 "Type 2") lab(6 "Type 3") lab(8 "Type 4") lab(10 "Type 5") lab(12 "Type 6"))
*/





