libname painel "C:\Users\Dix-Carneiro\Documents\Trade and Labor Dynamics\Data for Estimation New\Codes_7_Sectors\PanelRAIS" ;
libname reg "C:\Users\Dix-Carneiro\Documents\Trade and Labor Dynamics\Data for Estimation New\Codes_7_Sectors\PanelRAIS" ;

data new_painel (keep = pis educ_new idade_new sexo_new ano rwage_corr sector 
                        horas_contr choice0 choice1 choice2 choice3 choice4 choice5 choice6 choice7
                        log_w hourly_wage) ;
	set painel.panel9505 ;

	choice0 = 0 ;
	choice1 = 0 ;
	choice2 = 0 ;
	choice3 = 0 ;
	choice4 = 0 ;
	choice5 = 0 ;
	choice6 = 0 ;
	choice7 = 0 ;

	if sector = 0 then choice0 = 1 ;
	if sector = 1 then choice1 = 1 ;
	if sector = 2 then choice2 = 1 ;
	if sector = 3 then choice3 = 1 ;
	if sector = 4 then choice4 = 1 ;
	if sector = 5 then choice5 = 1 ;
	if sector = 6 then choice6 = 1 ;
	if sector = 7 then choice7 = 1 ;

	if (horas_contr ~= . and horas_contr) > 0 then hourly_wage = rwage_corr / (4.35*horas_contr) ;
	else hourly_wage = . ;

	log_w = log(hourly_wage) ;

run ;



data new_painel_1995 ;
	set new_painel ;

	if ano <= 1995 - 1 and ano >= 1995 - 9 ;

run ;



proc means data = new_painel_1995 noprint ;
	var choice1 choice2 choice3 choice4 choice5 choice6 choice7 ;
	class pis ;
	output out = Experience (drop = _type_ _freq_) sum = exper1 exper2 exper3 exper4 exper5 exper6 exper7 ;
run ;

data Experience ;
	set Experience ;

	ano = 1995 ; 

run ;



proc delete data = new_painel_1995 ; run ;



%macro macro1 ;

%do year = 1996 %to 2005 ;

data new_painel_&year. ;
	set new_painel ;

	if ano <= &year. - 1 and ano >= &year. - 9 ;

run ;

proc means data = new_painel_&year. noprint ;
	var choice1 choice2 choice3 choice4 choice5 choice6 choice7 ;
	class pis ;
	output out = Experience_&year. (drop = _type_ _freq_) sum = exper1 exper2 exper3 exper4 exper5 exper6 exper7 ;
run ;

data Experience_&year. ;
	set Experience_&year. ;

	ano = &year. ;

run ;

data experience ;
	set experience Experience_&year. ;
run ;

proc delete data = new_painel_&year. ; run ;
proc delete data = Experience_&year. ; run ;

%end ;

%mend ;

%macro1 ;



proc sort data = experience ; by pis ; run ;



data experience ;
	set experience ;

	if compress(pis) ~= "" ;

run ;




data painel_lag_sector (keep = pis ano lag_sector lag_log_w) ;
	set new_painel ;

	ano             = ano + 1 ;
	lag_sector      = sector ;
	lag_log_w       = log_w ;

run ;


data painel_lag2_sector (keep = pis ano lag2_sector) ;
	set new_painel ;

	ano         = ano + 2 ;
	lag2_sector = sector ;

run ;



proc sort data = experience          ; by pis ano ; run ;
proc sort data = new_painel          ; by pis ano ; run ;
proc sort data = painel_lag_sector   ; by pis ano ; run ;
proc sort data = painel_lag2_sector  ; by pis ano ; run ;


data new_painel (drop = sexo_new horas_contr) ;
	merge new_painel (in = a) experience painel_lag_sector painel_lag2_sector ;
	by pis ano ;
	if a ;

	if ano >= 1995 ;

    gender = (sexo_new = 2) ;


    dummy_educ2 = (educ_new = 2) ;
	dummy_educ3 = (educ_new = 3) ;
	dummy_educ4 = (educ_new = 4) ;

	educ        = (educ_new >= 3) ;


	age   = idade_new - 25 ;
	age_2 = (idade_new - 25)**2 ;	


	if idade_new >= 25 and idade_new <= 60 ;

	if (gender = . or educ_new = . or age = .) then delete ;

	year95 = (ano = 1995) ;
	year96 = (ano = 1996) ;
	year97 = (ano = 1997) ;
	year98 = (ano = 1998) ;
	year99 = (ano = 1999) ;
	year00 = (ano = 2000) ;
	year01 = (ano = 2001) ;
	year02 = (ano = 2002) ;
	year03 = (ano = 2003) ;
	year04 = (ano = 2004) ;
	year05 = (ano = 2005) ;

	sector0 = (sector = 0) ;
	sector1 = (sector = 1) ;
	sector2 = (sector = 2) ;
	sector3 = (sector = 3) ;
	sector4 = (sector = 4) ;
	sector5 = (sector = 5) ;
	sector6 = (sector = 6) ;
	sector7 = (sector = 7) ;

	lag_sector0 = (lag_sector = 0) ;
	lag_sector1 = (lag_sector = 1) ;
	lag_sector2 = (lag_sector = 2) ;
	lag_sector3 = (lag_sector = 3) ;
	lag_sector4 = (lag_sector = 4) ;
	lag_sector5 = (lag_sector = 5) ;
	lag_sector6 = (lag_sector = 6) ;
	lag_sector7 = (lag_sector = 7) ;

	sector00 = (lag_sector = 0 and sector = 0) ;
	sector01 = (lag_sector = 0 and sector = 1) ;
	sector02 = (lag_sector = 0 and sector = 2) ;
	sector03 = (lag_sector = 0 and sector = 3) ;
	sector04 = (lag_sector = 0 and sector = 4) ;
	sector05 = (lag_sector = 0 and sector = 5) ;
	sector06 = (lag_sector = 0 and sector = 6) ;
	sector07 = (lag_sector = 0 and sector = 7) ;

	sector10 = (lag_sector = 1 and sector = 0) ;
	sector11 = (lag_sector = 1 and sector = 1) ;
	sector12 = (lag_sector = 1 and sector = 2) ;
	sector13 = (lag_sector = 1 and sector = 3) ;
	sector14 = (lag_sector = 1 and sector = 4) ;
	sector15 = (lag_sector = 1 and sector = 5) ;
	sector16 = (lag_sector = 1 and sector = 6) ;
	sector17 = (lag_sector = 1 and sector = 7) ;

	sector20 = (lag_sector = 2 and sector = 0) ;
	sector21 = (lag_sector = 2 and sector = 1) ;
	sector22 = (lag_sector = 2 and sector = 2) ;
	sector23 = (lag_sector = 2 and sector = 3) ;
	sector24 = (lag_sector = 2 and sector = 4) ;
	sector25 = (lag_sector = 2 and sector = 5) ;
	sector26 = (lag_sector = 2 and sector = 6) ;
	sector27 = (lag_sector = 2 and sector = 7) ;

	sector30 = (lag_sector = 3 and sector = 0) ;
	sector31 = (lag_sector = 3 and sector = 1) ;
	sector32 = (lag_sector = 3 and sector = 2) ;
	sector33 = (lag_sector = 3 and sector = 3) ;
	sector34 = (lag_sector = 3 and sector = 4) ;
	sector35 = (lag_sector = 3 and sector = 5) ;
	sector36 = (lag_sector = 3 and sector = 6) ;
	sector37 = (lag_sector = 3 and sector = 7) ;

	sector40 = (lag_sector = 4 and sector = 0) ;
	sector41 = (lag_sector = 4 and sector = 1) ;
	sector42 = (lag_sector = 4 and sector = 2) ;
	sector43 = (lag_sector = 4 and sector = 3) ;
	sector44 = (lag_sector = 4 and sector = 4) ;
	sector45 = (lag_sector = 4 and sector = 5) ;
	sector46 = (lag_sector = 4 and sector = 6) ;
	sector47 = (lag_sector = 4 and sector = 7) ;

	sector50 = (lag_sector = 5 and sector = 0) ;
	sector51 = (lag_sector = 5 and sector = 1) ;
	sector52 = (lag_sector = 5 and sector = 2) ;
	sector53 = (lag_sector = 5 and sector = 3) ;
	sector54 = (lag_sector = 5 and sector = 4) ;
	sector55 = (lag_sector = 5 and sector = 5) ;
	sector56 = (lag_sector = 5 and sector = 6) ;
	sector57 = (lag_sector = 5 and sector = 7) ;

	sector60 = (lag_sector = 6 and sector = 0) ;
	sector61 = (lag_sector = 6 and sector = 1) ;
	sector62 = (lag_sector = 6 and sector = 2) ;
	sector63 = (lag_sector = 6 and sector = 3) ;
	sector64 = (lag_sector = 6 and sector = 4) ;
	sector65 = (lag_sector = 6 and sector = 5) ;
	sector66 = (lag_sector = 6 and sector = 6) ;
	sector67 = (lag_sector = 6 and sector = 7) ;

	sector70 = (lag_sector = 7 and sector = 0) ;
	sector71 = (lag_sector = 7 and sector = 1) ;
	sector72 = (lag_sector = 7 and sector = 2) ;
	sector73 = (lag_sector = 7 and sector = 3) ;
	sector74 = (lag_sector = 7 and sector = 4) ;
	sector75 = (lag_sector = 7 and sector = 5) ;
	sector76 = (lag_sector = 7 and sector = 6) ;
	sector77 = (lag_sector = 7 and sector = 7) ;

	wage_dif = . ;
	if (ano > 1995 & lag_sector ~= . & sector ~= . & lag_sector = sector & idade_new >= 26) then 
    wage_dif = log_w - lag_log_w ;

	born = ano - idade_new ;
	generation = 0 ;
	if (born >= 1995 - 32 and  born <= 1995 - 28) then generation = 1 ;
	if (born >= 1995 - 47 and  born <= 1995 - 43) then generation = 2 ;

	pre_exper4 = exper4 + exper5 + exper6 + exper7 ;

run ;

data painel.new_painel (keep = pis ano gender educ_new idade_new sector lag_sector lag2_sector log_w 
                               exper1 exper2 exper3 exper4 exper5 exper6 exper7 
                               year95 year96 year97 year98 year99 year00 year01
                               year02 year03 year04 year05) ;
	set new_painel ;
run ;


data choices (keep = pis ano sector) ;
	set new_painel ;
run ;

proc transpose data = choices out = choices_wide prefix = sector ;
    by pis ;
    id ano ;
    var sector ;
run;


data frequencies (keep = pis freq0 freq1 freq2 freq3 freq4 freq5 freq6 freq7 switches) ;
	set choices_wide ;

	freq0 = 0 ;
	freq1 = 0 ;
	freq2 = 0 ;
	freq3 = 0 ;
	freq4 = 0 ;
	freq5 = 0 ;
	freq6 = 0 ;
	freq7 = 0 ;
	switches = 0 ;

	array vars {11} sector1995 sector1996 sector1997 sector1998 sector1999 
                    sector2000 sector2001 sector2002 sector2003 sector2004 sector2005 ;
	do i = 1 to 11 ;
		if vars(i) = 0 then freq0 = freq0 + 1 ;
		if vars(i) = 1 then freq1 = freq1 + 1 ;
		if vars(i) = 2 then freq2 = freq2 + 1 ;
		if vars(i) = 3 then freq3 = freq3 + 1 ;
		if vars(i) = 4 then freq4 = freq4 + 1 ;
		if vars(i) = 5 then freq5 = freq5 + 1 ;
		if vars(i) = 6 then freq6 = freq6 + 1 ;
		if vars(i) = 7 then freq7 = freq7 + 1 ;
	end ;

	do i = 2 to 11 ;
		if vars(i) ~= vars(i-1) then switches = switches + 1 ;
	end ;

run ;


data choices_wide (keep = pis sector0_1998 sector1_1998 sector2_1998 sector3_1998 sector4_1998 sector5_1998 sector6_1998 sector7_1998
                              sector0_2000 sector1_2000 sector2_2000 sector3_2000 sector4_2000 sector5_2000 sector6_2000 sector7_2000
                              sector0_2005 sector1_2005 sector2_2005 sector3_2005 sector4_2005 sector5_2005 sector6_2005 sector7_2005) ;
	set choices_wide ;

	sector0_1998 = (sector1998 = 0) ;
	sector1_1998 = (sector1998 = 1) ;
	sector2_1998 = (sector1998 = 2) ;
	sector3_1998 = (sector1998 = 3) ;
	sector4_1998 = (sector1998 = 4) ;
	sector5_1998 = (sector1998 = 5) ;
	sector6_1998 = (sector1998 = 6) ;
	sector7_1998 = (sector1998 = 7) ;

	sector0_2000 = (sector2000 = 0) ;
	sector1_2000 = (sector2000 = 1) ;
	sector2_2000 = (sector2000 = 2) ;
	sector3_2000 = (sector2000 = 3) ;
	sector4_2000 = (sector2000 = 4) ;
	sector5_2000 = (sector2000 = 5) ;
	sector6_2000 = (sector2000 = 6) ;
	sector7_2000 = (sector2000 = 7) ;

	sector0_2005 = (sector2005 = 0) ;
	sector1_2005 = (sector2005 = 1) ;
	sector2_2005 = (sector2005 = 2) ;
	sector3_2005 = (sector2005 = 3) ;
	sector4_2005 = (sector2005 = 4) ;
	sector5_2005 = (sector2005 = 5) ;
	sector6_2005 = (sector2005 = 6) ;
	sector7_2005 = (sector2005 = 7) ;

run ;


proc sort nodupkey data = choices_wide ; by pis ; run ;
proc sort nodupkey data = frequencies  ; by pis ; run ;
proc sort data = new_painel ; by pis ; run ;


data new_painel ;
	merge new_painel (in = a) choices_wide frequencies ;
	by pis ;
	if a ;
run ;



proc delete data = Painel_lag_sector Experience ; run ;


proc means data = new_painel median mean ;
	var hourly_wage ;
run ;



/*Simple Means*/

proc means data = new_painel ;
	var hourly_wage log_w ;
	class sector ;
    output out = reg.WageMeans (drop = _type_ _freq_)
    mean = hourly_wage log_w ;
run ;

proc means data = new_painel ;
	var choice0 choice1 choice2 choice3 choice4 choice5 choice6 choice7 ;
    output out = reg.EmpMeans (drop = _type_ _freq_)
    mean = choice0 choice1 choice2 choice3 choice4 choice5 choice6 choice7 ;
run ;

proc tabulate data = new_painel out = reg.Transition (drop = _type_ _page_ _table_ rename = PctN_10 = transition);
	class lag_sector sector ;
	tables lag_sector, sector*(rowpctn) /printmiss ;
run ;



/*Total Hourly Wages by Sector*/

proc means data = new_painel (keep = ano sector educ_new hourly_wage rwage_corr) noprint ;
	var hourly_wage	rwage_corr ;
	class ano sector educ_new ;
	output out = reg.total_hr_wage (drop = _type_ _freq_)
	sum = hourly_wage rwage_corr
    n   = cellsize1 cellsize2 ;
run ;


data reg.total_hr_wage ;
	set	reg.total_hr_wage ;

	if (ano = . or sector = . or educ_new = . or sector = 0) then delete ;

run ;


%macro Return_reg ;

	%do s = 0 %to 7 ;

		proc reg data = new_painel noprint OUTSSCP = sscp outest = reg.Return&s. 
			(drop = _MODEL_ _TYPE_ _NAME_ _DEPVAR_ sector&s.) covout ;
			where lag2_sector = &s. and lag_sector ~= &s. and lag_sector ~= . and idade_new >= 25 and idade_new <= 60 ;
			model sector&s. = year95 year96 year97 year98 year99 year00 year01 year02 year03 year04 year05
                              gender dummy_educ2 dummy_educ3 dummy_educ4 age age_2 
                              exper1 exper2 exper3 exper4 exper5 exper6 exper7 /noint ;
		run ;
		quit ;

		data reg.CPReturn&s. (keep = Nobs) ;
			set sscp ;
			if _TYPE_ = "N" ;
			Nobs = Intercept ;
		run ;

		data reg.CPReturn&s. (keep = Nobs year95 year96 year97 year98 year99 year00 year01 year02 year03 year04 year05
    				            gender dummy_educ2 dummy_educ3 dummy_educ4 age age_2 
                                exper1 exper2 exper3 exper4 exper5 exper6 exper7 _NAME_) ;
			if _N_ = 1 then set reg.CPReturn&s. ;
			set sscp ;

			if _NAME_ = "Intercept" or _NAME_ = "sector&s." or _NAME_ = "" then delete ;

			year95       = year95       / Nobs ;
			year96       = year96       / Nobs ;
			year97       = year97       / Nobs ;
			year98       = year98       / Nobs ;
			year99       = year99       / Nobs ;
			year00       = year00       / Nobs ;
			year01       = year01       / Nobs ;
			year02       = year02       / Nobs ;
			year03       = year03       / Nobs ;
			year04       = year04       / Nobs ;
			year05       = year05       / Nobs ;
			gender       = gender       / Nobs ;
			dummy_educ2  = dummy_educ2  / Nobs ;
			dummy_educ3  = dummy_educ3  / Nobs ;
			dummy_educ4  = dummy_educ4  / Nobs ;
			age          = age          / Nobs ;
			age_2        = age_2        / Nobs ;
			exper1       = exper1       / Nobs ;
			exper2       = exper2       / Nobs ;
			exper3       = exper3       / Nobs ;
			exper4       = exper4       / Nobs ;
			exper5       = exper5       / Nobs ;
			exper6       = exper6       / Nobs ;
			exper7       = exper7       / Nobs ;

		run ;

	%end ;

%mend ;

%Return_reg ;


%macro Wage_reg ;

	%do s = 1 %to 7 ;

/*Simple Wage Regressions, with no lag sector*/

		proc reg data = new_painel noprint OUTSSCP = sscp outest = reg.WageSector&s. 
		(drop = _MODEL_ _TYPE_ _NAME_ _DEPVAR_ log_w) covout ;
		where sector = &s. and idade_new >= 25 and idade_new <= 60 ;
		model log_w = year95 year96 year97 year98 year99 year00 year01 year02 year03 year04 year05
		              gender dummy_educ2 dummy_educ3 dummy_educ4 age age_2 
                      exper1 exper2 exper3 exper4 exper5 exper6 exper7 / noint ;
		OUTPUT OUT = residuals RESIDUAL = resid ;
		run ;
		quit ;


		data reg.CPWageSector&s. (keep = Nobs) ;
			set sscp ;
			if _TYPE_ = "N" ;
			Nobs = Intercept ;
		run ;


        data reg.CPWageSector&s. (keep = Nobs year95 year96 year97 year98 year99 year00 year01 year02 year03 year04 year05
		                          gender dummy_educ2 dummy_educ3 dummy_educ4 age age_2 
                                  exper1 exper2 exper3 exper4 exper5 exper6 exper7 _NAME_) ;
			if _N_ = 1 then set reg.CPWageSector&s. ;
			set sscp ;

			if _NAME_ = "Intercept" or _NAME_ = "log_w" or _NAME_ = "" then delete ;

			year95               = year95               / Nobs ;
			year96               = year96               / Nobs ;
			year97               = year97               / Nobs ;
			year98               = year98               / Nobs ;
			year99               = year99               / Nobs ;
			year00               = year00               / Nobs ;
			year01               = year01               / Nobs ;
			year02               = year02               / Nobs ;
			year03               = year03               / Nobs ;
			year04               = year04               / Nobs ;
			year05               = year05               / Nobs ;
			gender               = gender               / Nobs ;
			dummy_educ2          = dummy_educ2          / Nobs ;
			dummy_educ3          = dummy_educ3          / Nobs ;
			dummy_educ4          = dummy_educ4          / Nobs ;
			age                  = age                  / Nobs ;
			age_2                = age_2                / Nobs ;
			exper1               = exper1               / Nobs ;
			exper2               = exper2               / Nobs ;
			exper3               = exper3               / Nobs ;
			exper4               = exper4               / Nobs ;
			exper5               = exper5               / Nobs ;
			exper6               = exper6               / Nobs ;
			exper7               = exper7               / Nobs ;
	
		run ;


/*Variance of Variance Residuals*/

		data residuals ;
			set residuals ;
	
			if sector = &s. and idade_new >= 25 and idade_new <= 60 ;
			resid_2 = resid**2 ;

		run ;       


		proc means data = residuals noprint ;
			class sector ;
			var resid_2 ;
			output out = var_residuals&s. (drop = _type_ _freq_)
			var = var_resid 
    		n   = nobs ;
		run ;

		data var_residuals&s. ;
			retain sector ;
			set var_residuals&s. ;

			if sector = . then delete ;

		run ;


/*Wage regressions with lag sectoral choice dummies        */

	    proc reg data = new_painel noprint OUTSSCP = sscp outest = reg.WageSector&s._lag 
		(drop = _MODEL_ _TYPE_ _NAME_ _DEPVAR_ log_w) covout ;
		where sector = &s. and idade_new >= 25 and idade_new <= 60 ;
		model log_w = year95 year96 year97 year98 year99 year00 year01 year02 year03 year04 year05
		              lag_sector1 lag_sector2 lag_sector3 lag_sector4 lag_sector5 lag_sector6 lag_sector7
    		          gender dummy_educ2 dummy_educ3 dummy_educ4 age age_2 
                      exper1 exper2 exper3 exper4 exper5 exper6 exper7 / noint ;
		OUTPUT OUT = residuals RESIDUAL = resid ;
		run ;
		quit ;


		data reg.CPWageSector&s._lag (keep = Nobs) ;
			set sscp ;
			if _TYPE_ = "N" ;
			Nobs = Intercept ;
		run ;


        data reg.CPWageSector&s._lag (keep = Nobs year95 year96 year97 year98 year99 year00 year01 year02 year03 year04 year05
		                          lag_sector1 lag_sector2 lag_sector3 lag_sector4 lag_sector5 lag_sector6 lag_sector7
              					  gender dummy_educ2 dummy_educ3 dummy_educ4 age age_2 exper1 exper2 exper3 exper4 exper5 exper6 exper7 _NAME_) ;
			if _N_ = 1 then set reg.CPWageSector&s._lag ;
			set sscp ;

			if _NAME_ = "Intercept" or _NAME_ = "log_w" or _NAME_ = "" then delete ;

			year95               = year95               / Nobs ;
			year96               = year96               / Nobs ;
			year97               = year97               / Nobs ;
			year98               = year98               / Nobs ;
			year99               = year99               / Nobs ;
			year00               = year00               / Nobs ;
			year01               = year01               / Nobs ;
			year02               = year02               / Nobs ;
			year03               = year03               / Nobs ;
			year04               = year04               / Nobs ;
			year05               = year05               / Nobs ;
			lag_sector1          = lag_sector1          / Nobs ;
			lag_sector2          = lag_sector2          / Nobs ;
			lag_sector3          = lag_sector3          / Nobs ;
			lag_sector4          = lag_sector4          / Nobs ;
			lag_sector5          = lag_sector5          / Nobs ;
			lag_sector6          = lag_sector6          / Nobs ;
			lag_sector7          = lag_sector7          / Nobs ;
			gender               = gender               / Nobs ;
			dummy_educ2          = dummy_educ2          / Nobs ;
			dummy_educ3          = dummy_educ3          / Nobs ;
			dummy_educ4          = dummy_educ4          / Nobs ;
			age                  = age                  / Nobs ;
			age_2                = age_2                / Nobs ;
			exper1               = exper1               / Nobs ;
			exper2               = exper2               / Nobs ;
			exper3               = exper3               / Nobs ;
			exper4               = exper4               / Nobs ;
			exper5               = exper5               / Nobs ;
			exper6               = exper6               / Nobs ;
			exper7               = exper7               / Nobs ;
	
		run ;


/*Variance of Variance Residuals*/

		data residuals ;
			set residuals ;
	
			if sector = &s. and idade_new >= 25 and idade_new <= 60 ;
			resid_2 = resid**2 ;

		run ;       


		proc means data = residuals noprint ;
			class sector ;
			var resid_2 ;
			output out = var_residuals&s._lag (drop = _type_ _freq_)
			var = var_resid 
    		n   = nobs ;
		run ;

		data var_residuals&s._lag ;
			retain sector ;
			set var_residuals&s._lag ;

			if sector = . then delete ;

		run ;


/*Variance of Wage Changes - Within Individuals*/

		proc reg data = new_painel noprint OUTSSCP = sscp outest = reg.wage_dif&s. 
		(drop = _MODEL_ _TYPE_ _NAME_ _DEPVAR_ wage_dif) covout ;
		where (ano > 1995 & lag_sector ~= . & sector ~= . & lag_sector = sector & age >= 1) and
		sector = &s. ; 
		model wage_dif = year96 year97 year98 year99 year00 year01 year02 year03 year04 year05 age / noint ;
		OUTPUT OUT = residuals_wagedif RESIDUAL = resid ;
		run ;
		quit ;

		data reg.CPWageDif&s. (keep = Nobs) ;
			set sscp ;
			if _TYPE_ = "N" ;
			Nobs = Intercept ;
		run ;

		data reg.CPWageDif&s. (keep = Nobs year96 year97 year98 year99 year00 year01 year02 year03 year04 year05 age _NAME_) ;
			if _N_ = 1 then set reg.CPWageDif&s. ;
			set sscp ;

			if _NAME_ = "Intercept" or _NAME_ = "wage_dif" or _NAME_ = "" then delete ;

			year96       = year96       / Nobs ;
			year97       = year97       / Nobs ;
			year98       = year98       / Nobs ;
			year99       = year99       / Nobs ;
			year00       = year00       / Nobs ;
			year01       = year01       / Nobs ;
			year02       = year02       / Nobs ;
			year03       = year03       / Nobs ;
			year04       = year04       / Nobs ;
			year05       = year05       / Nobs ;
			age          = age          / Nobs ;
	
		run ;

/*Variance of Variance Wage Changes - Within Individuals*/

		data residuals_wagedif ;
			set residuals_wagedif ;
	
			if (ano > 1995 & lag_sector ~= . & sector ~= . & lag_sector = sector & age >= 1) and sector = &s. ;
			resid_2 = resid**2 ;

		run ;


		proc means data = residuals_wagedif noprint ;
			class sector ;
			var resid_2 ;
			output out = var_residuals_wagedif&s. (drop = _type_ _freq_)
			var = var_resid 
    		n   = nobs ;
		run ;

		data var_residuals_wagedif&s. ;
			retain sector ;
			set var_residuals_wagedif&s. ;

			if sector = . then delete ;

		run ;


/*Dynamic Wage regressions */

	    proc reg data = new_painel noprint OUTSSCP = sscp outest = reg.DynWageSector&s. 
		(drop = _MODEL_ _TYPE_ _NAME_ _DEPVAR_ log_w) covout ;
		where (ano > 1995 & lag_sector ~= . & sector ~= . & lag_sector = sector & age >= 1) and
		sector = &s. ;
		model log_w = lag_log_w year96 year97 year98 year99 year00 year01 year02 year03 year04 year05
    		          gender dummy_educ2 dummy_educ3 dummy_educ4 age age_2 
                      exper1 exper2 exper3 exper4 exper5 exper6 exper7 / noint ;
		run ;
		quit ;


		data reg.CPDynWageSector&s. (keep = Nobs) ;
			set sscp ;
			if _TYPE_ = "N" ;
			Nobs = Intercept ;
		run ;


        data reg.CPDynWageSector&s. (keep = Nobs lag_log_w year96 year97 year98 year99 year00 year01 year02 year03 year04 year05
              					     gender dummy_educ2 dummy_educ3 dummy_educ4 age age_2 
                                     exper1 exper2 exper3 exper4 exper5 exper6 exper7 _NAME_) ;
			if _N_ = 1 then set reg.CPDynWageSector&s. ;
			set sscp ;

			if _NAME_ = "Intercept" or _NAME_ = "log_w" or _NAME_ = "" then delete ;

			lag_log_w            = lag_log_w            / Nobs ;
			year96               = year96               / Nobs ;
			year97               = year97               / Nobs ;
			year98               = year98               / Nobs ;
			year99               = year99               / Nobs ;
			year00               = year00               / Nobs ;
			year01               = year01               / Nobs ;
			year02               = year02               / Nobs ;
			year03               = year03               / Nobs ;
			year04               = year04               / Nobs ;
			year05               = year05               / Nobs ;
			gender               = gender               / Nobs ;
			dummy_educ2          = dummy_educ2          / Nobs ;
			dummy_educ3          = dummy_educ3          / Nobs ;
			dummy_educ4          = dummy_educ4          / Nobs ;
			age                  = age                  / Nobs ;
			age_2                = age_2                / Nobs ;
			exper1               = exper1               / Nobs ;
			exper2               = exper2               / Nobs ;
			exper3               = exper3               / Nobs ;
			exper4               = exper4               / Nobs ;
			exper5               = exper5               / Nobs ;
			exper6               = exper6               / Nobs ;
			exper7               = exper7               / Nobs ;
			
	
		run ;

	%end ;

%mend ;

%Wage_reg ;




data reg.var_residuals ;
	set var_residuals1 var_residuals2 var_residuals3 var_residuals4 
        var_residuals5 var_residuals6 var_residuals7 ;
run ;

data reg.var_residuals_lag ;
	set var_residuals1_lag var_residuals2_lag var_residuals3_lag var_residuals4_lag 
        var_residuals5_lag var_residuals6_lag var_residuals7_lag ;
run ;

data reg.var_residuals_wagedif ;
	set var_residuals_wagedif1 var_residuals_wagedif2 var_residuals_wagedif3 var_residuals_wagedif4 
        var_residuals_wagedif5 var_residuals_wagedif6 var_residuals_wagedif7 ;
run ;

proc delete data = var_residuals1 var_residuals2 var_residuals3 var_residuals4
                   var_residuals5 var_residuals6 var_residuals7
                   var_residuals_wagedif1 var_residuals_wagedif2 var_residuals_wagedif3 var_residuals_wagedif4 
				   var_residuals_wagedif5 var_residuals_wagedif6 var_residuals_wagedif7
                   var_residuals1_lag var_residuals2_lag var_residuals3_lag var_residuals4_lag 
                   var_residuals5_lag var_residuals6_lag var_residuals7_lag ; 
run ;





/*Employment Regressions*/

%macro Emp_reg ;

	%do s = 0 %to 7 ;

		proc reg data = new_painel noprint OUTSSCP = sscp outest = reg.Emp&s.
		(drop = _MODEL_ _TYPE_ _NAME_ _DEPVAR_ sector&s.) covout ;
		where idade_new >= 25 and idade_new <= 60 ;
		model sector&s. = year95 year96 year97 year98 year99 year00 year01 year02 year03 year04 year05
                		  gender dummy_educ2 dummy_educ3 dummy_educ4 age age_2 
                          exper1 exper2 exper3 exper4 exper5 exper6 exper7 / noint ;
		run ;
		quit ;

		data reg.CPEmp&s. (keep = Nobs) ;
			set sscp ;
			if _TYPE_ = "N" ;
			Nobs = Intercept ;
		run ;

		data reg.CPEmp&s. (keep = Nobs year95 year96 year97 year98 year99 year00 year01 year02 year03 year04 year05
    				            gender dummy_educ2 dummy_educ3 dummy_educ4 age age_2 
                                exper1 exper2 exper3 exper4 exper5 exper6 exper7 _NAME_) ;
			if _N_ = 1 then set reg.CPEmp&s. ;
			set sscp ;

			if _NAME_ = "Intercept" or _NAME_ = "sector&s." or _NAME_ = "" then delete ;

			year95       = year95       / Nobs ;
			year96       = year96       / Nobs ;
			year97       = year97       / Nobs ;
			year98       = year98       / Nobs ;
			year99       = year99       / Nobs ;
			year00       = year00       / Nobs ;
			year01       = year01       / Nobs ;
			year02       = year02       / Nobs ;
			year03       = year03       / Nobs ;
			year04       = year04       / Nobs ;
			year05       = year05       / Nobs ;
			gender       = gender       / Nobs ;
			dummy_educ2  = dummy_educ2  / Nobs ;
			dummy_educ3  = dummy_educ3  / Nobs ;
			dummy_educ4  = dummy_educ4  / Nobs ;
			age          = age          / Nobs ;
			age_2        = age_2        / Nobs ;
			exper1       = exper1       / Nobs ;
			exper2       = exper2       / Nobs ;
			exper3       = exper3       / Nobs ;
			exper4       = exper4       / Nobs ;
			exper5       = exper5       / Nobs ;
			exper6       = exper6       / Nobs ;
			exper7       = exper7       / Nobs ;

		run ;

	%end ;

%mend ;	


%Emp_reg ;






/*Transition Regressions*/

%macro Tr_reg ;

	%do s1 = 0 %to 7 ;
		%do s2 = 0 %to 7 ;

			proc reg data = new_painel noprint OUTSSCP = sscp outest = reg.Tr&s1.&s2. 
			(drop = _MODEL_ _TYPE_ _NAME_ _DEPVAR_ sector&s1.&s2.) covout ;
			where idade_new >= 25 and idade_new <= 60 and lag_sector = &s1. ;
			model sector&s1.&s2. = year95 year96 year97 year98 year99 year00 year01 year02 year03 year04 year05
                 			       gender dummy_educ2 dummy_educ3 dummy_educ4 age age_2 
                                   exper1 exper2 exper3 exper4 exper5 exper6 exper7 / noint ;
			run ;
			quit ;

			data reg.CPTr&s1.&s2. (keep = Nobs) ;
				set sscp ;
				if _TYPE_ = "N" ;
				Nobs = Intercept ;
			run ;

			data reg.CPTr&s1.&s2. (keep = Nobs year95 year96 year97 year98 year99 year00 year01 year02 year03 year04 year05
                 			       gender dummy_educ2 dummy_educ3 dummy_educ4 age age_2 
                                   exper1 exper2 exper3 exper4 exper5 exper6 exper7 _NAME_) ;
				if _N_ = 1 then set reg.CPTr&s1.&s2. ;
				set sscp ;
				if _NAME_ = "Intercept" or _NAME_ = "sector&s1.&s2." or _NAME_ = "" then delete ;
	
				year95       = year95       / Nobs ;
				year96       = year96       / Nobs ;
				year97       = year97       / Nobs ;
				year98       = year98       / Nobs ;
				year99       = year99       / Nobs ;
				year00       = year00       / Nobs ;
				year01       = year01       / Nobs ;
				year02       = year02       / Nobs ;
				year03       = year03       / Nobs ;
				year04       = year04       / Nobs ;
				year05       = year05       / Nobs ;
				gender       = gender       / Nobs ;
				dummy_educ2  = dummy_educ2  / Nobs ;
				dummy_educ3  = dummy_educ3  / Nobs ;
				dummy_educ4  = dummy_educ4  / Nobs ;
				age          = age          / Nobs ;
				age_2        = age_2        / Nobs ;
				exper1       = exper1       / Nobs ;
				exper2       = exper2       / Nobs ;
				exper3       = exper3       / Nobs ;
				exper4       = exper4       / Nobs ;
				exper5       = exper5       / Nobs ;
				exper6       = exper6       / Nobs ;
				exper7       = exper7       / Nobs ;

			run ;

		%end ;
	%end ;

%mend ;

%Tr_reg ;



/*%macro Frequencies ;*/
/**/
/*	%do s = 0 %to 7 ;*/
/**/
/*		proc freq data = new_painel noprint ;*/
/*    		where ano = 1995 and idade_new <= 60 - 10 and idade_new >= 25 ;*/
/*				tables freq&s. / out = painel.freq&s. ;*/
/*		run ;*/
/**/
/*		data painel.freq&s. ;*/
/*			set painel.freq&s. ;*/
/*			PERCENT = PERCENT/100.0 ;*/
/*			variance = (PERCENT)*(1 - PERCENT) ;*/
/*        run ;*/
/**/
/*	%end ;*/
/**/
/*	proc freq data = new_painel noprint ;*/
/*    	where ano = 1995 and idade_new <= 60 - 10 and idade_new >= 25 ;*/
/*			tables switches / out = painel.switches ;*/
/*	run ;*/
/**/
/*%mend ;*/
/**/
/*%Frequencies ;*/




%macro Freq_reg ;

	%do s = 0 %to 7 ;

		proc reg data = new_painel noprint OUTSSCP = sscp outest = reg.freq&s.
		(drop = _MODEL_ _TYPE_ _NAME_ _DEPVAR_ freq&s.) covout ;
		where idade_new >= 25 and idade_new <= 50 and ano = 1995 ;
		model freq&s. = lag_sector0 lag_sector1 lag_sector2 lag_sector3 lag_sector4 lag_sector5 lag_sector6 lag_sector7
                        gender dummy_educ2 dummy_educ3 dummy_educ4 age age_2 
                        exper1 exper2 exper3 exper4 exper5 exper6 exper7 / noint ;
		run ;
		quit ;

		data reg.CPfreq&s. (keep = Nobs) ;
			set sscp ;
			if _TYPE_ = "N" ;
			Nobs = Intercept ;
		run ;

		data reg.CPfreq&s. (keep = Nobs lag_sector0 lag_sector1 lag_sector2 lag_sector3 lag_sector4 lag_sector5 lag_sector6 lag_sector7
        gender dummy_educ2 dummy_educ3 dummy_educ4 age age_2 exper1 exper2 exper3 exper4 exper5 exper6 exper7 _NAME_) ;
			if _N_ = 1 then set reg.CPfreq&s. ;
			set sscp ;

			if _NAME_ = "Intercept" or _NAME_ = "freq&s." or _NAME_ = "" then delete ;

			lag_sector0  = lag_sector0  / Nobs ;
			lag_sector1  = lag_sector1  / Nobs ;
			lag_sector2  = lag_sector2  / Nobs ;
			lag_sector3  = lag_sector3  / Nobs ;
			lag_sector4  = lag_sector4  / Nobs ;
			lag_sector5  = lag_sector5  / Nobs ;
			lag_sector6  = lag_sector6  / Nobs ;
			lag_sector7  = lag_sector7  / Nobs ;
			gender       = gender       / Nobs ;
			dummy_educ2  = dummy_educ2  / Nobs ;
			dummy_educ3  = dummy_educ3  / Nobs ;
			dummy_educ4  = dummy_educ4  / Nobs ;
			age          = age          / Nobs ;
			age_2        = age_2        / Nobs ;
			exper1       = exper1       / Nobs ;
			exper2       = exper2       / Nobs ;
			exper3       = exper3       / Nobs ;
			exper4       = exper4       / Nobs ;
			exper5       = exper5       / Nobs ;
			exper6       = exper6       / Nobs ;
			exper7       = exper7       / Nobs ;

		run ;

	%end ;

%mend ;	


%Freq_reg ;






/*Persistence Regressions*/

%macro Pers_reg1998 ;

	%do s = 0 %to 7 ;

		proc reg data = new_painel noprint OUTSSCP = sscp outest = reg.Persistence&s._1998
		(drop = _MODEL_ _TYPE_ _NAME_ _DEPVAR_ sector&s._1998) covout ;
		where idade_new >= 25 and idade_new <= 60-3 and ano = 1995 ;
		model sector&s._1998 = lag_sector0 lag_sector1 lag_sector2 lag_sector3 lag_sector4 lag_sector5 lag_sector6 lag_sector7
        gender dummy_educ2 dummy_educ3 dummy_educ4 age age_2 
        exper1 exper2 exper3 exper4 exper5 exper6 exper7 / noint ;
		run ;
		quit ;

		data reg.CPPersistence&s._1998 (keep = Nobs) ;
			set sscp ;
			if _TYPE_ = "N" ;
			Nobs = Intercept ;
		run ;

		data reg.CPPersistence&s._1998 (keep = Nobs lag_sector0 lag_sector1 lag_sector2 lag_sector3 lag_sector4 lag_sector5 lag_sector6 lag_sector7
        gender dummy_educ2 dummy_educ3 dummy_educ4 age age_2 exper1 exper2 exper3 exper4 exper5 exper6 exper7 _NAME_) ;
			if _N_ = 1 then set reg.CPPersistence&s._1998 ;
			set sscp ;

			if _NAME_ = "Intercept" or _NAME_ = "sector&s._1998" or _NAME_ = "" then delete ;

			lag_sector0  = lag_sector0  / Nobs ;
			lag_sector1  = lag_sector1  / Nobs ;
			lag_sector2  = lag_sector2  / Nobs ;
			lag_sector3  = lag_sector3  / Nobs ;
			lag_sector4  = lag_sector4  / Nobs ;
			lag_sector5  = lag_sector5  / Nobs ;
			lag_sector6  = lag_sector6  / Nobs ;
			lag_sector7  = lag_sector7  / Nobs ;
			gender       = gender       / Nobs ;
			dummy_educ2  = dummy_educ2  / Nobs ;
			dummy_educ3  = dummy_educ3  / Nobs ;
			dummy_educ4  = dummy_educ4  / Nobs ;
			age          = age          / Nobs ;
			age_2        = age_2        / Nobs ;
			exper1       = exper1       / Nobs ;
			exper2       = exper2       / Nobs ;
			exper3       = exper3       / Nobs ;
			exper4       = exper4       / Nobs ;
			exper5       = exper5       / Nobs ;
			exper6       = exper6       / Nobs ;
			exper7       = exper7       / Nobs ;

		run ;

	%end ;

%mend ;	


%Pers_reg1998 ;


%macro Pers_reg2000 ;

	%do s = 0 %to 7 ;

		proc reg data = new_painel noprint OUTSSCP = sscp outest = reg.Persistence&s._2000
		(drop = _MODEL_ _TYPE_ _NAME_ _DEPVAR_ sector&s._2000) covout ;
		where idade_new >= 25 and idade_new <= 60-5 and ano = 1995 ;
		model sector&s._2000 = lag_sector0 lag_sector1 lag_sector2 lag_sector3 lag_sector4 lag_sector5 lag_sector6 lag_sector7 gender 
        dummy_educ2 dummy_educ3 dummy_educ4 age age_2 exper1 exper2 exper3 exper4 exper5 exper6 exper7 / noint ;
		run ;
		quit ;

		data reg.CPPersistence&s._2000 (keep = Nobs) ;
			set sscp ;
			if _TYPE_ = "N" ;
			Nobs = Intercept ;
		run ;

		data reg.CPPersistence&s._2000 (keep = Nobs lag_sector0 lag_sector1 lag_sector2 lag_sector3 lag_sector4 lag_sector5 lag_sector6 lag_sector7
        gender dummy_educ2 dummy_educ3 dummy_educ4 age age_2 exper1 exper2 exper3 exper4 exper5 exper6 exper7 _NAME_) ;
			if _N_ = 1 then set reg.CPPersistence&s._2000 ;
			set sscp ;

			if _NAME_ = "Intercept" or _NAME_ = "sector&s._2000" or _NAME_ = "" then delete ;

			lag_sector0  = lag_sector0  / Nobs ;
			lag_sector1  = lag_sector1  / Nobs ;
			lag_sector2  = lag_sector2  / Nobs ;
			lag_sector3  = lag_sector3  / Nobs ;
			lag_sector4  = lag_sector4  / Nobs ;
			lag_sector5  = lag_sector5  / Nobs ;
			lag_sector6  = lag_sector6  / Nobs ;
			lag_sector7  = lag_sector7  / Nobs ;
			gender       = gender       / Nobs ;
			dummy_educ2  = dummy_educ2  / Nobs ;
			dummy_educ3  = dummy_educ3  / Nobs ;
			dummy_educ4  = dummy_educ4  / Nobs ;
			age          = age          / Nobs ;
			age_2        = age_2        / Nobs ;
			exper1       = exper1       / Nobs ;
			exper2       = exper2       / Nobs ;
			exper3       = exper3       / Nobs ;
			exper4       = exper4       / Nobs ;
			exper5       = exper5       / Nobs ;
			exper6       = exper6       / Nobs ;
			exper7       = exper7       / Nobs ;

		run ;

	%end ;

%mend ;	


%Pers_reg2000 ;




%macro Pers_reg2005 ;

	%do s = 0 %to 7 ;

		proc reg data = new_painel noprint OUTSSCP = sscp outest = reg.Persistence&s._2005
		(drop = _MODEL_ _TYPE_ _NAME_ _DEPVAR_ sector&s._2005) covout ;
		where idade_new >= 25 and idade_new <= 60-10 and ano = 1995 ;
		model sector&s._2005 = lag_sector0 lag_sector1 lag_sector2 lag_sector3 lag_sector4 lag_sector5 lag_sector6 lag_sector7 gender 
        dummy_educ2 dummy_educ3 dummy_educ4 age age_2 exper1 exper2 exper3 exper4 exper5 exper6 exper7 / noint ;
		run ;
		quit ;

		data reg.CPPersistence&s._2005 (keep = Nobs) ;
			set sscp ;
			if _TYPE_ = "N" ;
			Nobs = Intercept ;
		run ;

		data reg.CPPersistence&s._2005 (keep = Nobs lag_sector0 lag_sector1 lag_sector2 lag_sector3 lag_sector4 lag_sector5 lag_sector6 lag_sector7
        gender dummy_educ2 dummy_educ3 dummy_educ4 age age_2 exper1 exper2 exper3 exper4 exper5 exper6 exper7 _NAME_) ;
			if _N_ = 1 then set reg.CPPersistence&s._2005 ;
			set sscp ;

			if _NAME_ = "Intercept" or _NAME_ = "sector&s._2005" or _NAME_ = "" then delete ;

			lag_sector0  = lag_sector0  / Nobs ;
			lag_sector1  = lag_sector1  / Nobs ;
			lag_sector2  = lag_sector2  / Nobs ;
			lag_sector3  = lag_sector3  / Nobs ;
			lag_sector4  = lag_sector4  / Nobs ;
			lag_sector5  = lag_sector5  / Nobs ;
			lag_sector6  = lag_sector6  / Nobs ;
			lag_sector7  = lag_sector7  / Nobs ;
			gender       = gender       / Nobs ;
			dummy_educ2  = dummy_educ2  / Nobs ;
			dummy_educ3  = dummy_educ3  / Nobs ;
			dummy_educ4  = dummy_educ4  / Nobs ;
			age          = age          / Nobs ;
			age_2        = age_2        / Nobs ;
			exper1       = exper1       / Nobs ;
			exper2       = exper2       / Nobs ;
			exper3       = exper3       / Nobs ;
			exper4       = exper4       / Nobs ;
			exper5       = exper5       / Nobs ;
			exper6       = exper6       / Nobs ;
			exper7       = exper7       / Nobs ;

		run ;

	%end ;

%mend ;	


%Pers_reg2005 ;



proc sort data = new_painel ; by lag2_sector lag_sector sector ; run ;

proc freq data = new_painel noprint ;
	where ano >= 1996 ;
	tables lag2_sector * lag_sector * sector / out = painel.Twoyrtransitions ;
run ;



proc delete data = new_painel sscp ; run ;

proc delete data = choices choices_wide ; run ;

proc delete data = residuals residuals_wagedif sscp Frequencies; run ;
