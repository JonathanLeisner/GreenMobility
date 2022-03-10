libname painel "C:\Users\Dix-Carneiro\Documents\Trade and Labor Dynamics\Data for Estimation New\Codes_7_Sectors\PanelRAIS" ;
libname boot "C:\Users\Dix-Carneiro\Documents\Trade and Labor Dynamics\Data for Estimation New\Codes_7_Sectors\PanelRAIS\Bootstrap" ;


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






/*Calculate Cohort Sizes*/

data Year_95 (keep = pis born educ) ;
	set painel.panel9505 ;
	if ano = 1995 ;
	if idade_new >= 25 and idade_new <= 60 ;
	born = 1995 - idade_new ;
	educ = (educ_new >= 3) ;
run ;

proc summary data = Year_95 noprint missing ;
	var born educ ;	
   	class born educ ;
   	output out = Year_95 n = total ;
run ;

data Year_95 (keep = born educ total) ;
	set Year_95 ;
	if born ~ = . ;
	if educ  ~= . ;
run ;


data Cohort_Sizes ;
	set Year_95 ;
run ;

proc delete data = Year_95 ; run ;

%macro coh_size ;

%do yr = 1996 %to 2005 ; 

%let ano = %substr(&yr.,3,2) ;

data Year_&ano. (keep = pis born educ) ;
	set painel.panel9505 ;
	if ano = &yr. ;
	if idade_new = 25 ;
	born = &yr. - 25 ;
	educ = (educ_new >= 3) ;
run ;

proc summary data = Year_&ano. noprint missing ;
	var born educ ;	
   	class born educ ;
   	output out = Year_&ano. n = total ;
run ;

data Year_&ano. (keep = born educ total) ;
	set Year_&ano. ;
	if born ~ = . ;
	if educ  ~= . ;
run ;


data Cohort_Sizes ;
	set Cohort_Sizes Year_&ano. ;
run ;

proc delete data = Year_&ano. ; run ;

%end ;

%mend ;

%coh_size ;



/*Get distribution of initial conditions for each cohort/educ pair*/

data Year_95_Age25_educ0 (keep = pis born ano_init) ;
	set painel.panel9505 ;
	if ano = 1995 ;
	ano_init = 1995 ;
	if idade_new = 25 ;
	born = 1995 - 25 ;
	if idade_new = . or sexo_new = . or educ_new = . then delete ;
	educ = (educ_new >= 3) ;
	if educ = 0 ;
run ;

proc surveyselect data = Year_95_Age25_educ0 out = Year_95_Age25_educ0 method = urs sampsize = 1000 seed = 2499 ; run ;

data PIS_list ;
	set Year_95_Age25_educ0 ;
run ;

proc delete data = Year_95_Age25_educ0 ; run ;


data Year_95_Age25_educ1 (keep = pis born ano_init) ;
	set painel.panel9505 ;
	if ano = 1995 ;
	ano_init = 1995 ;
	if idade_new = 25 ;
	born = 1995 - 25 ;
	if idade_new = . or sexo_new = . or educ_new = . then delete ;
	educ = (educ_new >= 3) ;
	if educ = 1 ;
run ;

proc surveyselect data = Year_95_Age25_educ1 out = Year_95_Age25_educ1 method = urs sampsize = 1000 seed = 2499 ; run ;

data PIS_list ;
	set PIS_list Year_95_Age25_educ1 ;
run ;

proc delete data = Year_95_Age25_educ1 ; run ;



%macro macro2 ;

%do age = 26 %to 60 ; 

	%do ed = 0 %to 1 ;

		data Year_95_Age&age._educ&ed. (keep = pis born ano_init) ;
			set painel.panel9505 ;
			if ano = 1995 ;
			ano_init = 1995 ;
			if idade_new = &age. ;
			born = 1995 - &age. ;
			if idade_new = . or sexo_new = . or educ_new = . then delete ;
			educ = (educ_new >= 3) ;
			if educ = &ed. ;
		run ;

		proc surveyselect data = Year_95_Age&age._educ&ed. out = Year_95_Age&age._educ&ed. 
                          method = urs sampsize = 1000 seed = 2499 ; run ;

		data PIS_list ;
			set PIS_list Year_95_Age&age._educ&ed. ;
		run ;

		proc delete data = Year_95_Age&age._educ&ed. ; run ;

	%end ;

%end ;

%mend ;

%macro2 ;




%macro macro3 ;

%do yr = 1996 %to 2005 ; 

	%let ano = %substr(&yr.,3,2) ;

	%do ed = 0 %to 1 ;

		data Year_&ano._educ&ed. (keep = pis born ano_init) ;
			set painel.panel9505 ;
			if idade_new = 25 ;
			if ano = &yr. ;
			ano_init = &yr. ;
			born = &yr. - 25 ;
			if idade_new = . or sexo_new = . or educ_new = . then delete ;
			educ = (educ_new >= 3) ;
			if educ = &ed. ;
		run ;

		proc surveyselect data = Year_&ano._educ&ed. out = Year_&ano._educ&ed.
                          method = urs sampsize = 1000 seed = 2499 ; run ;

		data PIS_list ;
			set PIS_list Year_&ano._educ&ed. ;
		run ;

		proc delete data = Year_&ano._educ&ed. ; run ;

	%end ;

%end ;

%mend ;

%macro3 ;


proc sort data = PIS_list ; by pis ; run ;
proc sort data = new_painel ; by pis ; run ;




data Sample ;
	merge PIS_list (in = a) new_painel ;
	by pis ;
	if a ;
	if idade_new >= 25 and idade_new <= 60 ;
	if ano >= 1995 ;
run ;

proc sort data = Sample ; by pis ano ; run ;


data Cohort_Weights (keep = first_coh);
	set cohort_sizes ;
	if born = 1935 and educ = 0 ;
	first_coh = total ;
run ;


data Cohort_Weights (keep = born educ CohWgt) ;
	if _N_ = 1 then set Cohort_Weights ;
	set cohort_sizes ;
	CohWgt = total / first_coh ;
run ;

proc sort data = Cohort_Weights ; by born educ ; run ;
proc sort data = Sample ; by born educ pis ; run ;

data Sample ;
	merge Sample Cohort_Weights ;
	by born educ ;
run ;


proc sort data = Sample ; by pis ano ; run ;




data Sample ;
  set Sample ;
  do i = 1 to NumberHits ;
    output ;
  end ;
run ;

proc sort data = Sample (drop = NumberHits) ; by pis i ano ; run ;





%macro bootstrap ;

%do s = 1 %to 7 ;

	data WageSector&s. ;
	run ;

%end ;

%do s = 0 %to 7 ;

	data Emp&s. ;
	run ;

%end ;

%do s1 = 0 %to 7 ;
	%do s2 = 0 %to 7 ;
	
		data Tr&s1.&s2. ;
		run ;

	%end ;
%end ;

%do s = 1 %to 7 ;

	data RMSE&s. ;
	run ;

%end ;

%do s = 1 %to 7 ;

	data WageDifRMSE&s. ;
	run ;

%end ;

%do s = 0 %to 7 ;

	data Pers&s._1998 ;
	run ;

%end ;

%do s = 0 %to 7 ;

	data Pers&s._2000 ;
	run ;

%end ;

%do s = 0 %to 7 ;

	data Pers&s._2005 ;
	run ;

%end ;

%do s = 0 %to 7 ;

	data Freq&s. ;
	run ;

%end ;

%do s = 0 %to 7 ;

	data Return&s. ;
	run ;

%end ;

%do i = 1 %to 10 ;

	proc surveyselect data = Sample out = bootsample
		method = urs samprate = 1 ;
	run ;

	data bootsample ;
		set bootsample ;
		Wgt = CohWgt*NumberHits ;
	run ;

	%do s = 1 %to 7 ;


     	proc reg data = bootsample noprint outest = coef_out (drop = _MODEL_ _TYPE_ _DEPVAR_ _RMSE_ log_w) ;
		where sector = &s and idade_new >= 25 and idade_new <= 60 ;
		model log_w = year95 year96 year97 year98 year99 year00 year01 year02 year03 year04 year05
              		  gender dummy_educ2 dummy_educ3 dummy_educ4 age age_2 
                      exper1 exper2 exper3 exper4 exper5 exper6 exper7 / noint ;
		weight Wgt ;
		OUTPUT OUT = residuals RESIDUAL = resid ;
		run ;
		quit ;

		data WageSector&s. ;
			set WageSector&s. coef_out ;
		run ;

		data residuals (keep = Wgt ones idade_new sector resid_2) ;
			set residuals ;
			resid_2 = resid**2 ;
			ones = 1 ;
		run ;

		proc reg data = residuals noprint outest = coef_out (drop = _MODEL_ _TYPE_ _DEPVAR_ _RMSE_ resid_2) ;
		where sector = &s. and idade_new >= 25 and idade_new <= 60 ;
		model resid_2 = ones / noint ;
		weight Wgt ; 
		run ;
		quit ;

		data coef_out (keep = RMSE&s.) ;
			set coef_out ;
			RMSE&s. = sqrt(ones) ;
		run ;

		data RMSE&s. ;
			set RMSE&s. coef_out ;
		run ;

	%end ;

	%do s = 0 %to 7 ;

		proc reg data = bootsample noprint outest = coef_out 
		(drop = _MODEL_ _TYPE_ _DEPVAR_ sector&s.) ;
		where idade_new >= 25 and idade_new <= 60 ;
		model sector&s. = year95 year96 year97 year98 year99 year00 year01 year02 year03 year04 year05
                		  gender dummy_educ2 dummy_educ3 dummy_educ4 age age_2 
						  exper1 exper2 exper3 exper4 exper5 exper6 exper7 / noint ;
		weight Wgt ;
		run ;
		quit ;

		data Emp&s. ;
			set Emp&s. coef_out ;
		run ;

	%end ;

	%do s1 = 0 %to 7 ;
		%do s2 = 0 %to 7 ;

			proc reg data = bootsample noprint outest = coef_out 
			(drop = _MODEL_ _TYPE_  _DEPVAR_ sector&s2.) ;
			where idade_new >= 25 and idade_new <= 60 and lag_sector = &s1 ;
			model sector&s2. = year95 year96 year97 year98 year99 year00 year01 year02 year03 year04 year05
                 			       gender dummy_educ2 dummy_educ3 dummy_educ4 age age_2 
								   exper1 exper2 exper3 exper4 exper5 exper6 exper7 / noint ;
			weight Wgt ;
			run ;
			quit ;

			data Tr&s1.&s2. ;
				set Tr&s1.&s2. coef_out ;
			run ;

		%end ;
	%end ;

	%do s = 1 %to 7 ;

		proc reg data = bootsample noprint outest = coef_out 
			(drop = _MODEL_ _TYPE_  _DEPVAR_ wage_dif) ;
		where (ano > 1995 & lag_sector ~= . & sector ~= . & lag_sector = sector & age >= 1) and
		sector = &s. ; 
		model wage_dif = year96 year97 year98 year99 year00 year01 year02 year03 year04 year05 age / noint ;
		weight Wgt ;
		OUTPUT OUT = residuals RESIDUAL = resid ;
		run ;
		quit ;

		data residuals (keep = Wgt ones idade_new sector resid_2) ;
			set residuals ;
			resid_2 = resid**2 ;
			ones = 1 ;
		run ;

		proc reg data = residuals noprint outest = coef_out (drop = _MODEL_ _TYPE_ _DEPVAR_ _RMSE_ resid_2) ;
		where sector = &s. and idade_new >= 25 and idade_new <= 60 ;
		model resid_2 = ones / noint ;
		weight Wgt ; 
		run ;
		quit ;

		data coef_out (keep = WageDifRMSE&s.) ;
			set coef_out ;
			WageDifRMSE&s. = sqrt(ones) ;
		run ;

		data WageDifRMSE&s. ;
			set WageDifRMSE&s. coef_out ;
		run ;		
	
	%end ;


	%do s = 0 %to 7 ;

		proc reg data = bootsample noprint outest = coef_out
		(drop = _MODEL_ _TYPE_ _DEPVAR_ sector&s._1998) ;
		where idade_new >= 25 and idade_new <= 60-3 and ano = 1995 ;
		model sector&s._1998 = lag_sector0 lag_sector1 lag_sector2 lag_sector3 lag_sector4 lag_sector5 lag_sector6 lag_sector7
        gender dummy_educ2 dummy_educ3 dummy_educ4 age age_2 exper1 exper2 exper3 exper4 exper5 exper6 exper7 / noint ;
		weight Wgt ;
		run ;
		quit ;

		data Pers&s._1998 ;
			set Pers&s._1998 coef_out ;
		run ;

	%end ;


	%do s = 0 %to 7 ;

		proc reg data = bootsample noprint outest = coef_out
		(drop = _MODEL_ _TYPE_ _DEPVAR_ sector&s._2000) ;
		where idade_new >= 25 and idade_new <= 60-5 and ano = 1995 ;
		model sector&s._2000 = lag_sector0 lag_sector1 lag_sector2 lag_sector3 lag_sector4 lag_sector5 lag_sector6 lag_sector7
        gender dummy_educ2 dummy_educ3 dummy_educ4 age age_2 exper1 exper2 exper3 exper4 exper5 exper6 exper7 / noint ;
		weight Wgt ;
		run ;
		quit ;

		data Pers&s._2000 ;
			set Pers&s._2000 coef_out ;
		run ;

	%end ;


	%do s = 0 %to 7 ;

		proc reg data = bootsample noprint outest = coef_out
		(drop = _MODEL_ _TYPE_ _DEPVAR_ sector&s._2005) ;
		where idade_new >= 25 and idade_new <= 60-10 and ano = 1995 ;
		model sector&s._2005 = lag_sector0 lag_sector1 lag_sector2 lag_sector3 lag_sector4 lag_sector5 lag_sector6 lag_sector7
        gender dummy_educ2 dummy_educ3 dummy_educ4 age age_2 exper1 exper2 exper3 exper4 exper5 exper6 exper7 / noint ;
		weight Wgt ;
		run ;
		quit ;

		data Pers&s._2005 ;
			set Pers&s._2005 coef_out ;
		run ;

	%end ;


	%do s = 0 %to 7 ;

		proc reg data = bootsample noprint outest = coef_out
		(drop = _MODEL_ _TYPE_ _DEPVAR_ freq&s.) ;
		where idade_new >= 25 and idade_new <= 50 and ano = 1995 ;
		model freq&s. = lag_sector0 lag_sector1 lag_sector2 lag_sector3 lag_sector4 lag_sector5 lag_sector6 lag_sector7
        gender dummy_educ2 dummy_educ3 dummy_educ4 age age_2 exper1 exper2 exper3 exper4 exper5 exper6 exper7 / noint ;
		weight Wgt ;
		run ;
		quit ;

		data Freq&s. ;
			set Freq&s. coef_out ;
		run ;

	%end ;


	%do s = 0 %to 7 ;

		proc reg data = bootsample noprint outest = coef_out 
		(drop = _MODEL_ _TYPE_ _DEPVAR_ sector&s.) ;
		where lag2_sector = &s. and lag_sector ~= &s. and lag_sector ~= . and idade_new >= 25 and idade_new <= 60 ;
		model sector&s. = year95 year96 year97 year98 year99 year00 year01 year02 year03 year04 year05
                		  gender dummy_educ2 dummy_educ3 dummy_educ4 age age_2 
						  exper1 exper2 exper3 exper4 exper5 exper6 exper7 / noint ;
		weight Wgt ;
		run ;
		quit ;

		data Return&s. ;
			set Return&s. coef_out ;
		run ;

	%end ;

%end ;

%do s = 1 %to 7 ;

	data WageSector&s. ;
		set WageSector&s. ;
		if year95 ~= . ;
	run ;

%end ;


%do s = 0 %to 7 ;

	data Emp&s. ;
		set Emp&s. ;
		if year95 ~= . ;
	run ;

%end ;

%do s1 = 0 %to 7 ;
	%do s2 = 0 %to 7 ;
	
		data Tr&s1.&s2. ;
			set Tr&s1.&s2. ;
			if year95 ~= . ;
		run ;

	%end ;
%end ;

%do s = 1 %to 7 ;

	data RMSE&s. ;
		set RMSE&s. ;
		if RMSE&s. ~= . ;
	run ;

%end ;

%do s = 1 %to 7 ;

	data WageDifRMSE&s. ;
		set WageDifRMSE&s. ;
		if WageDifRMSE&s. ~= . ;
	run ;

%end ;




%do s = 0 %to 7 ;

	data Pers&s._1998 ;
		set Pers&s._1998 ;
		if lag_sector0 ~= . ;
	run ;

%end ;


%do s = 0 %to 7 ;

	data Pers&s._2000 ;
		set Pers&s._2000 ;
		if lag_sector0 ~= . ;
	run ;

%end ;


%do s = 0 %to 7 ;

	data Pers&s._2005 ;
		set Pers&s._2005 ;
		if lag_sector0 ~= . ;
	run ;

%end ;


%do s = 0 %to 7 ;

	data Freq&s. ;
		set Freq&s. ;
		if lag_sector0 ~= . ;
	run ;

%end ;


%do s = 0 %to 7 ;

	data Return&s. ;
		set Return&s. ;
		if year95 ~= . ;
	run ;

%end ;

%mend ; 

%bootstrap ;





proc delete data = bootsample coef_out New_painel Residuals ; run ;




%macro rename ;

%do s = 1 %to 7 ;

	data Wagesector&s. ;
		retain obs ;
		set Wagesector&s. ;

		obs = _n_ ;

		rename _RMSE_      = rmse_&s.    ;
    	rename year95      = beta_&s._1  ;
		rename year96      = beta_&s._2  ;
		rename year97      = beta_&s._3  ;
		rename year98      = beta_&s._4  ;
		rename year99      = beta_&s._5  ;
		rename year00      = beta_&s._6  ;
		rename year01      = beta_&s._7  ;
		rename year02      = beta_&s._8  ;
		rename year03      = beta_&s._9  ;
		rename year04      = beta_&s._10 ;
		rename year05      = beta_&s._11 ;
		rename gender      = beta_&s._12 ;
		rename dummy_educ2 = beta_&s._13 ;
		rename dummy_educ3 = beta_&s._14 ;
		rename dummy_educ4 = beta_&s._15 ;
		rename age         = beta_&s._16 ;
		rename age_2       = beta_&s._17 ;
		rename exper1      = beta_&s._18 ;
		rename exper2      = beta_&s._19 ;
		rename exper3      = beta_&s._20 ;
		rename exper4      = beta_&s._21 ;
		rename exper5      = beta_&s._22 ;
		rename exper6      = beta_&s._23 ;
		rename exper7      = beta_&s._24 ;

	run ;

%end ;



%do s = 1 %to 7 ;

	data Wagesector&s. (drop = rmse_&s.) ;
		retain obs ;
		set Wagesector&s. ;

	run ;

	proc sort data = Wagesector&s. ; by obs ; run ;

%end ;


%do s = 0 %to 7 ;

	data Emp&s. (drop = _RMSE_) ; 
		retain obs ;
		set Emp&s. ;

		obs = _n_ ;

		rename year95      = gamma_&s._1  ;
		rename year96      = gamma_&s._2  ;
		rename year97      = gamma_&s._3  ;
		rename year98      = gamma_&s._4  ;
		rename year99      = gamma_&s._5  ;
		rename year00      = gamma_&s._6  ;
		rename year01      = gamma_&s._7  ;
		rename year02      = gamma_&s._8  ;
		rename year03      = gamma_&s._9  ;
		rename year04      = gamma_&s._10 ;
		rename year05      = gamma_&s._11 ;
		rename gender      = gamma_&s._12 ;
		rename dummy_educ2 = gamma_&s._13 ;
		rename dummy_educ3 = gamma_&s._14 ;
		rename dummy_educ4 = gamma_&s._15 ;
		rename age         = gamma_&s._16 ;
		rename age_2       = gamma_&s._17 ;
		rename exper1      = gamma_&s._18 ;
		rename exper2      = gamma_&s._19 ;
		rename exper3      = gamma_&s._20 ;
		rename exper4      = gamma_&s._21 ;
		rename exper5      = gamma_&s._22 ;
		rename exper6      = gamma_&s._23 ;
		rename exper7      = gamma_&s._24 ;

	run ;

	proc sort data = Emp&s. ; by obs ; run ;

%end ;


%do s1 = 0 %to 7 ;
	%do s2 = 0 %to 7 ;

		data Tr&s1.&s2. (drop = _RMSE_) ;
			retain obs ;
			set Tr&s1.&s2. ;

			obs = _n_ ;

			rename year95      = phi_&s1.&s2._1  ;
			rename year96      = phi_&s1.&s2._2  ;
			rename year97      = phi_&s1.&s2._3  ;
			rename year98      = phi_&s1.&s2._4  ;
			rename year99      = phi_&s1.&s2._5  ;
			rename year00      = phi_&s1.&s2._6  ;
			rename year01      = phi_&s1.&s2._7  ;
			rename year02      = phi_&s1.&s2._8  ;
			rename year03      = phi_&s1.&s2._9  ;
			rename year04      = phi_&s1.&s2._10 ;
			rename year05      = phi_&s1.&s2._11 ;
			rename gender      = phi_&s1.&s2._12 ;
			rename dummy_educ2 = phi_&s1.&s2._13 ;
			rename dummy_educ3 = phi_&s1.&s2._14 ;
			rename dummy_educ4 = phi_&s1.&s2._15 ;
			rename age         = phi_&s1.&s2._16 ;
			rename age_2       = phi_&s1.&s2._17 ;
			rename exper1      = phi_&s1.&s2._18 ;
			rename exper2      = phi_&s1.&s2._19 ;
			rename exper3      = phi_&s1.&s2._20 ;
			rename exper4      = phi_&s1.&s2._21 ;
			rename exper5      = phi_&s1.&s2._22 ;
			rename exper6      = phi_&s1.&s2._23 ;
			rename exper7      = phi_&s1.&s2._24 ;

		run ;

		proc sort data = Tr&s1.&s2. ; by obs ; run ;

	%end ;

%end ; 	


%do s = 1 %to 7 ;

	data RMSE&s. (keep = obs RMSE&s.) ;
		retain obs ;
		set RMSE&s. ;
		obs = _n_ ;
	run ;

	proc sort data = RMSE&s. ; by obs ; run ;

%end ;


%do s = 1 %to 7 ;

	data WageDifRMSE&s. (keep = obs WageDifRMSE&s.) ;
		retain obs ;
		set WageDifRMSE&s. ;
		obs = _n_ ;
	run ;

	proc sort data = WageDifRMSE&s. ; by obs ; run ;

%end ;



%do s = 0 %to 7 ;

	data Pers&s._1998 (drop = _RMSE_) ;
		retain obs ;
		set Pers&s._1998 ;

		obs = _n_ ;

		rename lag_sector0  = xsi1998_&s._1 ;
		rename lag_sector1  = xsi1998_&s._2 ;
		rename lag_sector2  = xsi1998_&s._3 ;
		rename lag_sector3  = xsi1998_&s._4 ;
		rename lag_sector4  = xsi1998_&s._5 ;
		rename lag_sector5  = xsi1998_&s._6 ;
		rename lag_sector6  = xsi1998_&s._7 ;
		rename lag_sector7  = xsi1998_&s._8 ;
		rename gender       = xsi1998_&s._9 ;
		rename dummy_educ2  = xsi1998_&s._10 ;
		rename dummy_educ3  = xsi1998_&s._11 ;
		rename dummy_educ4  = xsi1998_&s._12 ;
		rename age          = xsi1998_&s._13 ;
		rename age_2        = xsi1998_&s._14 ;
		rename exper1       = xsi1998_&s._15 ;
		rename exper2       = xsi1998_&s._16 ;
		rename exper3       = xsi1998_&s._17 ;
		rename exper4       = xsi1998_&s._18 ;
		rename exper5       = xsi1998_&s._19 ;
		rename exper6       = xsi1998_&s._20 ;
		rename exper7       = xsi1998_&s._21 ;

	run ;

	proc sort data = Pers&s._1998 ; by obs ; run ;

%end ;



%do s = 0 %to 7 ;

	data Pers&s._2000 (drop = _RMSE_) ;
		retain obs ;
		set Pers&s._2000 ;

		obs = _n_ ;

		rename lag_sector0  = xsi2000_&s._1 ;
		rename lag_sector1  = xsi2000_&s._2 ;
		rename lag_sector2  = xsi2000_&s._3 ;
		rename lag_sector3  = xsi2000_&s._4 ;
		rename lag_sector4  = xsi2000_&s._5 ;
		rename lag_sector5  = xsi2000_&s._6 ;
		rename lag_sector6  = xsi2000_&s._7 ;
		rename lag_sector7  = xsi2000_&s._8 ;
		rename gender       = xsi2000_&s._9 ;
		rename dummy_educ2  = xsi2000_&s._10 ;
		rename dummy_educ3  = xsi2000_&s._11 ;
		rename dummy_educ4  = xsi2000_&s._12 ;
		rename age          = xsi2000_&s._13 ;
		rename age_2        = xsi2000_&s._14 ;
		rename exper1       = xsi2000_&s._15 ;
		rename exper2       = xsi2000_&s._16 ;
		rename exper3       = xsi2000_&s._17 ;
		rename exper4       = xsi2000_&s._18 ;
		rename exper5       = xsi2000_&s._19 ;
		rename exper6       = xsi2000_&s._20 ;
		rename exper7       = xsi2000_&s._21 ;

	run ;

	proc sort data = Pers&s._2000 ; by obs ; run ;

%end ;



%do s = 0 %to 7 ;

	data Pers&s._2005 (drop = _RMSE_) ;
		retain obs ;
		set Pers&s._2005 ;

		obs = _n_ ;

		rename lag_sector0  = xsi2005_&s._1 ;
		rename lag_sector1  = xsi2005_&s._2 ;
		rename lag_sector2  = xsi2005_&s._3 ;
		rename lag_sector3  = xsi2005_&s._4 ;
		rename lag_sector4  = xsi2005_&s._5 ;
		rename lag_sector5  = xsi2005_&s._6 ;
		rename lag_sector6  = xsi2005_&s._7 ;
		rename lag_sector7  = xsi2005_&s._8 ;
		rename gender       = xsi2005_&s._9 ;
		rename dummy_educ2  = xsi2005_&s._10 ;
		rename dummy_educ3  = xsi2005_&s._11 ;
		rename dummy_educ4  = xsi2005_&s._12 ;
		rename age          = xsi2005_&s._13 ;
		rename age_2        = xsi2005_&s._14 ;
		rename exper1       = xsi2005_&s._15 ;
		rename exper2       = xsi2005_&s._16 ;
		rename exper3       = xsi2005_&s._17 ;
		rename exper4       = xsi2005_&s._18 ;
		rename exper5       = xsi2005_&s._19 ;
		rename exper6       = xsi2005_&s._20 ;
		rename exper7       = xsi2005_&s._21 ;

	run ;

	proc sort data = Pers&s._2005 ; by obs ; run ;

%end ;



%do s = 0 %to 7 ;

	data Freq&s. (drop = _RMSE_) ;
		retain obs ;
		set Freq&s. ;

		obs = _n_ ;

		rename lag_sector0  = eta_&s._1 ;
		rename lag_sector1  = eta_&s._2 ;
		rename lag_sector2  = eta_&s._3 ;
		rename lag_sector3  = eta_&s._4 ;
		rename lag_sector4  = eta_&s._5 ;
		rename lag_sector5  = eta_&s._6 ;
		rename lag_sector6  = eta_&s._7 ;
		rename lag_sector7  = eta_&s._8 ;
		rename gender       = eta_&s._9 ;
		rename dummy_educ2  = eta_&s._10 ;
		rename dummy_educ3  = eta_&s._11 ;
		rename dummy_educ4  = eta_&s._12 ;
		rename age          = eta_&s._13 ;
		rename age_2        = eta_&s._14 ;
		rename exper1       = eta_&s._15 ;
		rename exper2       = eta_&s._16 ;
		rename exper3       = eta_&s._17 ;
		rename exper4       = eta_&s._18 ;
		rename exper5       = eta_&s._19 ;
		rename exper6       = eta_&s._20 ;
		rename exper7       = eta_&s._21 ;

	run ;

	proc sort data = Freq&s. ; by obs ; run ;

%end ;


%do s = 0 %to 7 ;

	data Return&s. (drop = _RMSE_) ; 
		retain obs ;
		set Return&s. ;

		obs = _n_ ;

		rename year95      = rho_&s._1  ;
		rename year96      = rho_&s._2  ;
		rename year97      = rho_&s._3  ;
		rename year98      = rho_&s._4  ;
		rename year99      = rho_&s._5  ;
		rename year00      = rho_&s._6  ;
		rename year01      = rho_&s._7  ;
		rename year02      = rho_&s._8  ;
		rename year03      = rho_&s._9  ;
		rename year04      = rho_&s._10 ;
		rename year05      = rho_&s._11 ;
		rename gender      = rho_&s._12 ;
		rename dummy_educ2 = rho_&s._13 ;
		rename dummy_educ3 = rho_&s._14 ;
		rename dummy_educ4 = rho_&s._15 ;
		rename age         = rho_&s._16 ;
		rename age_2       = rho_&s._17 ;
		rename exper1      = rho_&s._18 ;
		rename exper2      = rho_&s._19 ;
		rename exper3      = rho_&s._20 ;
		rename exper4      = rho_&s._21 ;
		rename exper5      = rho_&s._22 ;
		rename exper6      = rho_&s._23 ;
		rename exper7      = rho_&s._24 ;

	run ;

	proc sort data = Return&s. ; by obs ; run ;

%end ;




data boot.Bootstrap ;
	set Wagesector1 ;
run ;

proc delete data = Wagesector1 ; run ;

%do s = 2 %to 7 ;

	data boot.Bootstrap ;
		merge boot.Bootstrap Wagesector&s. ;
		by obs ;
	run ;

	proc delete data = Wagesector&s. ; run ;

%end ;


%do s = 0 %to 7 ;

	data boot.Bootstrap ;
		merge boot.Bootstrap Emp&s. ;
		by obs ;
	run ;

	proc delete data = Emp&s. ; run ;

%end ;


%do s1 = 0 %to 7 ;
	%do s2 = 0 %to 7 ;

		data boot.Bootstrap ;
			merge boot.Bootstrap Tr&s1.&s2. ;
		by obs ;
		run ;

		proc delete data = Tr&s1.&s2. ; run ;

	%end ;
%end ;


%do s = 1 %to 7 ;

	data boot.Bootstrap ;
		merge boot.Bootstrap RMSE&s. ;
		by obs ;
	run ;

	proc delete data = RMSE&s. ; run ;

%end ;


%do s = 1 %to 7 ;

	data boot.Bootstrap ;
		merge boot.Bootstrap WageDifRMSE&s. ;
		by obs ;
	run ;

	proc delete data = WageDifRMSE&s. ; run ;

%end ;


%do s = 0 %to 7 ;

	data boot.Bootstrap ;
		merge boot.Bootstrap Pers&s._1998 ;
		by obs ;
	run ;

	proc delete data = Pers&s._1998 ; run ;

%end ;


%do s = 0 %to 7 ;

	data boot.Bootstrap ;
		merge boot.Bootstrap Pers&s._2000 ;
		by obs ;
	run ;

	proc delete data = Pers&s._2000 ; run ;

%end ;


%do s = 0 %to 7 ;

	data boot.Bootstrap ;
		merge boot.Bootstrap Pers&s._2005 ;
		by obs ;
	run ;

	proc delete data = Pers&s._2005 ; run ;

%end ;


%do s = 0 %to 7 ;

	data boot.Bootstrap ;
		merge boot.Bootstrap Freq&s. ;
		by obs ;
	run ;

	proc delete data = Freq&s. ; run ;

%end ;


%do s = 0 %to 7 ;

	data boot.Bootstrap ;
		merge boot.Bootstrap Return&s. ;
		by obs ;
	run ;

	proc delete data = Return&s. ; run ;

%end ;



%mend rename ;

%rename ;



/*proc delete data = Sample bootsample coef_out Cohort_weights New_Painel PIS_List Residuals Cohort_sizes ; run ;*/
/**/
/**/
/**/
/*PROC CORR DATA = boot.Bootstrap noprint COV NOCORR OUTP = boot.BOOT_COVAR_SIM ;*/
/*	VAR _ALL_ ;*/
/*RUN ;*/
/**/
/*data boot.BOOT_COVAR_SIM (drop = _TYPE_ obs _NAME_) ;*/
/*	set boot.BOOT_COVAR_SIM ;*/
/*	if _NAME_ = "obs" or _NAME_ = "" then delete ;*/
/*run ;*/
