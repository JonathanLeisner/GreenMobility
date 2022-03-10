/*Initial conditions are given by sampling, with replacement, 2000 individuals per cohort. */
/*1000 of them with less than a high school degree and 1000 with at least */
/*a high school degree.*/

libname painel "C:\Users\Dix-Carneiro\Documents\Trade and Labor Dynamics\Data for Estimation New\Codes_7_Sectors\PanelRAIS" ;
libname init   "C:\Users\Dix-Carneiro\Documents\Trade and Labor Dynamics\Data for Estimation New\Codes_7_Sectors\PanelRAIS" ;


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


data init.Cohort_Sizes ;
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


data init.Cohort_Sizes ;
	set init.Cohort_Sizes Year_&ano. ;
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

/*proc delete data = Year_95_Age25_educ0 ; run ;*/


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

/*proc delete data = Year_95_Age25_educ1 ; run ;*/



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

/*		proc delete data = Year_95_Age&age._educ&ed. ; run ;*/

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

/*		proc delete data = Year_&ano._educ&ed. ; run ;*/

	%end ;

%end ;

%mend ;

%macro3 ;


proc sort data = PIS_list sortsize = max ; by pis ; run ;
proc sort data = painel.panel9505 sortsize = max ; by pis ; run ;


data init.Initial_Conditions (keep = pis ano_init ano born sexo_new educ_new idade_new sector rwage_corr horas_contr NumberHits) ;

	merge PIS_list (in = a) painel.panel9505 ;

	by pis ;

	if a ;

	if ano = ano_init ;

run ;

proc sort data = init.Initial_Conditions ; by ano_init born ; run ;





/*We keep the whole history of choices up to 9 years*/


proc sort data = PIS_List ; by pis ; run ;
proc sort data = painel.panel9505 ; by pis ; run ;

data Past_Choices (keep = pis ano ano_init sector) ;
	merge PIS_List (in = a) painel.panel9505 ;
	by pis ;
	if a ;
run ;


proc transpose data = Past_Choices out = Past_Choices2 (drop = _NAME_) prefix = sector ; 
	by pis ; 
	id ano ;
	var sector ;
run ;


data Past_Choices2 (keep = pis sector_1-sector_9 exper1 exper2 exper3 exper4 exper5 exper6 exper7) ;
	merge PIS_List (in = a) Past_Choices2 ;
	by pis ;

	if ano_init = 1995 then
	do ;
		sector_1 = sector1994 ;
		sector_2 = sector1993 ;
		sector_3 = sector1992 ;
		sector_4 = sector1991 ;
		sector_5 = sector1990 ;
		sector_6 = sector1989 ;
		sector_7 = sector1988 ;
		sector_8 = sector1987 ;
		sector_9 = sector1986 ;
	end ;

	if ano_init = 1996 then
	do ;
		sector_1 = sector1995 ;
		sector_2 = sector1994 ;
		sector_3 = sector1993 ;
		sector_4 = sector1992 ;
		sector_5 = sector1991 ;
		sector_6 = sector1990 ;
		sector_7 = sector1989 ;
		sector_8 = sector1988 ;
		sector_9 = sector1987 ;
	end ;


	if ano_init = 1997 then
	do ;
		sector_1 = sector1996 ;
		sector_2 = sector1995 ;
		sector_3 = sector1994 ;
		sector_4 = sector1993 ;
		sector_5 = sector1992 ;
		sector_6 = sector1991 ;
		sector_7 = sector1990 ;
		sector_8 = sector1989 ;
		sector_9 = sector1988 ;
	end ;


	if ano_init = 1998 then
	do ;
		sector_1 = sector1997 ;
		sector_2 = sector1996 ;
		sector_3 = sector1995 ;
		sector_4 = sector1994 ;
		sector_5 = sector1993 ;
		sector_6 = sector1992 ;
		sector_7 = sector1991 ;
		sector_8 = sector1990 ;
		sector_9 = sector1989 ;
	end ;

	if ano_init = 1999 then
	do ;
		sector_1 = sector1998 ;
		sector_2 = sector1997 ;
		sector_3 = sector1996 ;
		sector_4 = sector1995 ;
		sector_5 = sector1994 ;
		sector_6 = sector1993 ;
		sector_7 = sector1992 ;
		sector_8 = sector1991 ;
		sector_9 = sector1990 ;
	end ;

	if ano_init = 2000 then
	do ;
		sector_1 = sector1999 ;
		sector_2 = sector1998 ;
		sector_3 = sector1997 ;
		sector_4 = sector1996 ;
		sector_5 = sector1995 ;
		sector_6 = sector1994 ;
		sector_7 = sector1993 ;
		sector_8 = sector1992 ;
		sector_9 = sector1991 ;
	end ;


	if ano_init = 2001 then
	do ;
		sector_1 = sector2000 ;
		sector_2 = sector1999 ;
		sector_3 = sector1998 ;
		sector_4 = sector1997 ;
		sector_5 = sector1996 ;
		sector_6 = sector1995 ;
		sector_7 = sector1994 ;
		sector_8 = sector1993 ;
		sector_9 = sector1992 ;
	end ;


	if ano_init = 2002 then
	do ;
		sector_1 = sector2001 ; 
		sector_2 = sector2000 ;
		sector_3 = sector1999 ;
		sector_4 = sector1998 ;
		sector_5 = sector1997 ;
		sector_6 = sector1996 ;
		sector_7 = sector1995 ;
		sector_8 = sector1994 ;
		sector_9 = sector1993 ;
	end ;


	if ano_init = 2003 then
	do ;
		sector_1 = sector2002 ;
		sector_2 = sector2001 ; 
		sector_3 = sector2000 ;
		sector_4 = sector1999 ;
		sector_5 = sector1998 ;
		sector_6 = sector1997 ;
		sector_7 = sector1996 ;
		sector_8 = sector1995 ;
		sector_9 = sector1994 ;
	end ;


	if ano_init = 2004 then
	do ;
		sector_1 = sector2003 ; 
		sector_2 = sector2002 ;
		sector_3 = sector2001 ; 
		sector_4 = sector2000 ;
		sector_5 = sector1999 ;
		sector_6 = sector1998 ;
		sector_7 = sector1997 ;
		sector_8 = sector1996 ;
		sector_9 = sector1995 ;
	end ;


	if ano_init = 2005 then
	do ;
		sector_1 = sector2004 ;
		sector_2 = sector2003 ; 
		sector_3 = sector2002 ;
		sector_4 = sector2001 ; 
		sector_5 = sector2000 ;
		sector_6 = sector1999 ;
		sector_7 = sector1998 ;
		sector_8 = sector1997 ;
		sector_9 = sector1996 ;
	end ;

	exper1 = 0 ;
	exper2 = 0 ;
	exper3 = 0 ;
	exper4 = 0 ;
	exper5 = 0 ;
	exper6 = 0 ;
	exper7 = 0 ;

	array vars {9} sector_1 sector_2 sector_3 sector_4 sector_5 
                   sector_6 sector_7 sector_8 sector_9 ;
	do i = 1 to 9 ;
		if vars(i) = 1 then exper1 = exper1 + 1 ;
		if vars(i) = 2 then exper2 = exper2 + 1 ;
		if vars(i) = 3 then exper3 = exper3 + 1 ;
		if vars(i) = 4 then exper4 = exper4 + 1 ;
		if vars(i) = 5 then exper5 = exper5 + 1 ;
		if vars(i) = 6 then exper6 = exper6 + 1 ;
		if vars(i) = 7 then exper7 = exper7 + 1 ;
	end ;

run ;


proc sort data = Past_Choices2 ; by pis ; run ;
proc sort data = init.Initial_Conditions ; by pis ; run ;

data init.Initial_Conditions (drop = ano) ;
	merge init.Initial_Conditions (in = a) Past_Choices2 ;
	by pis ;
	if a ;
run ;



data init.Initial_Conditions ;
  set init.Initial_Conditions ;
  do i = 1 to NumberHits ;
    output ;
  end ;
  drop i ;
run ;


proc sort data = init.Initial_Conditions ; by ano_init descending idade_new ; run ;


data init.Initial_Conditions (drop = pis rwage_corr horas_contr exper1 exper2 exper3 exper4 exper5 exper6 exper7 NumberHits) ;
	length ID 8 ;
	length born 8 ;
	set init.Initial_Conditions ;
	ID = _n_ - floor((_n_-0.01)/2000)*2000 ;
	if horas_contr > 0 then hourly_wage = rwage_corr / (4.35*horas_contr) ;
	else hourly_wage = -999 ;
	if hourly_wage = . or hourly_wage <= 0.0 then hourly_wage = -999 ;
run ;



proc delete data = Pis_List ; run ;
proc delete data = Past_Choices ; run ;
proc delete data = Past_Choices2 ; run ;
