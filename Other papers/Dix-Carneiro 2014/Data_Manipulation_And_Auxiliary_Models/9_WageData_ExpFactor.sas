/*This code follows the individuals drawn for the initial conditions, in order to*/
/*record their wages for the computation of the wage expansion factor used in the*/
/*estimation*/


libname painel "C:\Users\Dix-Carneiro\Documents\Trade and Labor Dynamics\Data for Estimation New\Codes_7_Sectors\PanelRAIS" ;
libname init   "C:\Users\Dix-Carneiro\Documents\Trade and Labor Dynamics\Data for Estimation New\Codes_7_Sectors\PanelRAIS" ;


/*libname painel "C:\Documents and Settings\Rafael Dix Carneiro\My Documents\Thesis\Data for Estimation\Panel of Workers\Data" ;*/
/*libname init   "C:\Documents and Settings\Rafael Dix Carneiro\My Documents\Thesis\Data for Estimation\Panel of Workers\Data" ;*/



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



proc sort data = PIS_list (keep = pis NumberHits) ; by pis ; run ;








/*Now that I have PIS_list, that is, the id's of individuals that were drawn for the initial*/
/*conditions, I merge this list with the panel of workers */

proc sort data = painel.panel9505 ; by pis ; run ;

data init.Wagedata_expfactor ;
	merge PIS_list (in = a) painel.panel9505 ;
	by pis ;
	if a ;
	if idade_new >= 25 and idade_new <= 60 ;
	if ano >= 1995 ;
	if (horas_contr ~= . and horas_contr) > 0 then hourly_wage = rwage_corr / (4.35*horas_contr) ;
	else hourly_wage = . ;
run ;

data init.Wagedata_expfactor ;
  set init.Wagedata_expfactor ;
  do i = 1 to NumberHits ;
    output ;
	replication = i ;
  end ;
  drop i ;
run ;

data init.Wagedata_expfactor ;
	set init.Wagedata_expfactor ;

	if replication = . then replication = 0 ;
	replication = replication + 1 ;
  
run ;

proc sort data = init.Wagedata_expfactor ; by pis ano ; run ;

proc sort data = init.Wagedata_expfactor ; by pis replication ano ; run ;



data id (keep = pis coh NumberHits) ;
	set init.Wagedata_expfactor ;
	by pis ;
	if first.pis ;
	coh = ano - idade_new ;
run ;


data id ;
  set id ;
  do i = 1 to NumberHits ;
    output ;
	replication = i ;
  end ;
  drop i ;
run ;


data id ;
	set id ;

	if replication = . then replication = 0 ;
	replication = replication + 1 ;
  
run ;


proc sort data = id ; by coh pis replication ; run ;


data id (keep = pis replication id coh) ;
	set id ;
	
	obs = _n_ ;
	gen_number = coh - 1935 ;
	id = obs - gen_number*2000 ;
		
run ;


proc sort data = id ; by pis replication id coh ; run ;
proc sort data = init.Wagedata_expfactor ; by pis replication ano ; run ;

data init.wagedata_expfactor (keep = ID coh ano sector hourly_wage educ_new);
	retain ID coh ano sector ;
	merge init.Wagedata_expfactor (in = a) id ;
	by pis replication ;

	if hourly_wage = . then hourly_wage = -999 ;

run ;

proc sort data = init.Wagedata_expfactor ; by coh id ano ; run ;


proc delete data = PIS_LIST ; run ;
