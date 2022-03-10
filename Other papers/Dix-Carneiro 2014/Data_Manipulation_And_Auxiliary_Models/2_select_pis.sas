/*This code reads the am_pis SAS data file with a random sample of all PIS ID's that appear*/
/*In RAIS at least once from 1994 and 2005.*/
/*This code is written in order to track those PIS's over time from 1986 to 2005*/
/*At this point we have a panel of workers in the SAS data file am_rais_final*/
/*However, we still need to convert CNAE and SUBS_IBGE codes to the 7 sectors used in the paper*/
/*This is the second code to be run in order to construct the panel of workers*/



libname rais   "C:\Data\RAIS" ;
libname rafael "C:\Users\Dix-Carneiro\Documents\Trade and Labor Dynamics\Data for Estimation New\Codes_7_Sectors\PanelRAIS" ;



%macro select_pis1 ;

%do ano = 1986 %to 1993 ;

%let ano2 = %substr(&ano.,3,2) ;

data am_pis (rename = pis = pis&ano2.) ;
	set rafael.am_pis ;
run ;

data rafael.rais_aux (keep =  PIS&ano2. 
                              CNPJ&ano2. 
                              SUBS_IBGE&ano2.  
					          CODEMUN&ano2. 
					          GRAU_INSTR&ano2. 
                              SEXO&ano2.
					          REM_DEZ&ano2. 
					          TEMP_EMPR&ano2. 
							  MES_ADM&ano2.
                              MES_DESLIG&ano2. 
							  REM_MEDIA&ano2. ) ;
	set rais.brasil&ano2. ;

	if length(compress(pis&ano2.)) ne 11 then delete ;

run ;


proc sort data = rafael.rais_aux tagsort sortsize = max ; by pis&ano2. ; run ;

proc sort data = am_pis sortsize = max; by pis&ano2. ; run ;

data rafael.am_rais_&ano2. ;
	merge am_pis (in = a) rafael.rais_aux ; 
	by pis&ano2. ; 
	if a ; 
	if compress(MES_DESLIG&ano2.) in ('0','00','13','88') then MES_DESLIG&ano2. = '13' ;
	if compress(MES_ADM&ano2.) in ('0','00','13','88') then MES_ADM&ano2. = '00' ;
	if compress(MES_DESLIG&ano2.) = '13' 
    and (rem_dez&ano2. = 0 or rem_dez&ano2. = .) then rem_dez&ano2. = rem_media&ano2. ;
run ;



/* First we look at those who have jobs in December 31st */
/* Those are mostly workers with nonmissing and nonzero wage for December */
/* For those workers, I select the highest paying job */
/* If there are ties, I select the longest tenure job */

data rais_aux ;
	set rafael.am_rais_&ano2. ;
	if rem_dez&ano2. > 0 and rem_dez&ano2. ~= . ;
run ;

proc sort data = rais_aux tagsort sortsize = max ; 
	by pis&ano2. descending REM_DEZ&ano2. descending TEMP_EMPR&ano2.  ; 
run ;

data rais_aux ;
	set rais_aux ;
	by pis&ano2. ; 
	if first.pis&ano2. ;
run ;




/* Now we look at those who do not have jobs in December 31st */
/* Those are mostly workers with missing or zero wage for December */
/* I first look at whether the job spell lasted for 4 months or more */
/* If it did not, I flag this job with temp = 0 */
/* If it did, I flag this job with temp = 1 */
/* I then select the hihest paying job with temp = 1 */
/* If there are ties, then I select the job with longest tenure */

data rais_aux2 ;
	set rafael.am_rais_&ano2. ;

	if (rem_dez&ano2. <= 0 or rem_dez&ano2. = .) ;

	adm = input(MES_ADM&ano2., 5.) ;
	if adm = 0  then adm = 1 ;
	if adm = 99 then adm = . ;
	deslig = input(MES_DESLIG&ano2., 5.) ;
	if deslig = 13 then deslig = 12 ;
	if deslig = 99 then deslig = . ;

	if deslig > 13 or deslig < 0 then deslig = . ;
	if adm > 13 or adm < 0 then adm = . ;

	temp = deslig - adm + 1 ;

	if temp <= 3 or temp = . then temp = 0 ;
	if temp >= 4 and temp ~= . then temp = 1 ;

run ; 

/* Among the jobs that lasted long enough (4 or more months), select the highest paying one */

proc sort data = rais_aux2 ; by pis&ano2. descending temp descending rem_media&ano2. descending TEMP_EMPR&ano2. ;

data rais_aux2 ;
	set rais_aux2 ;
	by pis&ano2. ;
	if first.pis&ano2. ;
run ;


/* We join those two datasets*/

data rafael.am_rais_&ano2. ;
	set rais_aux rais_aux2 ;
run ;


/* There may still be 2 observations per PIS */
/* We keep the one with positive rem_dez */

proc sort data = rafael.am_rais_&ano2. ; by pis&ano2. descending rem_dez&ano2. ; run ;


data rafael.am_rais_&ano2. ;
	set rafael.am_rais_&ano2. ;
	by pis&ano2. ;
	if first.pis&ano2 ;
	if grau_instr&ano2. < 1 or grau_instr&ano2. > 9 then grau_instr&ano2. = . ;
	if SEXO&ano2. = 3 then SEXO&ano2. = 2 ;
run ;


%end ;

%mend ;

%select_pis1 ;




%macro select_pis2 ;

%do ano = 1994 %to 2005 ;

%let ano2 = %substr(&ano.,3,2) ;

data am_pis (rename = pis = pis&ano2.) ;
	set rafael.am_pis ;
run ;

data rafael.rais_aux (keep =  PIS&ano2. 
                              CNPJ&ano2. 
							  SUBS_IBGE&ano2.
                              CLAS_CNAE&ano2.  
					          CODEMUN&ano2. 
					          GRAU_INSTR&ano2. 
                              SEXO&ano2.
					          REM_DEZ&ano2. 
					          TEMP_EMPR&ano2.
                              HORAS_CONTR&ano2. 
                              IDADE&ano2. 
							  MES_ADM&ano2.
                              MES_DESLIG&ano2. 
							  REM_MEDIA&ano2. ) ;
	set rais.brasil&ano2. ;

	if length(compress(pis&ano2.)) ne 11 then delete ;

run ;


proc sort data = rafael.rais_aux tagsort sortsize = max ; by pis&ano2. ; run ;

proc sort data = am_pis sortsize = max; by pis&ano2. ; run ;

data rafael.am_rais_&ano2. ;
	merge am_pis (in = a) rafael.rais_aux ; 
	by pis&ano2. ; 
	if a ; 
	if compress(MES_DESLIG&ano2.) in ('0','00','13','88') then MES_DESLIG&ano2. = '13' ;
	if compress(MES_ADM&ano2.) in ('0','00','13','88') then MES_ADM&ano2. = '00' ;
	if compress(MES_DESLIG&ano2.) = '13' 
    and (rem_dez&ano2. = 0 or rem_dez&ano2. = .) then rem_dez&ano2. = rem_media&ano2. ;
run ;



/* First we look at those who have jobs in December 31st */
/* Those are mostly workers with nonmissing and nonzero wage for December */
/* For those workers, I select the highest paying job */
/* If there are ties, I select the longest tenure job */

data rais_aux ;
	set rafael.am_rais_&ano2. ;
	if rem_dez&ano2. > 0 and rem_dez&ano2. ~= . ;
run ;

proc sort data = rais_aux tagsort sortsize = max ; 
	by pis&ano2. descending REM_DEZ&ano2. descending horas_contr&ano2. descending TEMP_EMPR&ano2.  ; 
run ;

data rais_aux ;
	set rais_aux ;
	by pis&ano2. ; 
	if first.pis&ano2. ;
run ;




/* Now we look at those who do not have jobs in December 31st */
/* Those are mostly workers with missing or zero wage for December */
/* I first look at whether the job spell lasted for 4 months or more */
/* If it did not, I flag this job with temp = 0 */
/* If it did, I flag this job with temp = 1 */
/* I then select the hihest paying job with temp = 1 */
/* If there are ties, then I select the job with longest tenure */

data rais_aux2 ;
	set rafael.am_rais_&ano2. ;

	if (rem_dez&ano2. <= 0 or rem_dez&ano2. = .) ;

	adm = input(MES_ADM&ano2., 5.) ;
	if adm = 0  then adm = 1 ;
	if adm = 99 then adm = . ;
	deslig = input(MES_DESLIG&ano2., 5.) ;
	if deslig = 13 then deslig = 12 ;
	if deslig = 99 then deslig = . ;

	if deslig > 13 or deslig < 0 then deslig = . ;
	if adm > 13 or adm < 0 then adm = . ;

	temp = deslig - adm + 1 ;

	if temp <= 3 or temp = . then temp = 0 ;
	if temp >= 4 and temp ~= . then temp = 1 ;

run ;   

/* Among the jobs that lasted long enough (4 or more months), select the highest paying one */

proc sort data = rais_aux2 ; by pis&ano2. descending temp descending rem_media&ano2. descending TEMP_EMPR&ano2. ;

data rais_aux2 ;
	set rais_aux2 ;
	by pis&ano2. ;
	if first.pis&ano2. ;
run ;


/* We join those two datasets*/

data rafael.am_rais_&ano2. ;
	set rais_aux rais_aux2 ;
run ;


/* There may still be 2 observations per PIS */
/* We keep the one with positive rem_dez */

proc sort data = rafael.am_rais_&ano2. ; by pis&ano2. descending rem_dez&ano2. ; run ;


data rafael.am_rais_&ano2. ;
	set rafael.am_rais_&ano2. ;
	by pis&ano2. ;
	if first.pis&ano2 ;
	if grau_instr&ano2. < 1 or grau_instr&ano2. > 9 then grau_instr&ano2. = . ;
	if SEXO&ano2. = 3 then SEXO&ano2. = 2 ;
run ;


%end ;

%mend ;

%select_pis2 ;




%macro conserta1 ;

%do ano = 1986 %to 1993 ;

%let ano2 = %substr(&ano.,3,2) ;

data rafael.am_rais_&ano2. ;

	rename PIS&ano2. = PIS ;
	rename CNPJ&ano2. = CNPJ ;
	rename CODEMUN&ano2. = CODEMUN ;
	rename SUBS_IBGE&ano2.  = SUBS_IBGE ;
	rename GRAU_INSTR&ano2. = GRAU_INSTR ;
	rename SEXO&ano2. = SEXO ;
	rename REM_DEZ&ano2. = REM_DEZ ;
	rename TEMP_EMPR&ano2. = TEMP_EMPR ;
	rename MES_ADM&ano2. = MES_ADM ;
    rename MES_DESLIG&ano2. = MES_DESLIG ;
	rename REM_MEDIA&ano2. = REM_MEDIA ;

	set rafael.am_rais_&ano2. ;
	ano = &ano. ;

%end ;

%mend ;

%conserta1 ;




%macro conserta2 ;

%do ano = 1994 %to 2005 ;

%let ano2 = %substr(&ano.,3,2) ;

data rafael.am_rais_&ano2. ;

	rename PIS&ano2. = PIS ;
	rename CNPJ&ano2. = CNPJ ;
	rename CODEMUN&ano2. = CODEMUN ;
	rename SUBS_IBGE&ano2.  = SUBS_IBGE ;
	rename CLAS_CNAE&ano2.  = CLAS_CNAE ;
	rename GRAU_INSTR&ano2. = GRAU_INSTR ;
	rename SEXO&ano2. = SEXO ;
	rename REM_DEZ&ano2. = REM_DEZ ;
	rename TEMP_EMPR&ano2. = TEMP_EMPR ;
	rename HORAS_CONTR&ano2. = HORAS_CONTR ;
	rename IDADE&ano2. = IDADE ;
	rename MES_ADM&ano2. = MES_ADM ;
    rename MES_DESLIG&ano2. = MES_DESLIG ;
	rename REM_MEDIA&ano2. = REM_MEDIA ;

	set rafael.am_rais_&ano2. ;
	ano = &ano. ;

%end ;

%mend ;

%conserta2 ;



data rafael.am_rais_final ;
set rafael.am_rais_86
	rafael.am_rais_87
	rafael.am_rais_88
	rafael.am_rais_89
	rafael.am_rais_90
	rafael.am_rais_91
	rafael.am_rais_92 
 	rafael.am_rais_93
	rafael.am_rais_94
	rafael.am_rais_95
    rafael.am_rais_96 
    rafael.am_rais_97 
    rafael.am_rais_98 
    rafael.am_rais_99 
    rafael.am_rais_00 
    rafael.am_rais_01 
    rafael.am_rais_02 
    rafael.am_rais_03 
    rafael.am_rais_04 
    rafael.am_rais_05 ;
run ;


proc sort data = rafael.am_rais_final sortsize = max; by pis ano ;


proc delete data = rafael.rais_aux 
                   rafael.am_rais_86
				   rafael.am_rais_87
	               rafael.am_rais_88
	               rafael.am_rais_89
	               rafael.am_rais_90
	               rafael.am_rais_91
             	   rafael.am_rais_92 
 	               rafael.am_rais_93
	               rafael.am_rais_94
	               rafael.am_rais_95
                   rafael.am_rais_96 
                   rafael.am_rais_97 
                   rafael.am_rais_98 
                   rafael.am_rais_99 
                   rafael.am_rais_00 
                   rafael.am_rais_01 
                   rafael.am_rais_02 
                   rafael.am_rais_03 
                   rafael.am_rais_04 
                   rafael.am_rais_05 
                   rafael.am_pis 
                   rafael.pis        ; run ;
