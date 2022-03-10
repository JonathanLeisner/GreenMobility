/*This code is written in order to assign a single CNAE classification to each firm.*/
/*The objective of this code is to generate a list of firms and their respective classifications*/
/*so that we can impute and/or correct classifications for firms between 1986 and 1993 when the*/
/*economic ativity of firms were reported at the SUBS_IBGE level*/
/*Firms may have multiple plants with different CNAE classifications. However, this code*/
/*Assigns a single CNAE classification by getting the CNAE classification that is most */
/*frequently reported for each firm between 1994 and 1998. We get one single classification 
/*for each firm. This classification is constant over time*/
/*This is done with the objective of imputing a CNAE classification in case it is missing for some*/
/*observations*/
/*The key data file generated is CNPJ_CNAE9498*/
/*This is the fourth code to be run in order to construct the panel of workers*/





libname rais "C:\Data\RAIS" ;
libname CNPJ "C:\Users\Dix-Carneiro\Documents\Trade and Labor Dynamics\Data for Estimation New\Codes_7_Sectors\PanelRAIS" ;


/* Para cada CNPJ a 8 digitos selecionamos apenas uma classificacao CNAE a 2 digitos */
/* Se um determinado CNPJ8 aparece com mais de uma classificacao CNAE a 2 digitos, */

%macro CNPJ(year) ;

%let ano = %substr(&year,3,2) ;

data rais&ano. (keep = cnpj8 cnae2_&ano. ones) ;

	length cnae2_&ano. $ 3 ;

	set rais.brasil&ano. (keep = cnpj&ano. clas_cnae&ano.) ;

	cnae2_&ano.= substr(CLAS_CNAE&ano.,1,2) ;

	cnpj8 = substr(cnpj&ano.,1,8) ;

	ones = 1 ;

run ;


proc sort data = rais&ano. sortsize = max ;  by cnpj8 cnae2_&ano. ; run;

proc means noprint data = rais&ano. sumsize = max;
	by cnpj8 cnae2_&ano. ;
	var ones ;
	output out = CNPJ.CNPJ_CNAE&ano. (keep = cnpj8 cnae2_&ano. quant&ano.)
	sum(ones) = quant&ano. ;
run ;

proc delete data = rais&ano. ; run ;

proc sort data = CNPJ.CNPJ_CNAE&ano. sortsize = max ;  by cnpj8 descending quant&ano. ; run;

data CNPJ.CNPJ_CNAE&ano. ;
	set CNPJ.CNPJ_CNAE&ano. ;

	by cnpj8 ; 
	if first.cnpj8 ;

run ;

%mend ;

%CNPJ(1994) ;
%CNPJ(1995) ;
%CNPJ(1996) ;
%CNPJ(1997) ;
%CNPJ(1998) ;
/*%CNPJ(1999) ;*/
/*%CNPJ(2000) ;*/
/*%CNPJ(2001) ;*/
/*%CNPJ(2002) ;*/
/*%CNPJ(2003) ;*/
/*%CNPJ(2004) ;*/
/*%CNPJ(2005) ;*/


%macro merging(year1,year2) ;

%let ano1 = %substr(&year1,3,2) ;
%let ano2 = %substr(&year2,3,2) ;

data CNPJ.CNPJ_CNAE&ano2. ;
merge CNPJ.CNPJ_CNAE&ano1. CNPJ.CNPJ_CNAE&ano2. ; 
	by cnpj8 ; 
run ;

proc delete data = CNPJ.CNPJ_CNAE&ano1. ; run ;

proc sort data = CNPJ.CNPJ_CNAE&ano2. sortsize = max ;  by cnpj8 ; run;

%mend ;

%merging(1994,1995) ;
%merging(1995,1996) ;
%merging(1996,1997) ;
%merging(1997,1998) ;
/*%merging(1998,1999) ;*/
/*%merging(1999,2000) ;*/
/*%merging(2000,2001) ;*/
/*%merging(2001,2002) ;*/
/*%merging(2002,2003) ;*/
/*%merging(2003,2004) ;*/
/*%merging(2004,2005) ;*/



data CNPJ.CNPJ_CNAE98_aux ;
set CNPJ.CNPJ_CNAE98 ;
run ;

proc sort data = CNPJ.CNPJ_CNAE98_aux sortsize = max ;  by cnpj8 ; run;

proc transpose data = CNPJ.CNPJ_CNAE98_aux out = CNPJ.CNPJ_CNAE98_aux2 ; 
by cnpj8 ; 
var quant94 quant95 quant96 quant97 quant98 ;
/*quant99 quant00 quant01 quant02 quant03 quant04 quant05*/
run ;

data CNPJ.CNPJ_CNAE98_aux2 ;
set CNPJ.CNPJ_CNAE98_aux2 ;
if _NAME_ = "quant94" then year = 1994 ;
if _NAME_ = "quant95" then year = 1995 ;
if _NAME_ = "quant96" then year = 1996 ;
if _NAME_ = "quant97" then year = 1997 ;
if _NAME_ = "quant98" then year = 1998 ;
/*if _NAME_ = "quant99" then year = 1999 ;*/
/*if _NAME_ = "quant00" then year = 2000 ;*/
/*if _NAME_ = "quant01" then year = 2001 ;*/
/*if _NAME_ = "quant02" then year = 2002 ;*/
/*if _NAME_ = "quant03" then year = 2003 ;*/
/*if _NAME_ = "quant04" then year = 2004 ;*/
/*if _NAME_ = "quant05" then year = 2005 ;*/
drop _NAME_ ;
rename COL1 = quant ;
run ;

proc transpose data = CNPJ.CNPJ_CNAE98_aux out = CNPJ.CNPJ_CNAE98_aux3 ; 
by cnpj8 ; 
var cnae2_94 cnae2_95 cnae2_96 cnae2_97 cnae2_98 ;
/*cnae2_99 cnae2_00 cnae2_01 cnae2_02 cnae2_03 cnae2_04 cnae2_05*/
run ; 

data CNPJ.CNPJ_CNAE98_aux3 ;
set CNPJ.CNPJ_CNAE98_aux3 ;
if _NAME_ = "cnae2_94" then year = 1994 ;
if _NAME_ = "cnae2_95" then year = 1995 ;
if _NAME_ = "cnae2_96" then year = 1996 ;
if _NAME_ = "cnae2_97" then year = 1997 ;
if _NAME_ = "cnae2_98" then year = 1998 ;
drop _NAME_ ;
rename COL1 = CNAE ;
ones = 1 ;
run ;

proc sort data = CNPJ.CNPJ_CNAE98_aux2 sortsize = max ;  by cnpj8 year ; run;

proc sort data = CNPJ.CNPJ_CNAE98_aux3 sortsize = max ;  by cnpj8 year ; run;

data CNPJ.CNPJ_CNAE ;
	merge CNPJ.CNPJ_CNAE98_aux2 CNPJ.CNPJ_CNAE98_aux3 ;
	by cnpj8 year ;
	if cnae = "" or cnae = . then delete ;
run ;

proc sort data = CNPJ.CNPJ_CNAE sortsize = max ;  by cnpj8 year ; run ;

proc sort data = CNPJ.CNPJ_CNAE sortsize = max ;  by cnpj8 cnae ; run;

proc means noprint data = CNPJ.CNPJ_CNAE sumsize = max;
	by cnpj8 cnae ;
	var ones ;
	output out = CNPJ.CNPJ_CNAE (keep = cnpj8 cnae quant)
	sum(ones) = quant ;
run ;

proc sort data = CNPJ.CNPJ_CNAE sortsize = max ;  by cnpj8 descending quant ; run ;

data CNPJ.CNPJ_CNAE9498 (keep = cnpj8 cnae quant) ;
set CNPJ.CNPJ_CNAE ;
by cnpj8 ;
if first.cnpj8 ;
run ;

proc delete data = CNPJ.CNPJ_CNAE98_aux ; run ;

proc delete data = CNPJ.CNPJ_CNAE98_aux2 ; run ;

proc delete data = CNPJ.CNPJ_CNAE98_aux3 ; run ;

proc delete data = CNPJ.cnpj_cnae98 CNPJ.CNPJ_CNAE ; run ;
