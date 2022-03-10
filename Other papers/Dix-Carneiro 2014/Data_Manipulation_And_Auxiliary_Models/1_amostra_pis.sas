/*This code lists all the PIS ID's that appear at RAIS at least once between 1995 and 2005*/
/*In the end, it selects a 1% random sample of those ID's*/
/*The file am_pis contains this random sample of ID's*/
/*It is the first code to be run in order to obtain the panel of workers*/



libname rais   "C:\Data\RAIS" ;
libname rafael "C:\Users\Dix-Carneiro\Documents\Trade and Labor Dynamics\Data for Estimation New\Codes_7_Sectors\PanelRAIS" ;


%macro pis ;

%do ano = 1995 %to 2005 ;

%let ano2 = %substr(&ano,3,2) ;

data pis&ano2. (keep = pis&ano2. rename = pis&ano2. = pis) ;
	set rais.brasil&ano2. ;
	if rem_dez&ano2. = 0 or rem_dez&ano2. = . then delete ;
	if horas_contr&ano2. = 0 or horas_contr&ano2. = . then delete ;
	if length(compress(pis&ano2.)) ne 11 then delete ;
run ;

proc sort data = pis&ano2. sortsize = max nodup ;  by pis ; run;

%end ;

%mend ; 

%pis ;


data pis ;
merge pis95 pis96 pis97 pis98 pis99 pis00 
      pis01 pis02 pis03 pis04 pis05 ;
by pis ;
run ;


proc sort data = pis sortsize = max nodup ;  by pis ; run;


%macro delete_files ;

%do ano = 1995 %to 2005 ;

%let ano2 = %substr(&ano,3,2) ;

proc delete data = pis&ano2. ; run;

%end ;

%mend ; 

%delete_files ;



data teste_pis2 (compress=yes) ;
	set pis ;

		/*1° passo - qualquer PIS com tamanho diferente de 11 é eliminado*/	
		if length(compress(pis)) ne 11 then delete ; 

/*		programa para verificação da consistência do CPF - fórmula encotrada na internet e muito consistente*/
	    dig1=input(compress(substr(compress(pis),1,1)),1.) ;
		dig2=input(compress(substr(compress(pis),2,1)),1.) ;
		dig3=input(compress(substr(compress(pis),3,1)),1.) ;
		dig4=input(compress(substr(compress(pis),4,1)),1.) ;
		dig5=input(compress(substr(compress(pis),5,1)),1.) ;
		dig6=input(compress(substr(compress(pis),6,1)),1.) ;
		dig7=input(compress(substr(compress(pis),7,1)),1.) ;
		dig8=input(compress(substr(compress(pis),8,1)),1.) ;
		dig9=input(compress(substr(compress(pis),9,1)),1.) ;
		dig10=input(compress(substr(compress(pis),10,1)),1.) ;
		dig_pis=input(compress(substr(compress(pis),11,1)),1.) ;
		

		digito1=(11-mod(sum(dig1*3,dig2*2,dig3*9,dig4*8,dig5*7,dig6*6,dig7*5,dig8*4,dig9*3,dig10*2),11));
			if digito1=10 or digito1=11 then digito1=0;
				
		if dig_pis ne digito1 then pis_invalido=1;
			else pis_invalido=0;

/*		drop dig1 dig2 dig3 dig4 dig5 dig6 dig7 dig8 dig9 dig10 dig_pis94 ;*/
run;

proc freq data = teste_pis2 ;
	table pis_invalido / out = rafael.pis_invalido ;
run;


data rafael.pis (keep = pis) ;
	set teste_pis2 ;
	if pis_invalido = 1 then delete ;
run ;






%macro MTCNTOBS ;

/*Get the number of observations in rafael.pis and draws a 1% random sample of pis*/

%let DSID = %sysfunc(open(rafael.pis, IS));
%let anobs = %sysfunc(attrn(&DSID, ANOBS));
%let whstmt = %sysfunc(attrn(&DSID, WHSTMT));
%if &anobs = 1 & &whstmt = 0 %then
%do;
%let counted = %sysfunc(attrn(&DSID, NLOBS));
%end;

%let samplesize = %SYSEVALF(0.01*(&counted.),integer) ;

proc surveyselect data = rafael.pis out = rafael.am_pis sampsize = &samplesize seed = 2499; run;

%mend MTCNTOBS;

%MTCNTOBS;

