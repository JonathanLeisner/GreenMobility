/*This code is written in order to correct for inconsistencies in the report of*/
/*age, gender and education. We also impute some wages in case they are missing.*/
/*This is the sixth and last code to be run in order to generate the panel.*/



libname painel "C:\Users\Dix-Carneiro\Documents\Trade and Labor Dynamics\Data for Estimation New\Codes_7_Sectors\PanelRAIS" ;



proc sort data = painel.Painel9505 sortsize = max ; by pis ano ; run ;


data painel_aux_lag (keep = pis ano lag_cnpj8 lag_wage lag_mes_deslig lag_setor) ;
	set painel.Painel9505 ;
	ano = ano + 1 ;
	lag_wage = wage ;
	lag_cnpj8 = CNPJ8 ;
	lag_mes_deslig = MES_DESLIG ;
	lag_setor = setor ;
run ;


data painel_aux_lag2 (keep = pis ano lag2_cnpj8 lag2_wage lag2_mes_deslig lag2_setor) ;
	set painel.Painel9505 ;
	ano = ano + 2 ;
	lag2_wage = wage ;
	lag2_cnpj8 = CNPJ8 ;
	lag2_mes_deslig = MES_DESLIG ;
	lag2_setor = setor ;
run ;



data painel_aux_lead (keep = pis ano lead_cnpj8 lead_wage lead_mes_adm) ;
	set painel.Painel9505 ;
	ano = ano - 1 ;
	lead_wage = wage ;
	lead_cnpj8 = CNPJ8 ;
	lead_mes_adm = MES_ADM ;
run ;



data painel_aux_lead2 (keep = pis ano lead2_cnpj8 lead2_wage lead2_mes_adm) ;
	set painel.Painel9505 ;
	ano = ano - 2 ;
	lead2_wage = wage ;
	lead2_cnpj8 = CNPJ8 ;
	lead2_mes_adm = MES_ADM ;
run ;



data painel.painel_final ;
	merge painel.Painel9505 (in = a) painel_aux_lag ;
	by pis ano ;
	if a ;
run ;


data painel.painel_final ;
	merge painel.painel_final (in = a) painel_aux_lag2 ;
	by pis ano ;
	if a ;
run ;


data painel.painel_final ;
	merge painel.painel_final (in = a) painel_aux_lead ;
	by pis ano ;
	if a ;
run ;



data painel.painel_final ;
	merge painel.painel_final (in = a) painel_aux_lead2 ;
	by pis ano ;
	if a ;
run ;

proc delete data = painel_aux_lag 
                   painel_aux_lag2 
                   painel_aux_lead 
                   painel_aux_lead2 
                   painel.Painel9505 ; run ;





/* We impute wage values if the person has a missing wage but stayed in the same job during 2 or 3 consecutive years */
/*(drop = REM_MEDIA REM_DEZ TEMP_EMPR wage lead_cnpj8 lag_cnpj8)*/
data painel.painel_final  ;

	set painel.painel_final ;


	wage_corr = wage ;

	flag = 0 ;

	/* Imputing Wage */
	if compress(CNPJ8) = compress(lead_cnpj8) and compress(CNPJ8) ~= "" and 
       (lead_wage ~= 0 and lead_wage ~= .) and (wage = . or wage = 0) then 
	do ;
		wage_corr = lead_wage ;
		flag = 1 ;
	end ;
	if compress(CNPJ8) = compress(lag_cnpj8) and compress(CNPJ8) ~= "" and 
       (lag_wage ~= 0 and lag_wage ~= .) and (wage = . or wage = 0) then 
	do ;
		wage_corr = lag_wage ;
		flag = 2 ;
	end ;
	if compress(CNPJ8) = compress(lead_cnpj8) and compress(CNPJ8) = compress(lag_cnpj8) and 
	   compress(CNPJ8) ~= "" and (lag_wage ~= 0 and lag_wage ~= .) and (lead_wage ~= 0 and lead_wage ~= .) and
       (wage = . or wage = 0) then 
	do; 
		wage_corr = (lag_wage + lead_wage)/2 ;
		flag = 3 ;
	end ;
	idade95 = idade + (1995 - ano) ;



	/*Imputing Sector and Wage*/
	if compress(lag_CNPJ8) = compress(lead_cnpj8) and compress(lag_CNPJ8) ~= "" and 
       lag_mes_deslig = '13' and lead_mes_adm in ('0','00') and (compress(CNPJ8) = "" and wage = .) then 
	do ;
		setor = lag_setor ;
		cnpj8 = lag_cnpj8 ;
		wage_corr = (lag_wage + lead_wage) / 2 ;		
		temp = 1 ;
		flag = 4 ;
	end ;

	
	if compress(lag_CNPJ8) = compress(lead2_cnpj8) and lag_mes_deslig = '13' and lead2_mes_adm in ('00','0') and  
	   (compress(CNPJ8) = "" and wage = .) then
	do ;
		setor = lag_setor ;
		cnpj8 = lag_cnpj8 ;
		wage_corr = (lag_wage + lead2_wage) / 2 ;		
		temp = 1 ;
		flag = 5 ;
	end ;


	if compress(lag2_CNPJ8) = compress(lead_cnpj8) and lag2_mes_deslig = '13' and lead_mes_adm in ('00','0') and  
	   (compress(CNPJ8) = "" and wage = .) then
	do ;
		setor = lag2_setor ;
		cnpj8 = lag2_cnpj8 ;
		wage_corr = (lag2_wage + lead_wage) / 2 ;		
		temp = 1 ;
		flag = 6 ;
	end ;



	cat_educ = . ;
	if grau_instr = 1 or grau_instr = 2 or grau_instr = 3 then cat_educ = 1 ;
	if grau_instr = 4 or grau_instr = 5 or grau_instr = 6 then cat_educ = 2 ;
	if grau_instr = 7                                     then cat_educ = 3 ;
	if grau_instr = 8 or grau_instr = 9                   then cat_educ = 4 ;
	if grau_instr = . or grau_instr <= 0                  then cat_educ = . ;

	if temp = 0 then wage_corr = . ;
	if wage_corr <= 0 then wage_corr = . ;


run ;



data inpc_min_wage (drop = year) ;
	set painel.inpc_min_wage ;

	ano = year ;

run ;


proc sort data = painel.painel_final sortsize = max ; by ano ; run ;

proc sort data = inpc_min_wage sortsize = max ; by ano ; run ;

data painel.painel_final ;
	merge painel.painel_final (in = a) inpc_min_wage ;
	by ano ;
	if a ;
	rwage_corr = wage_corr*min_wage / inpc05 ;
run ;



proc delete data = inpc_min_wage ; run ;



proc sort data = painel.painel_final ; by pis ano ; run ;



data painel.painel_final (keep = pis ano rwage_corr horas_contr idade idade95 sexo 
                           cat_educ SUBS_IBGE CODEMUN Setor CNPJ cnpj8 CLAS_CNAE adm deslig temp) ;
	set painel.painel_final ;
run ;










/*Getting the mode of the "age distribution" for each individual */

data painel_idade ;
	set painel.painel_final ;
	if idade95 = . then delete ;
	if ano >= 1994 ;
run ;


proc sort data = painel_idade sortsize = max ; by pis idade95 ; run ;


proc means noprint data = painel_idade sumsize = max ;
	by pis idade95 ;
	var idade95 ;
	output out = mode_idade (keep = pis idade95 freq_idade)
	n(idade95) = freq_idade ;
run ;


proc delete data = painel_idade ; run ;


/* We can have 2 age modes. If that is the case, we will trust the latest gender report */

proc sort data = mode_idade sortsize = max ; by pis descending freq_idade descending idade95 ; run ;


data mode_idade1 (keep = pis idade_mode1) ;
	set mode_idade ;
	by pis ;
	if first.pis ;
	rename idade95 = idade_mode1 ;
run ;

proc sort data = mode_idade sortsize = max ; by pis descending freq_idade idade95 ; run ;

data mode_idade2 (keep = pis idade_mode2) ;
	set mode_idade ;
	by pis ;
	if first.pis ;
	rename idade95 = idade_mode2 ;
run ;



proc sort data = mode_idade1 sortsize = max ; by pis ; run ;
proc sort data = mode_idade2 sortsize = max ; by pis ; run ;


proc sort data = painel.painel_final sortsize = max ; by pis ; run ;

data painel.painel_final ;
	merge painel.painel_final (in = a) mode_idade1 ;
	by pis ;
	if a ;
run ;


data painel.painel_final ;
	merge painel.painel_final (in = a) mode_idade2 ;
	by pis ;
	if a ;
	flag = 0 ;
	if idade_mode1 ~= idade_mode2 then flag = 1 ;
run ;


proc delete data = mode_idade1 ; run ;
proc delete data = mode_idade2 ; run ;
proc delete data = mode_idade ; run ;


data painel_idade (keep = pis ano idade) ;
	set painel.painel_final ;
	if idade = . then delete  ;
	if ano >= 1994 ;
run ;


proc sort data = painel_idade ; by pis descending ano ; run ;



data painel_idade (keep = pis idade_latest) ;
	set painel_idade ;
	by pis ;
	if first.pis ;
	idade_latest = idade + (1995-ano) ;
run ;




proc sort data = painel_idade ; by pis ; run ;

proc sort data = painel.painel_final ; by pis ano ; run ;



data painel.painel_final (drop = idade_mode1 idade_mode2 flag idade_latest) ;
	merge painel.painel_final (in = a) painel_idade ;
	by pis ;
	idade_new = idade_mode1 + (ano-1995) ;
	if flag = 1 then idade_new = idade_latest + (ano-1995) ;	
	if idade_new < 25 then cat_educ = . ;
run ;


proc delete data = painel_idade ; run ;









/*Getting the mode of the "education distribution" for each individual */

data painel_educ ;
	set painel.painel_final ;
	if ano >= 1994 ;
run ;


proc sort data = painel_educ sortsize = max ; by pis cat_educ ; run ;

proc means noprint data = painel_educ sumsize = max ;
	by pis cat_educ ;
	var cat_educ ;
	output out = mode_educ (keep = pis cat_educ freq_educ)
	n(cat_educ) = freq_educ ;
run ;


proc sort data = mode_educ sortsize = max ; by pis descending freq_educ descending cat_educ ; run ;


data mode_educ1 (keep = pis educ_mode1) ;
	set mode_educ ;
	by pis ;
	if first.pis ;
	rename cat_educ = educ_mode1 ;
run ;


proc sort data = mode_educ sortsize = max ; by pis descending freq_educ cat_educ ; run ;

data mode_educ2 (keep = pis educ_mode2) ;
	set mode_educ ;
	by pis ;
	if first.pis ;
	rename cat_educ = educ_mode2 ;
run ;


proc sort data = mode_educ1 sortsize = max ; by pis ; run ;

proc sort data = mode_educ2 sortsize = max ; by pis ; run ;



proc sort data = painel.painel_final sortsize = max ; by pis ; run ;

data painel.painel_final ;
	merge painel.painel_final (in = a) mode_educ1 ;
	by pis ;
	if a ;
run ;


data painel.painel_final ;
	merge painel.painel_final (in = a) mode_educ2 ;
	by pis ;
	if a ;
	flag = 0 ;
	if educ_mode1 ~= educ_mode2 then flag = 1 ;
run ;


proc delete data = mode_educ1 ; run ;
proc delete data = mode_educ2 ; run ;
proc delete data = mode_educ  ; run ;





data painel_educ (keep = pis ano cat_educ) ;
	set painel.painel_final ;
	if cat_educ = . then delete  ;
	if ano >= 1994 ;
run ;


proc sort data = painel_educ ; by pis descending ano ; run ;



data painel_educ (keep = pis educ_latest) ;
	set painel_educ ;
	by pis ;
	if first.pis ;
	rename cat_educ = educ_latest ;
run ;




proc sort data = painel_educ ; by pis ; run ;

proc sort data = painel.painel_final ; by pis ; run ;



data painel.painel_final (drop = educ_mode1 educ_mode2 flag educ_latest) ;
	merge painel.painel_final (in = a) painel_educ ;
	by pis ;
	educ_new = educ_mode1 ;
	if flag = 1 then educ_new = educ_latest ;
run ;



proc delete data = painel_educ  ; run ;
















/****************************************************************************/
/****************************************************************************/

/*Getting the mode of the "gender distribution" for each individual */


data painel_sexo ;
	set painel.painel_final ;
	if ano >= 1994 ;
run ;


proc sort data = painel_sexo sortsize = max ; by pis sexo ; run ;

proc means noprint data = painel_sexo sumsize = max ;
	by pis sexo ;
	var sexo ;
	output out = mode_sexo (keep = pis sexo freq_sexo)
	n(sexo) = freq_sexo ;
run ;



/* We can have 2 gender modes. If that is the case, we will trust the latest gender report */

proc sort data = mode_sexo sortsize = max ; by pis descending freq_sexo descending sexo ; run ;


data mode_sexo1 (keep = pis sexo_mode1) ;
	set mode_sexo ;
	by pis ;
	if first.pis ;
	rename sexo = sexo_mode1 ;
run ;


proc sort data = mode_sexo sortsize = max ; by pis descending freq_sexo sexo ; run ;


data mode_sexo2 (keep = pis sexo_mode2) ;
	set mode_sexo ;
	by pis ;
	if first.pis ;
	rename sexo = sexo_mode2 ;
run ;




proc sort data = mode_sexo1 sortsize = max ; by pis ; run ;

proc sort data = mode_sexo2 sortsize = max ; by pis ; run ;






proc sort data = painel.painel_final sortsize = max ; by pis ; run ;

data painel.painel_final ;
	merge painel.painel_final (in = a) mode_sexo1 ;
	by pis ;
	if a ;
run ;

data painel.painel_final ;
	merge painel.painel_final (in = a) mode_sexo2 ;
	by pis ;
	if a ;
	flag = 0 ;
	if sexo_mode1 ~= sexo_mode2 then flag = 1 ;
run ;


proc delete data = mode_sexo1 ; run ;
proc delete data = mode_sexo2 ; run ;
proc delete data = mode_sexo ; run ;




data painel_sexo (keep = pis ano sexo) ;
	set painel.painel_final ;
	if sexo = . then delete  ;
	if ano >= 1994 ;
run ;


proc sort data = painel_sexo ; by pis descending ano ; run ;


data painel_sexo (keep = pis sexo_latest) ;
	set painel_sexo ;
	by pis ;
	if first.pis ;
	rename sexo = sexo_latest ;
run ;


proc sort data = painel_sexo ; by pis ; run ;

proc sort data = painel.painel_final ; by pis ; run ;



data painel.painel_final (drop = sexo_mode1 sexo_mode2 sexo_latest flag) ;
	merge painel.painel_final (in = a) painel_sexo ;
	by pis ;
	sexo_new = sexo_mode1 ;
	if flag = 1 then sexo_new = sexo_latest ;
run ;


proc delete data = painel_sexo ; run ;




/****************************************************************************/
/****************************************************************************/











data painel.painel_final (keep = pis ano CNPJ cnpj8 codemun CLAS_CNAE SUBS_IBGE idade95 cat_educ sexo idade_new sexo_new
                           educ_new sector rwage_corr horas_contr adm deslig temp) ;
	set painel.painel_final ;

	sector = 0 ;
	if compress(setor) = 'Primary' then sector = 1 ;
	if compress(setor) = 'Low'     then sector = 2 ;
	if compress(setor) = 'High'    then sector = 3 ;
	if compress(setor) = 'Const'   then sector = 4 ;
	if compress(setor) = 'Trade'   then sector = 5 ;
	if compress(setor) = 'Trans'   then sector = 6 ;
	if compress(setor) = 'Service' then sector = 7 ;
	if temp            = 0         then sector = 0 ;

run ;


proc sort data = painel.painel_final sortsize = max ; by pis ano ; run ;







/*If a worker is in the same establishment for two consecutive years but appears with different sectors,*/
/*Then I correct for his sector in the second year as being the same as in the first. That is, if a */
/*workers did not change employment, he did not change sector either */

data painel_aux_lag (keep = pis ano lag_cnpj lag_sector) ;
	set painel.painel_final ;
	ano = ano + 1 ;
	lag_sector = sector ;
	lag_cnpj = CNPJ ;
run ;


data painel.painel_final ;
	merge painel.painel_final (in = a) painel_aux_lag ;
	by pis ano ;
	if a ;
run ;

data painel ;
	set painel.painel_final ;
	flag = 0 ;
	if compress(lag_cnpj) = compress(CNPJ) and lag_sector ~= sector
    and sector ~= 0 and sector ~= . then flag = 1 ;
	if compress(lag_cnpj) = compress(CNPJ) and lag_sector ~= sector
    and sector ~= 0 and sector ~= . then sector = lag_sector ;
run ;


proc freq data = painel ; 
	tables flag / out = flag ; 
run;


%macro correct_sector ;

data painel_aux_lag (keep = pis ano lag_cnpj lag_sector) ;
	set painel ;
	ano = ano + 1 ;
	lag_sector = sector ;
	lag_cnpj = CNPJ ;
run ;


data painel ;
	merge painel (in = a) painel_aux_lag ;
	by pis ano ;
	if a ;
run ;

data painel ;
	set painel ;
	flag = 0 ;
	if compress(lag_cnpj) = compress(CNPJ) and lag_sector ~= sector
    and sector ~= 0 and sector ~= . then flag = 1 ;
	if compress(lag_cnpj) = compress(CNPJ) and lag_sector ~= sector
    and sector ~= 0 and sector ~= . then sector = lag_sector ;
run ;

proc freq data = painel ; 
	tables flag / out = flag ; 
run;

%mend ;

%correct_sector ;
%correct_sector ;
%correct_sector ;
%correct_sector ;
%correct_sector ;
%correct_sector ;
%correct_sector ;
%correct_sector ;
%correct_sector ;
%correct_sector ;
%correct_sector ;
%correct_sector ;
%correct_sector ;
%correct_sector ;
%correct_sector ;
%correct_sector ;
%correct_sector ;


data painel.painel_final (keep = pis ano CNPJ cnpj8 codemun CLAS_CNAE SUBS_IBGE idade95 cat_educ sexo idade_new sexo_new
                           educ_new sector rwage_corr horas_contr adm deslig temp) ;
	set painel ;
run ;



proc delete data = painel ; run ;
proc delete data = painel_aux_lag ; run ;
proc delete data = flag ; run ;
