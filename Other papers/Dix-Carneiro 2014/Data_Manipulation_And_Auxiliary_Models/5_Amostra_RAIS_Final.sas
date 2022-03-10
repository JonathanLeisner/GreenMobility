/*This code is written in order to assign a sector to each observation.*/
/*This is done using the CNPJ_CNAE9498 in order to impute and/or correct for */
/*classification from 1986 to 1994 and CNPJ_9405 in order to impute a classification*/
/*in case and individual is working in a firm, getting a wage but the classification for */
/*his firm is missing.*/
/*The data file Painel9505 is generated*/
/*This is the fifth code to be run in order to construct the panel of workers.*/




libname painel "C:\Users\Dix-Carneiro\Documents\Trade and Labor Dynamics\Data for Estimation New\Codes_7_Sectors\PanelRAIS" ;


data painel.Painel9505 ;
	set painel.am_rais_final ;
	wage = rem_media ;
	if (wage = . or wage = 0) and (rem_dez ~= . and rem_dez ~= 0) then wage = rem_dez ;
run ;



/* Corrige CNPJ em 1990 */
/* Em 1990 varios CNPJ's aparecem com formato errado */

data painel.Painel9505  (drop = CNPJ_new flag);
	length CNPJ8 $ 8 ;
	set painel.Painel9505 ;

	CNPJ_new = compress(CNPJ,".") ;

	flag = 0 ;
	
	if compare(CNPJ,CNPJ_new,':') ~= 0 and ano = 1990 then flag = 1 ;

	CNPJ8 = substr(CNPJ,1,8) ;
	if flag = 1 then CNPJ8 = substr(CNPJ_new,3,10) ;

run ;


/*Para anos antes de 1994 utilizamos informacao sobre a CNAE reportada                 */
/*entre 1994 e 2005 de forma a obter mais informacao sobre a atividade economica       */


data Cnpj_cnae (drop = quant) ;
	length CNPJ8 $ 8 ;
	set painel.Cnpj_cnae9405 ;
	rename CNAE = CNAE_post94 ;
run ;


proc sort data = cnpj_cnae sortsize = max ; by CNPJ8 ; run ;


proc sort data = painel.Painel9505 sortsize = max ; by CNPJ8 ; run ;


data painel.Painel9505 ;
	length CLAS_CNAE $ 5 ;
	length CNAE2 $ 2 ;
	length CNAE_post94 $ 2 ;
	merge painel.Painel9505 (in = a) cnpj_cnae ;
	by cnpj8 ;

	if a ;
	CNAE2 = substr(CLAS_CNAE,1,2) ;    
	if sexo = 3 then sexo = 2 ;

run ;

proc sort data = painel.Painel9505 ; by cnpj8 ; run ;



/* Corrige alguns SUBS_IBGE - Nao dah para corrigir todos  */
/* Pois se SUBS_IBGE = IND, nao sabemos separar entre alta tecnologia ou baixa tecnologia */
/* Utilizamos informacao pos-1994 de forma a classificar as empresas que declararam IND como subsetor IBGE */

data painel.Painel9505 ;
	length Setor $ 10 ;
	set painel.Painel9505 ;

	Setor = 'Missing' ;

	/* 1993 e antes */
	/* Utilizando apenas informacao do SUBS_IBGE */

	if compress(SUBS_IBGE) = 'EXTR' and ano <= 1993 then Setor = 'Primary' ;
	if compress(SUBS_IBGE) in ('PAPE','BOR') and ano <= 1993 then Setor = 'Low' ;
	if compress(SUBS_IBGE) = 'IND' and ano <= 1993 then Setor = 'Missing' ;
	if compress(SUBS_IBGE) in ('ADM','ALOJ','MED') and ano <= 1993 then Setor = 'Service' ;
	if compress(SUBS_IBGE) in ('CONS') and ano <= 1993 then Setor = 'Const' ;
	if compress(SUBS_IBGE) in ('TRAN','SER') and ano <= 1993 then Setor = 'Trans' ;


	if compress(SUBS_IBGE) in ('1101','4405') and ano <= 1993 then Setor = 'Primary' ;
	if compress(SUBS_IBGE) in ('4506','4507','4509','4511','4513','4514','4516','4517') and ano <= 1993 then Setor = 'Low' ;
	if compress(SUBS_IBGE) in ('4508','4510','4512','4515') and ano <= 1993 then Setor = 'High' ;
	if compress(SUBS_IBGE) in ('5719','5820','5821','5822','5823','5824') and ano <= 1993 then Setor = 'Service' ;
	if compress(SUBS_IBGE) in ('3304') and ano <= 1993 then Setor = 'Const' ;
	if compress(SUBS_IBGE) in ('5825','4618') and ano <= 1993 then Setor = 'Trans' ;
	if compress(SUBS_IBGE) in ('2202','2203') and ano <= 1993 then Setor = 'Trade' ;


	/* Utilizando informacao sobre a CNAE pos 1994 para imputar o setor em 1993 e antes */

	if compress(CNAE_post94) in ('23','33') and ano <= 1993 then Setor = 'High' ;
	if compress(CNAE_post94) in ('25') and ano <= 1993 then Setor = 'Low' ;

	if compress(SUBS_IBGE) = 'IND' and ano <= 1993 and 
   	compress(CNAE_post94) in ('23','24','29','30','31','32','33','34','35') then Setor = 'High' ;

	if compress(SUBS_IBGE) = 'IND' and ano <= 1993 and 
   	compress(CNAE_post94) in ('15','16','17','18','19','20','21','22','25','26','27','28','36','37') then Setor = 'Low' ;


	/*	Ainda teremos varias observacoes com setor "Missing" pois foram declarados*/
	/*	com SUBS_IBGE = 9999 ou sem qualquer numero em SUBS_IBGE*/
	/*	Tentamos imputar esses valores a seguir, utilizando a CNAE mais frequentemente */
	/*	reportada entre 1994 e 1998. Obviamente, esse procedimento nao consegue consertar*/
	/*	todos os missing para setor antes de 1993 pois muitas empresas presentes entre 1986 e 1993*/
	/*	podem ter saido da amostra*/


	/* A partir de 1994 a selecao eh feita a partir da CNAE2 */
	/* Ao contrario da OECD, produtos derivados do petroleo, 23, sao considerados */
	/* Alta tecnologia */

	if compress(CNAE2) in ('01','02','05','10','11','13','14') and ano >= 1994 then Setor = 'Primary' ;
	if compress(CNAE2) in ('15','16','17','18','19','20','21','22','25','26','27','28','36','37') 
                   and ano >= 1994 then Setor = 'Low' ;
	if compress(CNAE2) in ('23','24','29','30','31','32','33','34','35') 
                   and ano >= 1994 then Setor = 'High' ;


	if compress(CNAE2) in ('45') and ano >= 1994 then Setor = 'Const' ;
	if compress(CNAE2) in ('50','51','52') and ano >= 1994 then Setor = 'Trade' ;
	if compress(CNAE2) in ('40','41','60','61','62','63','64','90') and ano >= 1994 then Setor = 'Trans' ;
	if compress(CNAE2) in ('55','65','66','67','70','71','72','73','74',
                           '75','80','85','91','92','93','95','99') and ano >= 1994 then Setor = 'Service' ;



	/* Para imputar valores Missing ateh 1993 utiliza-se a variavel CNAE_post94, que eh o modo entre 1994 e 2005 */
	/* Para imputar valores Missing a partir de 1994 utiliza-se a variavel CNAE_post94, que eh o modo entre 1994 e 2005 */

	if Setor = 'Missing' and ano <= 1993 and compress(CNAE_post94) 
	in ('01','02','05','10','11','13','14') then Setor = 'Primary' ;
	if Setor = 'Missing' and ano <= 1993 and compress(CNAE_post94) 
	in ('15','16','17','18','19','20','21','22','25','26','27','28','36','37') then Setor = 'Low' ;
	if Setor = 'Missing' and ano <= 1993 and compress(CNAE_post94) 
	in ('23','24','29','30','31','32','33','34','35') then Setor = 'High' ;

	if Setor = 'Missing' and ano <= 1993 and compress(CNAE_post94) 
	in ('45') then Setor = 'Const' ;
	if Setor = 'Missing' and ano <= 1993 and compress(CNAE_post94) 
	in ('50','51','52') then Setor = 'Trade' ;
	if Setor = 'Missing' and ano <= 1993 and compress(CNAE_post94) 
	in ('40','41','60','61','62','63','64','90') then Setor = 'Trans' ;
	if Setor = 'Missing' and ano <= 1993 and compress(CNAE_post94) 
	in ('55','65','66','67','70','71','72','73','74',
                           '75','80','85','91','92','93','95','99') then Setor = 'Service' ;

run ;


proc delete data = painel.Cnpj_cnae9405 painel.CNPJ_CNAE9498 painel.am_rais_final ; run ;

