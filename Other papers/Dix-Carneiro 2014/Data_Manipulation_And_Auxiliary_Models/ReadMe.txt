For questions regarding these codes, please contact
Rafael Dix-Carneiro at rafael.dix.carneiro@duke.edu or 
dix-carneiro@econ.umd.edu

This folder contains the SAS codes used in order to construct the panel
of workers used in the paper. This panel of workers in then used in order
to estimate the auxiliary models used in the estimation procedure.

The following codes are executed, in this order:
1_amostra_pis.sas
2_select_pis.sas
3_CNPJ_CNAE05.sas
4_CNPJ_CNAE.sas
5_Amostra_RAIS_Final.sas
6_Painel_Final.sas
7_Regressions.sas
8_Initial_Conditions.sas
9_WageData_Expfactor.sas

The output of these codes is in folder Estimation_Data.

====
NOTE
====

The raw microdata employed in the paper comes from records of the
Relacao Anual de Informacoes Sociais (RAIS), collected from the
Brazilian Ministry of Labor. This paper uses rounds of the data 
from 1986 to 2005.

RAIS is not publicly available. I had access to it through Instituto
de Pesquisa Economica Aplicada (IPEA) in Brasilia. I am grateful to
Joao De Negri for granting me access. He can be reached at 
joao.denegri@ipea.gov.br

The estimation codes can be run using the results of the auxiliary models,
as well as the empirical distribution of initial conditions, which are 
made available in the Estimation_Data folder.
