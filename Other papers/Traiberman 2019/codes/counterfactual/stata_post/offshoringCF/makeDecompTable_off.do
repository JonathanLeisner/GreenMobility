global doFiles "/Users/sharont/Dropbox/DanishLEEDProject/structuralEstimation_Final/Estimation/doFiles"
global dataPath "/Users/sharont/Dropbox/DanishLEEDProject/structuralEstimation_Final/Estimation/inputs"
global outPath "/Users/sharont/Dropbox/DanishLEEDProject/structuralEstimation_Final/Estimation/outputs"
global matlabPath "/Users/sharont/Dropbox/DanishLEEDProject/structuralEstimation_Final/Estimation/matlabFiles"

global figurePath "/Users/sharont/Dropbox/DanishLEEDProject/structuralEstimation_Final/Estimation/figures"
global tablePath "/Users/sharont/Dropbox/DanishLEEDProject/structuralEstimation_Final/Estimation/tables"
global workerFiles /Users/sharont/Dropbox/DanishLEEDProject/structuralEstimation_Final/Calibration/workerFiles

global timingInputs ~/Dropbox/DanishLEEDProject/structuralEstimation_Final/counterfactualOutputs/timingCF
global offshoringInputs ~/Dropbox/DanishLEEDProject/structuralEstimation_Final/counterfactualOutputs/offshoringCF

insheet using $offshoringInputs/varDecomp_OFF.csv, clear c
format v* %4.1f

gen amp = "&"
gen lb = "\\"
gen id = _n

gen rowName = "Occs/Sectors" if id==1
replace rowName = "Occupations Only" if id==2
replace rowName = "Sectors Only" if id==3

order rowName v*

outsheet rowName amp v1 amp v2 amp v3 lb using $tablePath/cfIncomeVarDecomp_OFF.tex, noq replace non 
