% Paths for Estimated Labor Supply
%estOutPath = '/Users/sharontraiberman/Dropbox/DanishLEEDProject/structuralEstimation_Final/Estimation/matlabFiles/';
%stataOutPath = '/Users/sharontraiberman/Dropbox/DanishLEEDProject/structuralEstimation_Final/Estimation/outputs/';
%stataInPath = '/Users/sharontraiberman/Dropbox/DanishLEEDProject/structuralEstimation_Final/Estimation/inputs/';

% Paths for Calibrated Labor Demand
%outPathMatlab = '/Users/sharontraiberman/Dropbox/DanishLEEDProject/structuralEstimation_Final/Calibration/matlabFiles/';


% Paths for Estimated Labor Supply
estOutPath = '~/Dropbox/DanishLEEDProject/structuralEstimation_Final/Estimation/matlabFiles/';
stataOutPath = '~/Dropbox/DanishLEEDProject/structuralEstimation_Final/Estimation/outputs/';
stataInPath = '~/Dropbox/DanishLEEDProject/structuralEstimation_Final/Estimation/inputs/';

% Paths for Calibrated Labor Demand
outPathMatlab = '~/Dropbox/DanishLEEDProject/structuralEstimation_Final/Calibration/matlabFiles/';


GAMMA = csvread(strcat(estOutPath,'GAMMAQ.csv'));

ETA  = csvread(strcat(estOutPath,'ETAQ.csv'));

SEARCH = csvread(strcat(estOutPath,'SEARCHQ.csv'));

RHO = csvread(strcat(estOutPath,'RHOQ.csv'));

NONEMP = csvread(strcat(estOutPath,'NONEMPQ.csv'));

MINCER = csvread(strcat(stataOutPath,'MINCER.csv'));

INITIALEDF = csvread(strcat(outPathMatlab,'workerEDF1996.csv'));

INITIALEDF_U = csvread(strcat(outPathMatlab,'UworkerEDF1996.csv'));

%distMatrix = csvread('occDistMatrix.csv');
%distMatrix = csvread(strcat(stataInPath,'occDistMatrixS.csv'));

charMat = csvread(strcat(stataOutPath,'occDistMat.csv'),0,0);
distMatrix = charMat(:,2:end);
