%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This routine does two things:
%
% 1. Reads in all calibrated parameters and time series from Stata
% 2. Calculates Implied Productivity Series
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

outPathMatlab = '~/Dropbox/DanishLEEDProject/structuralEstimation_Final/Calibration/matlabFiles/';
outPathStata = '~/Dropbox/DanishLEEDProject/structuralEstimation_Final/Estimation/outputs/';
sigma = 4;
%sigma = 3;

ALPHA = csvread(strcat(outPathMatlab,'ALPHA_D.csv'),0,0);

BK_D = csvread(strcat(outPathMatlab,'BK_D.csv'),0,0);
BL_D = csvread(strcat(outPathMatlab,'BL_D.csv'),0,0);

BTT_D = csvread(strcat(outPathMatlab,'BTT_D.csv'),0,0);
BTN_D = csvread(strcat(outPathMatlab,'BTN_D.csv'),0,0);
BNT_D = csvread(strcat(outPathMatlab,'BNT_D.csv'),0,0);
BNN_D = csvread(strcat(outPathMatlab,'BNN_D.csv'),0,0);

RELPRICES_D = csvread(strcat(outPathMatlab,'RELPRICES_D.csv'),0,0);

DELPF_D = csvread(strcat(outPathMatlab,'DELPF_D.csv'),0,0);
DELPF_D(:,3) = log(DELPF_D(:,3));
DELNTP_D = csvread(strcat(outPathMatlab,'DELNTP_D.csv'),0,0);
DELNTP_D(:,3) = log(DELNTP_D(:,3));
DELR_D = csvread(strcat(outPathMatlab,'DELR_D.csv'),0,0);

CAPINC_D = csvread(strcat(outPathMatlab,'CAPINCOME_D.csv'),0,0);
EXPSHIFT_D = csvread(strcat(outPathMatlab,'EXPSHIFT_D.csv'),0,0);

DELRELDEMAND_D = csvread(strcat(outPathMatlab,'RELDEMAND_D.csv'),0,0);

BK = BK_D(:,2);
BL = BL_D(:,2:end);
BTT = BTT_D(:,2:end);
BTN = BTN_D(:,2:end);
BNT = BNT_D(:,2:end);
BNN = BNN_D(:,2:end);

BM = [BTT BTN; BNT BNN];

W = csvread(strcat(outPathMatlab,'WAGERATES.csv'),0,0);


COEFFICIENTS = [BK BL BM];

T = max(DELNTP_D(:,1));
NINDS = max(DELNTP_D(:,2));
NTRADABLES = max(DELPF_D(:,2));


DELZF = zeros(NINDS,T);

for t = 1:T
    r0 = RELPRICES_D(RELPRICES_D(:,1)==t-1,3)';
    delW = log(W(W(:,1)==t,3)) - log(W(W(:,1)==t-1,3));
    
    delRELP = log(RELPRICES_D(RELPRICES_D(:,1)==t,3)) - log(RELPRICES_D(RELPRICES_D(:,1)==t-1,3));
    delPK = DELR_D(DELR_D(:,1)==t,2);
    delPF1 = DELPF_D(DELPF_D(:,1)==t,3);
    delNT1 = DELNTP_D(DELNTP_D(:,1)==t,3);
    
    DELZF(:,t) = [BTT - eye(NTRADABLES); BNT]*delPF1 + BL*delW +  BK*delPK + ([bsxfun(@times,[BTT;BNT],r0.^(1-sigma)./(1+r0.^(1-sigma))) [BTN;BNN]] - eye(NINDS))*[delRELP;delNT1];
end

clear delRELP delPK delPF1 delNT1 delW r0