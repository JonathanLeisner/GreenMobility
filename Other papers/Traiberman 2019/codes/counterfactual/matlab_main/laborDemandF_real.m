function [HD, E] = laborDemandF_real(techParams,r1,INCOME,EXPSHIFT,w)

% Counts
NINDS = techParams.NINDS;
NTRADABLES = techParams.NTRADABLES;

% IO Matrix
BTT = techParams.BTT;
BTN = techParams.BTN;
BNT = techParams.BNT;
BNN = techParams.BNN;

BL = techParams.BL;
BK = techParams.BK;

% Consumer Parameters
ALPHA = techParams.ALPHA;
sigma = techParams.sigma;

% 
% % Step 1: Calculate Change in Prices from Log-Linearized System
% MIDMATBM0 = [bsxfun(@times,[BTT;BNT],r0.^(1-sigma)./(1+r0.^(1-sigma))) [BTN;BNN]];
% 
% X = (eye(NINDS) - MIDMATBM0)^(-1) * (-delZ + [BTT - eye(NTRADABLES); BNT]*delPF + BK*delPK + BL*delW);
% 
% r1 = r0.*exp(X(1:NTRADABLES))';


% Step 2: Calculate Expenditure
MIDMATBM1 = [bsxfun(@times,[BTT;BNT],r1.^(1-sigma)./(1+r1.^(1-sigma))) [BTN;BNN]];

MIDMATAL1 = [ALPHA(1:NTRADABLES,2).*r1'.^(1-sigma)./(1+r1'.^(1-sigma));
            ALPHA(NTRADABLES+1:end,2)];
           
MIDMATEX1 = [EXPSHIFT.*r1'.^(1-sigma);
            zeros(NINDS-NTRADABLES,1)];


E = (eye(NINDS) - MIDMATBM1')\(MIDMATAL1*(INCOME) + MIDMATEX1);

% Step 3: Calculate Labor Demand
HD = (BL'*E)./exp(w);