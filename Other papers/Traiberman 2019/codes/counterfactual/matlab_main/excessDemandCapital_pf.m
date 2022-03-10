function [out, HD, HS] = excessDemandCapital_pf(w1,pK1,w0,pK0,V,VU,P0,L0,LU0,n,un,b,params,techParams,r0,delPF,delZ,KINC,EXPSHIFT,TA,TA_U,cMat,hMat,state)

% Step 0: Unpack parameters
% Demand Parameters
sigma = techParams.sigma;
ALPHA = techParams.ALPHA(:,2);

% Programming Parameters
NTRADABLES = techParams.NTRADABLES;
NINDS = techParams.NINDS;

% IO Matrix
BTT = techParams.BTT;
BTN = techParams.BTN;
BNT = techParams.BNT;
BNN = techParams.BNN;

BL = techParams.BL;
BK = techParams.BK;


% Step 1: Given a guess of wages and exogeneous stuff, update the price level.
delW = w1 - w0;
delPK = log(pK1) - log(pK0);

MIDMATBM0 = [bsxfun(@times,[BTT;BNT],r0.^(1-sigma)./(1+r0.^(1-sigma))) [BTN;BNN]];

delPrices = (eye(NINDS) - MIDMATBM0)^(-1) * (-delZ + [BTT - eye(NTRADABLES); BNT]*delPF + BK*delPK + BL*delW);
r1 = r0.*exp(delPrices(1:NTRADABLES))';

adjTerm = 1/(1-sigma)*(log(1+r1.^(1-sigma)) - log(1+r0.^(1-sigma)));
delP = delPrices;
delP(1:NTRADABLES) = delPF + adjTerm';
P1 = P0*exp(dot(ALPHA,delP));





% Step 2: Use REAL wages to update labor supply
[HS, ~, ~, ~] = laborSupplyF_pf_new(w1-log(P1),V,VU,L0,LU0,n,un,b,params,TA,TA_U,cMat,hMat,state);



% Step 3: Get Change in NOMINAL Wages and Income
WINC = dot(HS,exp(w1));
KS = KINC*pK1;
INCOME = WINC + KS;

% Exports are NUMERAIRE (no adjustment)
EXPTERM = EXPSHIFT;

% Step 4: Solve Firm's Problem for Labor Demand
[HD, E] = laborDemandF_real(techParams,r1,INCOME,EXPTERM,w1);
KD = BK'*E;


out = [(HD - HS)./HD; (KD - KS)./KS];

