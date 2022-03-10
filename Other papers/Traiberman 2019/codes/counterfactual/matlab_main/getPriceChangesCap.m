function [Ppath] = getPriceChangesCap(wPath,pKPath,techParams,P0,r0,delPFpath,delZpath)

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


numT = size(wPath,2);
Ppath = zeros(1,numT);
Ppath(1) = P0;

for period = 2:numT
    
    % Change in nominal wages
    delW = wPath(:,period) - wPath(:,period-1);
    delPK = log(pKPath(period)) - log(pKPath(period-1));

    % Peel off current period changes in exogenous prices
    delPF = delPFpath(:,period-1);
    delZ = delZpath(:,period-1);
    
    
    MIDMATBM0 = [bsxfun(@times,[BTT;BNT],r0.^(1-sigma)./(1+r0.^(1-sigma))) [BTN;BNN]];
    
    % Change in relative prices of foreign goods + nontradable prices
    delPrices = (eye(NINDS) - MIDMATBM0)^(-1) * (-delZ + [BTT - eye(NTRADABLES); BNT]*delPF + BK*delPK + BL*delW);
    
    delP = delPrices;
    
    % Relative prices in levels now
    r1 = r0.*exp(delPrices(1:NTRADABLES))';

    % Term to adjust changes in tradables 
    adjTerm = 1/(1-sigma)*(log(1+r1.^(1-sigma)) - log(1+r0.^(1-sigma)));
    delP(1:NTRADABLES) = delPF + adjTerm';
    
    % Change in consumer price index
    Ppath(period) = P0*exp(dot(ALPHA,delP));
    P0 = Ppath(period);
    r0 = r1;
end
